#ifndef TAPAS_APPROACH_H
#define TAPAS_APPROACH_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <string>
#include <random>
#include <map>
#include "../common/TrafficAssignmentApproach.h"
#include "./components/TapasShiftMethod.h"
#include "../../core/Network.h"

namespace TrafficAssignment {
  template <typename T>
  class TapasApproach : public TrafficAssignmentApproach <T> {
  public:
    TapasApproach(Network<T>& network, T alpha = 1e-6):
      TrafficAssignmentApproach<T>::TrafficAssignmentApproach(network, alpha) {

      shift_method_ = nullptr;

      for (int i = 0; i < this->network_.number_of_zones(); i++) {
        if (this->network_.origin_info()[i].size() > 0) {
          origins_.push_back(i);
        }
      }
      number_of_origins_ = origins_.size();
      bush_links_flows_.resize(this->number_of_origins_);
      origin_link_corresponding_pas_.resize(this->number_of_origins_);
      origin_corresponding_pas_set_.resize(this->number_of_origins_);
      links_origin_least_cost_routes_tree_.resize(this->number_of_origins_);
      parent_origin_least_cost_routes_tree_.resize(this->number_of_origins_);
      origin_link_pair_pas_.resize(this->number_of_origins_);

      for (int origin_index = 0; origin_index < this->number_of_origins_; origin_index++) {
        bush_links_flows_[origin_index].resize(this->network_.number_of_links(), 0);
        origin_link_corresponding_pas_[origin_index].resize(this->network_.number_of_links(), -1);
        parent_origin_least_cost_routes_tree_[origin_index].resize(this->network_.number_of_nodes(), -1);
      }
    }

    ~TapasApproach() {
      delete shift_method_;
    }

    void ComputeTrafficFlows() {
      this->statistics_recorder_.StartRecording(this->GetApproachName());
      AllOrNothingAssignment();
      this->statistics_recorder_.RecordStatistics();
      T prev = 0, now = 0;
      for (int i = 0; i < 20; i++) {
        EquilibrationIteration(i + 1);
        prev = now;
        now = this->network_.ObjectiveFunction();
        //std::cout << prev - now << ' ' << reserved_pas_hash_.size() << '\n';
      }
      std::cout << 1 << '\n';
    }

  protected:
    int number_of_origins_;

    // stores all index of every origin 
    std::vector <int> origins_;

    // stores all destinations for every origin
    std::vector <std::vector <int>> destinations_for_origin_;

    // stores the amount of flow for each link separately for each bush
    std::vector <std::vector <T>> bush_links_flows_;

    // stores for each link for each origin individually corresponding pas
    std::vector <std::vector <int>> origin_link_corresponding_pas_;

    //vector <vector <int>> bush_links_pas_;
    // Stores paired alternate segments descriprions by storing links list
    //unordered_map <int, pair <vector <int>, vector <int>>> paired_alternate_segments_;
    // Stores all origins that will perform flow redistribution by using certain paired alternate segment
    // displays from pas hah function to origins array
    std::unordered_map <int, std::unordered_set <int>> pas_users_;

    // For every origin stores least cost routes tree by storing for every node it's parent
    std::vector <std::vector <int>> parent_origin_least_cost_routes_tree_;

    // For every origin stores links being used in least cost routes tree
    std::vector <std::unordered_set <int>> links_origin_least_cost_routes_tree_;

    // For every pas we are going to calculate a hash corresponding to it
    std::unordered_map <int, std::pair <std::vector <int>, std::vector <int>>> reserved_pas_hash_;

    // 
    std::unordered_map <int, T> pas_flow_shift_starting_point_;

    //
    std::vector <std::unordered_set <int>> origin_corresponding_pas_set_;

    //
    TapasShiftMethod <T>* shift_method_;

    // 
    const T computation_threshold_ = 1e-10;

    std::priority_queue <std::pair <T, int>> delta_pas_queue;

    int pas_processed_last_queue_update;

    std::vector <std::map <std::pair <int, int>, int>> origin_link_pair_pas_;

    void AllOrNothingAssignment() {
      for (int origin_index = 0; origin_index < this->number_of_origins_; origin_index++) {
        std::vector <std::pair <int, std::vector <int>>> shortest_routes = this->network_.ComputeSingleOriginBestRoutes(origins_[origin_index]);
        this->network().AddRoutes(shortest_routes);
        for (int route_index = 0; route_index < shortest_routes.size(); route_index++) {
          for (auto link_index : shortest_routes[route_index].second) {
            bush_links_flows_[origin_index][link_index] += this->network_.od_pairs()[shortest_routes[route_index].first].GetDemand();
          }
        }
      }
    }

    // TEMPORARY SOLUTION
    void PasOriginCheck(int pas_hash, int origin_index) {
      if (this->pas_users_[pas_hash].count(origin_index)) {
        return;
      }
      if (Link<T>::GetLinksDelay(this->network_.links(), this->reserved_pas_hash_[pas_hash].first) >
        Link<T>::GetLinksDelay(this->network_.links(), this->reserved_pas_hash_[pas_hash].second)) {
        for (auto now : this->reserved_pas_hash_[pas_hash].first) {
          if (this->bush_links_flows_[origin_index][now] == 0) {
            return;
          }
        }
      }
      else {
        for (auto now : this->reserved_pas_hash_[pas_hash].second) {
          if (this->bush_links_flows_[origin_index][now] == 0) {
            return;
          }
        }
      }
      this->pas_users_[pas_hash].insert(origin_index);
    }

    void AllPasOriginCheck() {
      for (auto& now : this->reserved_pas_hash_) {
        for (int origin_index = 0; origin_index < this->number_of_origins_; origin_index++) {
          PasOriginCheck(now.first, origin_index);
        }
      }
    }

    void RemoveCycleFlow(int origin_index, std::vector <int>& cycle) {
      T mn_flow = bush_links_flows_[origin_index][cycle[0]];
      for (auto link_id : cycle) {
        mn_flow = std::min(mn_flow, bush_links_flows_[origin_index][link_id]);
      }
      for (auto link_id : cycle) {
        bush_links_flows_[origin_index][link_id] -= mn_flow;
        this->network_.SetLinkFlow(link_id, -mn_flow);
      }
    }

    void DfsForCycleIdentification(int origin_index, int cur_node, std::vector <int>& link_used, std::vector <int>& state, bool& cycle_found) {
      if (cycle_found) {
        return;
      }
      if (state[cur_node] != 0) {
        return;
      }
      state[cur_node] = 1;
      for (const auto now : this->network_.adjacency()[cur_node]) {
        if (cycle_found) {
          return;
        }
        int next_node = this->network_.links()[now].term;
        //int next_node = this->links_[now].term;
        if ((bush_links_flows_[origin_index][now] < this->computation_threshold_) || (state[next_node] == 2)) {
          continue;
        }
        if (state[next_node] == 1) {
          cycle_found = true;
          std::vector <int> cycle;
          cycle.push_back(now);
          int temp = cur_node;
          while (temp != next_node) {
            cycle.push_back(link_used[temp]);
            temp = this->network_.links()[temp].init;
            //temp = this->links_[link_used[temp]].init;
          }
          RemoveCycleFlow(origin_index, cycle);
          return;
        }
        if (state[next_node] == 0) {
          link_used[next_node] = now;
          DfsForCycleIdentification(origin_index, next_node, link_used, state, cycle_found);
        }
      }
      state[cur_node] = 2;
    }

    // Perfoms dfs on links with positive flow
    // if finds backward link, that means that cycle is found
    // when cycle is found, it gets removed, then all process of finding a cycle starts again from scratch 
    void SingleOriginRemoveAllCyclicFlows(int origin_index) {
      std::vector <int> link_used(this->network_.number_of_nodes(), -1);
      std::vector <int> state(this->network_.number_of_nodes(), 0);
      bool cycle_found = false;
      DfsForCycleIdentification(origin_index, origins_[origin_index], link_used, state, cycle_found);
      while (cycle_found) {
        cycle_found = false;
        for (int i = 0; i < this->network_.number_of_nodes(); i++) {
          state[i] = 0;
        }
        DfsForCycleIdentification(origin_index, origins_[origin_index], link_used, state, cycle_found);
      }
    }

    void RemoveAllCyclicFlows() {
      for (int origin_index = 0; origin_index < this->number_of_origins_; origin_index++) {
        SingleOriginRemoveAllCyclicFlows(origin_index);
      }
    }

    // For single origin fills parent_origin_least_cost_routes_tree_ and links_origin_least_cost_routes_tree_
    void BuildSingleOriginLeastCostRoutesTree(int origin_index) {
      this->links_origin_least_cost_routes_tree_[origin_index].clear();
      int current_origin = this->origins_[origin_index];
      // Pairs in q contain delay from origin and node index  
      std::priority_queue <std::pair <T, int>, std::vector <std::pair <T, int>>, std::greater <std::pair <T, int>>> q;
      int u;
      int link_index;
      int origin = this->origins_[origin_index];
      T cur_delay = 0;
      q.push({ 0, -1 });
      std::unordered_set <int> processed;
      while (!q.empty()) {
        cur_delay = q.top().first;
        link_index = q.top().second;
        if (q.top().second != -1) {
          u = this->network_.links()[link_index].term;
          //u = this->links_[link_index].term;
        }
        else {
          u = origin;
        }
        q.pop();
        if (processed.count(u)) {
          continue;
        }
        if (u != current_origin) {
          this->parent_origin_least_cost_routes_tree_[origin_index][u] = this->network_.links()[link_index].init;
          //this->parent_origin_least_cost_routes_tree_[origin_index][u] = this->links_[link_index].init;
          this->links_origin_least_cost_routes_tree_[origin_index].insert(link_index);
        }

        processed.insert(u);
        for (const auto now : this->network_.adjacency()[u]) {
          if (!processed.count(this->network_.links()[now].term)) {
          //if (!processed.count(this->links_[now].term)) {
            q.push({ cur_delay + this->network_.links()[now].Delay(), now });
          }
        }
      }
    }

    void BuildLeastCostRoutesTrees() {
      for (int i = 0; i < this->number_of_origins_; i++) {
        this->BuildSingleOriginLeastCostRoutesTree(i);
      }
    }

    std::pair <std::vector <int>, std::vector <int>> PasReconstruction(int origin_index, int pas_init_node, int link_index,
      const std::unordered_set <int>& origin_term_route, const std::unordered_map <int, int>& used_link) {

      std::vector <int> pas_part_least_cost_routes_tree, pas_part_not_least_cost_routes_tree;
      // pas_part_not_least_cost_routes_tree building
      int temp = pas_init_node;
      while (temp != this->network_.links()[link_index].term) {
        pas_part_not_least_cost_routes_tree.push_back(used_link.at(temp));
        // ????
        temp = this->network_.links()[used_link.at(temp)].term;
      }
      // pas_part_least_cost_routes_tree building
      temp = this->network_.links()[link_index].term;
      int used_link_index = 0, prev_node = 0;
      while (temp != pas_init_node) {
        prev_node = this->parent_origin_least_cost_routes_tree_[origin_index][temp];
        for (const auto now : this->network_.adjacency()[prev_node]) {
          if (this->network_.links()[now].term == temp) {
            pas_part_least_cost_routes_tree.push_back(now);
            break;
          }
        }
        temp = prev_node;
      }
      std::reverse(pas_part_least_cost_routes_tree.begin(), pas_part_least_cost_routes_tree.end());
      return { pas_part_least_cost_routes_tree, pas_part_not_least_cost_routes_tree };
    }
    int HashCombine(int seed, int value) {
      seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      return seed;
    }
    int HashArray(const std::vector<int>& array) {
      std::hash<int> hasher;
      int result = 0;
      for (const int& element : array) {
        result = HashCombine(result, hasher(element));
      }
      return result;
    }
    bool CheckHashPresence(const std::pair <std::vector <int>, std::vector <int>>& pas) {
      std::vector <int> temp = pas.first;
      temp.insert(temp.end(), pas.second.begin(), pas.second.end());
      int hash = HashArray(temp);
      while (1) {
        if (!this->reserved_pas_hash_.count(hash)) {
          return false;
        }
        if (this->reserved_pas_hash_[hash] == pas) {
          return true;
        }
        hash++;
      }
    }
    int HashPas(const std::pair <std::vector <int>, std::vector <int>>& pas, bool skip_hash_presence = false) {
      if (!skip_hash_presence && CheckHashPresence({ pas.second, pas.first })) {
        return HashPas({ pas.second, pas.first }, true);
      }
      std::vector <int> temp = pas.first;
      temp.insert(temp.end(), pas.second.begin(), pas.second.end());
      int hash = HashArray(temp);
      bool hash_found = false;
      while (!hash_found) {
        if (!this->reserved_pas_hash_.count(hash) || this->reserved_pas_hash_[hash] == pas) {
          this->reserved_pas_hash_[hash] = pas;
          hash_found = true;
        }
        else {
          hash++;
        }
      }
      return hash;
    }
    void OriginLinkPasInitialization(int origin_index, const std::pair <std::vector <int>, std::vector <int>>& found_pas) {
      int pas_hash = HashPas(found_pas);
      this->pas_users_[pas_hash].insert(origin_index);
      std::pair <int, int> term_links = { found_pas.first[found_pas.first.size() - 1], found_pas.second[found_pas.second.size() - 1] };
      if (term_links.first > term_links.second) {
        std::swap(term_links.first, term_links.second);
      }
      origin_link_pair_pas_[origin_index][term_links] = pas_hash;
      if (this->pas_users_[pas_hash].size() == 1) {
        delta_pas_queue.push({ PasDelta(pas_hash), pas_hash });
      }

      origin_corresponding_pas_set_[origin_index].insert(pas_hash);
      for (auto pas_link : found_pas.first) {
        this->origin_link_corresponding_pas_[origin_index][pas_link] = pas_hash;
      }
      for (auto pas_link : found_pas.second) {
        this->origin_link_corresponding_pas_[origin_index][pas_link] = pas_hash;
      }
      // requires reconsideration
      this->pas_flow_shift_starting_point_[pas_hash] = 1;
    }

    // Pas building based on the origin and considered link
    std::pair <std::vector <int>, std::vector <int>> OriginLinkPasConstruction(int origin_index, int link_index) {
      // contains nodes that are going to be processed through bfs
      std::queue <int> q;
      // contains nodes that are/were in q
      std::unordered_set <int> processed_nodes;
      // contains route from origin to current link termination
      std::unordered_set <int> origin_term_route;
      // contains for every node a link that was used in order to get to that link in bfs
      std::unordered_map <int, int> used_link;
      // origin_term_route building
      int now = this->network_.links()[link_index].term;
      while (now != origins_[origin_index]) {
        now = this->parent_origin_least_cost_routes_tree_[origin_index][now];
        origin_term_route.insert(now);
      }
      now = this->network_.links()[link_index].init;
      used_link[now] = link_index;
      q.push(now);
      processed_nodes.insert(now);
      // bfs in order to build pas
      std::pair <std::vector <int>, std::vector <int>> found_pas;
      bool fl = false;
      while (!fl) {
        now = q.front();
        q.pop();
        if (origin_term_route.count(now)) {
          found_pas = PasReconstruction(origin_index, now, link_index, origin_term_route, used_link);
          fl = true;
        }
        else {
          for (auto cur_link : this->network_.reverse_adjacency()[now]) {
            if ((bush_links_flows_[origin_index][cur_link] > 0) && (!processed_nodes.count(this->network_.links()[cur_link].init))) {
              processed_nodes.insert(this->network_.links()[cur_link].init);
              q.push(this->network_.links()[cur_link].init);
              used_link[this->network_.links()[cur_link].init] = cur_link;
            }
          }
        }
      }
      return found_pas;
    }

    bool CheckLinkPairPas(int origin_index, int link_1, int link_2) {
      if (link_1 > link_2) {
        std::swap(link_1, link_2);
      }
      if (origin_link_pair_pas_[origin_index].count({ link_1, link_2 })) {
        return true;
      }
      else {
        return false;
      }
    }

    bool LinkPasConstructionCondition(int origin_index, int link_index) {
      if (this->links_origin_least_cost_routes_tree_[origin_index].count(link_index)) {
        return false;
      }
      if (this->bush_links_flows_[origin_index][link_index] < this->computation_threshold_) {
        return false;
      }
      int term = this->network_.links()[link_index].term;
      int least_cost_routes_link = -1;
      for (int link_index : this->network_.adjacency()[parent_origin_least_cost_routes_tree_[origin_index][term]]) {
        if (this->network_.links()[link_index].term == term) {
          int least_cost_routes_link = link_index;
        }
      }
      if (CheckLinkPairPas(origin_index, link_index, least_cost_routes_link)) {
        return false;
      }
      return true;
    }

    // For every link that doesn't have a PAS corresponding to it builds a new PAS
    void SingleOriginLinksPasConstruction(int origin_index) {
      for (int link_index = 0; link_index < this->network_.number_of_links(); link_index++) {
        //(this->origin_link_corresponding_pas_[origin_index][link_index] == -1) &&
        if (LinkPasConstructionCondition(origin_index, link_index)) {
          OriginLinkPasInitialization(origin_index, OriginLinkPasConstruction(origin_index, link_index));
        }
      }
    }

    // Returns information about how much flow on PAS we have on every origin that is using it specifically
    // pair reflects the information about two parts of PAS, unordered map indicates how much flow origin with origin_index has
    std::pair <std::unordered_map <int, T>, std::unordered_map <int, T>> PasOriginsFlows(int pas_hash) {
      std::pair <std::unordered_map <int, T>, std::unordered_map <int, T>> result;
      const std::pair <std::vector <int>, std::vector <int>> pas = this->reserved_pas_hash_[pas_hash];
      for (auto origin_index : this->pas_users_[pas_hash]) {
        T min_flow = -1;
        for (auto link_index : pas.first) {
          if (min_flow == -1) {
            min_flow = this->bush_links_flows_[origin_index][link_index];
          }
          else {
            min_flow = std::min(min_flow, this->bush_links_flows_[origin_index][link_index]);
          }
        }
        result.first[origin_index] = min_flow;
      }
      for (auto origin_index : this->pas_users_[pas_hash]) {
        T min_flow = -1;
        for (auto link_index : pas.second) {
          if (min_flow == -1) {
            min_flow = this->bush_links_flows_[origin_index][link_index];
          }
          else {
            min_flow = std::min(min_flow, this->bush_links_flows_[origin_index][link_index]);
          }
        }
        result.second[origin_index] = min_flow;
      }
      return result;
    }

    std::pair <T, T> PasTotalAmountOfFlow(int pas_hash) {
      std::pair <T, T> total_amount = { 0, 0 };
      const std::pair <std::vector <int>, std::vector <int>> pas = this->reserved_pas_hash_[pas_hash];
      for (auto origin_index : this->pas_users_[pas_hash]) {
        T min_flow = -1;
        for (auto link_index : pas.first) {
          if (min_flow == -1) {
            min_flow = this->bush_links_flows_[origin_index][link_index];
          }
          else {
            min_flow = std::min(min_flow, this->bush_links_flows_[origin_index][link_index]);
          }
        }
        total_amount.first += min_flow;
      }
      for (auto origin_index : this->pas_users_[pas_hash]) {
        T min_flow = -1;
        for (auto link_index : pas.second) {
          if (min_flow == -1) {
            min_flow = this->bush_links_flows_[origin_index][link_index];
          }
          else {
            min_flow = std::min(min_flow, this->bush_links_flows_[origin_index][link_index]);
          }
        }
        total_amount.second += min_flow;
      }
      return total_amount;
    }

    // Returns true if delay with flow_shift of the first part of the PAS will be more then the second part
    bool PasFlowShiftResult(const int pas_hash, const std::pair <T, T> flow_shift) {
      std::pair <T, T> pas_delay = { 0, 0 };
      for (auto link_index : this->reserved_pas_hash_[pas_hash].first) {
        pas_delay.first += this->network_.links()[link_index].Delay(this->network_.links()[link_index].flow + flow_shift.first);
      }
      for (auto link_index : this->reserved_pas_hash_[pas_hash].second) {
        pas_delay.second += this->network_.links()[link_index].Delay(this->network_.links()[link_index].flow + flow_shift.second);
      }
      return (pas_delay.first > pas_delay.second);
    }

    std::unordered_map <int, T> PasProportion(const int pas_hash, const std::pair <T, T> flow_shift, const std::pair <T, T> total_pas_flow) {
      std::pair <std::unordered_map <int, T>, std::unordered_map <int, T>> pas_origins_flows = this->PasOriginsFlows(pas_hash);
      std::unordered_map <int, T> proportion;
      if (flow_shift.first < 0) {
        for (auto origin_index : this->pas_users_[pas_hash]) {
          proportion[origin_index] = pas_origins_flows.first[origin_index] / total_pas_flow.first;
        }
      }
      else {
        for (auto origin_index : this->pas_users_[pas_hash]) {
          proportion[origin_index] = pas_origins_flows.second[origin_index] / total_pas_flow.second;
        }
      }
      return proportion;
    }

    // For PAS with pas_hash shifts the amount of flow equal flow_shift
    // the first value of flow_shift reflects how much flow is going to to be shifted (if value is negative) / received (if positive)
    // from the first part of PAS
    // the second value of flow_shift reflects how much flow is going to to be received (if positive) / shifted (if value is negative)
    // to the second part of PAS
    // values in this pair are supposed to have a format {x, -x}
    void ImplementPasFlowShift(const int pas_hash, const std::pair <T, T> flow_shift, std::pair <T, T> total_pas_flow) {
      if (std::abs(flow_shift.first) < this->computation_threshold_) {
        return;
      }
      std::unordered_map <int, T> proportion = PasProportion(pas_hash, flow_shift, total_pas_flow);
      std::pair <std::vector <int>, std::vector <int>> cur_pas = this->reserved_pas_hash_[pas_hash];
      for (auto origin_index : this->pas_users_[pas_hash]) {
        for (auto link_index : cur_pas.first) {
          this->bush_links_flows_[origin_index][link_index] += proportion[origin_index] * flow_shift.first;
        }
      }
      for (auto origin_index : this->pas_users_[pas_hash]) {
        for (auto link_index : cur_pas.second) {
          this->bush_links_flows_[origin_index][link_index] += proportion[origin_index] * flow_shift.second;
        }
      }
      for (auto link_id : cur_pas.first) {
        this->network_.SetLinkFlow(link_id, flow_shift.first);
        //this->links_[link_index].flow += flow_shift.first;
      }
      for (auto link_id : cur_pas.second) {
        this->network_.SetLinkFlow(link_id, flow_shift.second);
        //this->links_[link_index].flow += flow_shift.second;
      }
    }

    void EliminatePas(const int pas_hash) {
      //std::cout << 1 << '\n';
      std::pair <std::vector <int>, std::vector <int>> cur_pas = this->reserved_pas_hash_[pas_hash];
      std::pair <int, int> term_links = { cur_pas.first[cur_pas.first.size() - 1], cur_pas.second[cur_pas.second.size() - 1] };
      if (term_links.first > term_links.second) {
        std::swap(term_links.first, term_links.second);
      }
      for (auto origin_index : pas_users_[pas_hash]) {
        for (auto link_index : cur_pas.first) {
          origin_link_corresponding_pas_[origin_index][link_index] = -1;
        }
        origin_corresponding_pas_set_[origin_index].erase(pas_hash);
        origin_link_pair_pas_[origin_index].erase(term_links);
      }
      for (auto origin_index : pas_users_[pas_hash]) {
        for (auto link_index : cur_pas.second) {
          origin_link_corresponding_pas_[origin_index][link_index] = -1;
        }
      }
      pas_users_.erase(pas_hash);
      // REQUIRES CONSIDERATION
      reserved_pas_hash_.erase(pas_hash);
    }

    // If direction is true than shift is performed from the first to the second part of PAS
    void ImplementFullPasFlowShift(const int pas_hash, const bool direction, const std::pair <T, T> total_pas_flow, bool elimination_flag) {
      std::pair <std::vector <int>, std::vector <int>> cur_pas = this->reserved_pas_hash_[pas_hash];
      std::pair <T, T> flow_shift;
      if (direction) {
        flow_shift = { -total_pas_flow.first, total_pas_flow.first };
      }
      else {
        flow_shift = { total_pas_flow.second, -total_pas_flow.second };
      }
      if (std::abs(flow_shift.first)) {
        ImplementPasFlowShift(pas_hash, flow_shift, total_pas_flow);
      }
      if (elimination_flag) {
        EliminatePas(pas_hash);
      }
    }

    bool TryFullPasFlowShift(const int pas_hash, std::pair <T, T> total_flow, bool elimination_flag) {
      if (PasFlowShiftResult(pas_hash, { -total_flow.first, total_flow.first })) {
        ImplementFullPasFlowShift(pas_hash, true, total_flow, elimination_flag);
        return true;
      }
      if (!PasFlowShiftResult(pas_hash, { total_flow.second, -total_flow.second })) {
        ImplementFullPasFlowShift(pas_hash, false, total_flow, elimination_flag);
        return true;
      }
      return false;
    }

    // Provides another format for PasFlowShiftResult
    bool DirectionPasFlowShiftResult(const int pas_hash, const T flow_shift, bool direction) {
      if (direction) {
        return PasFlowShiftResult(pas_hash, { -flow_shift, flow_shift });
      }
      else {
        return PasFlowShiftResult(pas_hash, { flow_shift, -flow_shift });
      }
    }

    // 
    std::unordered_set <int> PasUsersToIgnore(const int pas_hash) {
      std::unordered_set <int> users_to_ignore;
      std::vector <int> pas_part;
      if (PasFlowShiftResult(pas_hash, { 0, 0 })) {
        pas_part = this->reserved_pas_hash_[pas_hash].first;
      }
      else {
        pas_part = this->reserved_pas_hash_[pas_hash].second;
      }
      for (auto origin_index : this->pas_users_[pas_hash]) {
        for (auto link_index : pas_part) {
          if (this->bush_links_flows_[origin_index][link_index] < computation_threshold_) {
            users_to_ignore.insert(origin_index);
            break;
          }
        }
      }
      return users_to_ignore;
    }

    // Flow shift under chosen PAS with "optimal" step size
    void PasFLowShift(const int pas_hash, bool elimination_flag) {
      //unordered_map <int> pas_users_to_ignore = PasUsersToIgnore(pas_hash);
      //if (elimination_flag) {
      //
      //}
      std::pair <T, T> total_pas_flow = this->PasTotalAmountOfFlow(pas_hash);
      if (TryFullPasFlowShift(pas_hash, total_pas_flow, elimination_flag)) {
        return;
      }
      if (shift_method_ == nullptr) {
        throw std::runtime_error("shift_method_ is not initialized");
      }
      std::pair <T, T> flow_shift = shift_method_->FlowShift(this->pas_flow_shift_starting_point_[pas_hash], this->reserved_pas_hash_[pas_hash], total_pas_flow);
      //cout << flow_shift.first << ' ' << flow_shift.second << '\n';
      this->pas_flow_shift_starting_point_[pas_hash] = abs(flow_shift.first);
      ImplementPasFlowShift(pas_hash, flow_shift, total_pas_flow);
    }

    void EliminatePasOriginUser(int pas_hash, int origin_index) {

    }

    // requires update
    // TOO SLOW SOLUTION
    void EliminationPasNonLeastCostRoutesTree() {
      for (int origin_index = 0; origin_index < this->number_of_origins_; origin_index++) {
        for (int link_index = 0; link_index < this->number_of_links_; link_index++) {
          origin_link_corresponding_pas_[origin_index][link_index] = -1;  
        }
        origin_corresponding_pas_set_[origin_index].clear();
      }
      pas_users_.clear();
      reserved_pas_hash_.clear();
      pas_flow_shift_starting_point_.clear();
    }

    std::vector <int> GetPasSubset() {
      std::vector <int> pas_subset;
      for (const auto& [pas_hash, pas] : this->reserved_pas_hash_) {
        pas_subset.push_back(pas_hash);
      }
      std::random_device rd;
      std::mt19937 g(rd());
      std::shuffle(pas_subset.begin(), pas_subset.end(), g);
      return pas_subset;
    }

    void DeltaPasQueuePreparation() {
      while (!delta_pas_queue.empty()) {
        delta_pas_queue.pop();
      }
      pas_processed_last_queue_update = 0;
    }

    T PasDelta(int pas) {
      return std::abs(Link<T>::GetLinksDelay(this->network_.links(), reserved_pas_hash_[pas].first) - Link<T>::GetLinksDelay(this->network_.links(), reserved_pas_hash_[pas].second));
    }

    void UpdateDeltaPasQueue() {
      while (!delta_pas_queue.empty()) {
        delta_pas_queue.pop();
      }
      pas_processed_last_queue_update = 0;
      for (const auto& now : reserved_pas_hash_) {
        delta_pas_queue.push({ PasDelta(now.first), now.first });
      }
    }

    bool DeltaPasQueuePushRequirement(int pas_hash) {
      if (!reserved_pas_hash_.count(pas_hash)) {
        return false;
      }
      std::pair <T, T> pas_flow = PasTotalAmountOfFlow(pas_hash);
      if (pas_flow.first < computation_threshold_ || pas_flow.second < computation_threshold_) {
        return false;
      }
      return true;
    }

    void PasProcessing(bool elimination_flag) {
      if (pas_processed_last_queue_update > reserved_pas_hash_.size() / 2) {
        UpdateDeltaPasQueue();
      }
      if (delta_pas_queue.empty()) {
        return;
      }
      pas_processed_last_queue_update++;
      int pas_hash = delta_pas_queue.top().second;
      delta_pas_queue.pop();
      PasFLowShift(pas_hash, elimination_flag);
      if (DeltaPasQueuePushRequirement(pas_hash)) {
        delta_pas_queue.push({ PasDelta(pas_hash), pas_hash });
      }
    }

    void EquilibrationIteration(int iteration_number) {
      //EliminationPasNonLeastCostRoutesTree();
      UpdateDeltaPasQueue();
      std::vector <int> origin_order(number_of_origins_);
      for (int origin_index = 0; origin_index < this->number_of_origins_; origin_index++) {
        origin_order[origin_index] = origin_index;
      }
      std::random_device rd;
      std::mt19937 g(rd());
      std::shuffle(origin_order.begin(), origin_order.end(), g);
      for (int origin_index : origin_order) {
        SingleOriginRemoveAllCyclicFlows(origin_index);
        BuildSingleOriginLeastCostRoutesTree(origin_index);
        SingleOriginLinksPasConstruction(origin_index);
        for (int cnt = 0; cnt < 1; cnt++) {
          PasProcessing(false);
          //this->AllPasOriginCheck();
          //std::vector <int> pas_subset = GetPasSubset();
          //for (auto pas_hash : origin_corresponding_pas_set_[origin_index]) {
          //  pas_subset.push_back(pas_hash);
          //}
          //std::vector <int> pas_subset = GetPasSubset();
          //for (auto pas_hash : pas_subset) {
          //for (auto pas_hash : origin_corresponding_pas_set_[origin_index]) {
          //  PasFLowShift(pas_hash, false);
          //}
          this->statistics_recorder_.RecordStatistics();
        }
        //cout << this->reserved_pas_hash_.size() << ' ' << this->ObjectiveFunction() << '\n';
        this->statistics_recorder_.RecordStatistics();
      }
      for (int cnt = 0; cnt < iteration_number * 5; cnt++) {
        //std::vector <int> pas_subset = GetPasSubset();
        for (int i = 0; i < reserved_pas_hash_.size(); i++) {
        //for (int i = 0; i < pas_subset.size(); i++) {
          //PasFLowShift(pas_subset[i], true);
          PasProcessing(true);
        }
        this->statistics_recorder_.RecordStatistics();
      }
      // cout << this->ObjectiveFunction() << '\n';
      this->statistics_recorder_.RecordStatistics();
    }

    const std::vector <Link<T>>& GetLinksRef() {
      return this->network_.links();
    }
  };
}

#endif
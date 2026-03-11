#ifndef TAPAS_APPROACH_H
#define TAPAS_APPROACH_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <string>
#include <random>
#include <algorithm>
#include <map>
#include "../common/TrafficAssignmentApproach.h"
#include "./components/TapasShiftMethod.h"
#include "../../core/Network.h"
#include "TapasShiftMethodFactory.h"

namespace TrafficAssignment {
  /**
   * @brief TAPAS (Traffic Assignment by Paired Alternative Segments) algorithm.
   *
   * Implements the bush-based TAPAS algorithm from Bar-Gera (2010) "Traffic assignment
   * by paired alternative segments", Fig. 7. Uses origin decomposition: maintains
   * per-origin bush flows and constructs Paired Alternative Segments (PAS) to
   * redistribute flow toward user equilibrium.
   *
   * Key algorithmic structure (per iteration):
   * 1. For each origin (random order): build Dijkstra least-cost tree, discover new
   *    PAS via backward BFS (Fig. 6), perform local PAS flow shifts
   * 2. Global PAS processing: shift flows on all active PAS using the configured
   *    shift method (Table 3), eliminate converged PAS
   *
   * Achieves very high precision (~1e-15 RGAP). Includes stall detection: after 3
   * non-improving RGAP checks, restarts bush flows via AON and eliminates all PAS.
   *
   * @tparam T Numeric type for flow computations.
   */
  template <typename T>
  class TapasApproach : public TrafficAssignmentApproach <T> {
  public:
    TapasApproach(Network<T>& network,
                  T alpha = 1e-6,
                  std::string shift_method_name = "NewtonStep",
                  int max_iterations = 200,
                  int pas_per_origin = 1,
                  int pas_multiplier = 5,
                  int rgap_check_interval = 5)
      : TrafficAssignmentApproach<T>::TrafficAssignmentApproach(network, alpha),
        max_iterations_(max_iterations),
        pas_per_origin_(pas_per_origin),
        pas_multiplier_(pas_multiplier),
        rgap_check_interval_(rgap_check_interval),
        rng_(std::random_device{}()) {

      shift_method_ = TapasShiftMethodFactory<T>::GetInstance().Create(shift_method_name, this->GetLinksRef(), computation_threshold_);
      this->shift_method_name_ = shift_method_name;

      for (int i = 0; i < this->network_.number_of_zones(); i++) {
        if (this->network_.origin_info()[i].size() > 0) {
          origins_.push_back(i);
        }
      }
      number_of_origins_ = origins_.size();
      bush_links_flows_.resize(this->number_of_origins_);
      origin_link_corresponding_pas_.resize(this->number_of_origins_);
      origin_corresponding_pas_set_.resize(this->number_of_origins_);
      links_origin_least_cost_routes_tree_.resize(this->number_of_origins_,
          std::vector<bool>(this->network_.number_of_links(), false));
      parent_origin_least_cost_routes_tree_.resize(this->number_of_origins_);
      origin_link_pair_pas_.resize(this->number_of_origins_);

      for (int origin_index = 0; origin_index < this->number_of_origins_; origin_index++) {
        bush_links_flows_[origin_index].resize(this->network_.number_of_links(), 0);
        origin_link_corresponding_pas_[origin_index].resize(this->network_.number_of_links(), -1);
        parent_origin_least_cost_routes_tree_[origin_index].resize(this->network_.number_of_nodes(), -1);
      }
    }

    ~TapasApproach() = default;

    void Reset() override {
      TrafficAssignmentApproach<T>::Reset();
      ClearInternalState();
    }

    /**
     * @brief Main TAPAS loop: AON init, then equilibration iterations with stall detection.
     *
     * After AON initialization, runs up to max_iterations equilibration iterations.
     * Every rgap_check_interval iterations, checks RGAP for convergence. If RGAP stalls
     * (< 0.1% improvement for 3 consecutive checks), restarts bush flows and PAS.
     */
    void ComputeTrafficFlows(bool statistics_recording = false) override {
      this->network_.Reset();
      ClearInternalState();
      if (statistics_recording) {
        this->statistics_recorder_.StartRecording(this->GetApproachName());
      }
      AllOrNothingAssignment();
      if (statistics_recording) {
        this->statistics_recorder_.RecordStatistics();
      }
      T best_rgap = T(1);
      int stall_count = 0;
      for (int i = 0; i < max_iterations_; i++) {
        EquilibrationIteration(i + 1, statistics_recording);
        if ((i + 1) % rgap_check_interval_ == 0) {
          T rgap = this->network_.RelativeGap();
          if (rgap >= 0 && rgap <= this->alpha_) {
            break;
          }
          // Detect stalling: rgap hasn't improved by at least 0.1%
          if (rgap < best_rgap * T(0.999)) {
            best_rgap = rgap;
            stall_count = 0;
          } else {
            stall_count++;
          }
          // After 3 stalled checks, restart bush flows + PAS
          if (stall_count >= 3) {
            RestartBushFlows();
            EliminateAllPas();
            stall_count = 0;
          }
        }
      }
    }
    
    std::string GetApproachName() override {
      return "Tapas" + shift_method_name_;
    }

  protected:
    std::string shift_method_name_;
    int max_iterations_;
    int pas_per_origin_;
    int pas_multiplier_;
    int rgap_check_interval_;

    int number_of_origins_;

    std::mt19937 rng_; ///< RNG for randomizing origin processing order.

    /// @brief Indices of active origins (zones with nonzero demand).
    std::vector <int> origins_;

    std::vector <std::vector <int>> destinations_for_origin_;

    /// @brief Origin-decomposed link flows: bush_links_flows_[origin][link] = f_{pa}.
    /// Each origin maintains its own bush (subgraph of links with positive flow).
    /// See Bar-Gera (2010) Section 5.
    std::vector <std::vector <T>> bush_links_flows_;

    /// @brief Maps each (origin, link) to the PAS hash that covers it, or -1 if none.
    std::vector <std::vector <int>> origin_link_corresponding_pas_;

    /// @brief Maps PAS hash -> set of origin indices that use this PAS for flow redistribution.
    std::unordered_map <int, std::unordered_set <int>> pas_users_;

    /// @brief Per-origin Dijkstra least-cost tree: parent_[origin][node] = parent node.
    std::vector <std::vector <int>> parent_origin_least_cost_routes_tree_;

    /// @brief Per-origin least-cost tree membership: links_[origin][link] = true if link is in tree.
    std::vector <std::vector<bool>> links_origin_least_cost_routes_tree_;

    /// @brief PAS hash table: maps hash -> (segment_1 links, segment_2 links).
    /// Uses open addressing with linear probing for collision resolution.
    std::unordered_map <int, std::pair <std::vector <int>, std::vector <int>>> reserved_pas_hash_;

    /// @brief Remembers last shift amount per PAS (warm-start hint for line search).
    std::unordered_map <int, T> pas_flow_shift_starting_point_;

    /// @brief Per-origin set of PAS hashes associated with this origin.
    std::vector <std::unordered_set <int>> origin_corresponding_pas_set_;

    std::unique_ptr <TapasShiftMethod <T>> shift_method_;

    const T computation_threshold_ = 1e-10;
    static constexpr int max_hash_probes_ = 10000;

    /// @brief Priority queue of PAS ordered by cost difference (largest first).
    /// Entries may be stale (PAS eliminated); checked before processing.
    std::priority_queue <std::pair <T, int>> delta_pas_queue;

    /// @brief Counter for queue staleness detection (triggers rebuild when > size/2).
    int pas_processed_last_queue_update;

    /// @brief Per-origin map from terminal link pair -> PAS hash (deduplication).
    std::vector <std::map <std::pair <int, int>, int>> origin_link_pair_pas_;

    void ClearInternalState() {
      for (int i = 0; i < number_of_origins_; i++) {
        std::fill(bush_links_flows_[i].begin(), bush_links_flows_[i].end(), T(0));
        std::fill(origin_link_corresponding_pas_[i].begin(), origin_link_corresponding_pas_[i].end(), -1);
        origin_corresponding_pas_set_[i].clear();
        std::fill(links_origin_least_cost_routes_tree_[i].begin(), links_origin_least_cost_routes_tree_[i].end(), false);
        std::fill(parent_origin_least_cost_routes_tree_[i].begin(), parent_origin_least_cost_routes_tree_[i].end(), -1);
        origin_link_pair_pas_[i].clear();
      }
      reserved_pas_hash_.clear();
      pas_users_.clear();
      pas_flow_shift_starting_point_.clear();
      while (!delta_pas_queue.empty()) {
        delta_pas_queue.pop();
      }
    }

    /// @brief All-or-nothing (AON) initialization: assigns full demand to shortest paths
    /// and populates per-origin bush link flows.
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
      auto users_it = this->pas_users_.find(pas_hash);
      if (users_it != this->pas_users_.end() && users_it->second.count(origin_index)) {
        return;
      }
      auto pas_it = this->reserved_pas_hash_.find(pas_hash);
      if (pas_it == this->reserved_pas_hash_.end()) {
        return;
      }
      const auto& pas = pas_it->second;
      if (Link<T>::GetLinksDelay(this->network_.links(), pas.first) >
        Link<T>::GetLinksDelay(this->network_.links(), pas.second)) {
        for (auto now : pas.first) {
          if (this->bush_links_flows_[origin_index][now] == 0) {
            return;
          }
        }
      }
      else {
        for (auto now : pas.second) {
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
            temp = this->network_.links()[link_used[temp]].init;
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

    /**
     * @brief Removes all cyclic flows from a single origin's bush via DFS.
     *
     * Performs DFS on links with positive bush flow. When a back edge is found
     * (indicating a cycle), removes the minimum flow around the cycle. Repeats
     * until no cycles remain. Ensures bush remains a DAG.
     */
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

    /**
     * @brief Builds Dijkstra least-cost tree for a single origin.
     *
     * Populates parent_origin_least_cost_routes_tree_[origin] with parent nodes
     * and links_origin_least_cost_routes_tree_[origin] with tree membership flags.
     * This tree is used for PAS construction: links not in the tree but carrying
     * bush flow are candidates for new PAS. See Bar-Gera (2010) Section 6.
     */
    void BuildSingleOriginLeastCostRoutesTree(int origin_index) {
      std::fill(this->links_origin_least_cost_routes_tree_[origin_index].begin(),
                this->links_origin_least_cost_routes_tree_[origin_index].end(), false);
      int current_origin = this->origins_[origin_index];
      // Pairs in q contain delay from origin and node index
      std::priority_queue <std::pair <T, int>, std::vector <std::pair <T, int>>, std::greater <std::pair <T, int>>> q;
      int u;
      int link_index;
      int origin = this->origins_[origin_index];
      T cur_delay = 0;
      q.push({ 0, -1 });
      std::vector<bool> processed(this->network_.number_of_nodes(), false);
      while (!q.empty()) {
        cur_delay = q.top().first;
        link_index = q.top().second;
        if (q.top().second != -1) {
          u = this->network_.links()[link_index].term;
        }
        else {
          u = origin;
        }
        q.pop();
        if (processed[u]) {
          continue;
        }
        if (u != current_origin) {
          this->parent_origin_least_cost_routes_tree_[origin_index][u] = this->network_.links()[link_index].init;
          this->links_origin_least_cost_routes_tree_[origin_index][link_index] = true;
        }

        processed[u] = true;
        for (const auto now : this->network_.adjacency()[u]) {
          if (!processed[this->network_.links()[now].term]) {
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
        if (prev_node < 0) {
          return {{}, {}};
        }
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
      for (int probe = 0; probe < max_hash_probes_; ++probe) {
        if (!this->reserved_pas_hash_.count(hash)) {
          return false;
        }
        if (this->reserved_pas_hash_.at(hash) == pas) {
          return true;
        }
        hash++;
      }
      throw std::runtime_error("CheckHashPresence: exceeded max probe limit (" + std::to_string(max_hash_probes_) + ")");
    }
    int HashPas(const std::pair <std::vector <int>, std::vector <int>>& pas, bool skip_hash_presence = false) {
      if (!skip_hash_presence && CheckHashPresence({ pas.second, pas.first })) {
        return HashPas({ pas.second, pas.first }, true);
      }
      std::vector <int> temp = pas.first;
      temp.insert(temp.end(), pas.second.begin(), pas.second.end());
      int hash = HashArray(temp);
      for (int probe = 0; probe < max_hash_probes_; ++probe) {
        if (!this->reserved_pas_hash_.count(hash) || this->reserved_pas_hash_.at(hash) == pas) {
          this->reserved_pas_hash_[hash] = pas;
          return hash;
        }
        hash++;
      }
      throw std::runtime_error("HashPas: exceeded max probe limit (" + std::to_string(max_hash_probes_) + ")");
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

    /**
     * @brief Constructs a PAS for a given origin and non-tree link via backward BFS.
     *
     * Starting from the non-tree link's init node, performs BFS backward along bush
     * links (positive origin flow) until reaching a node on the least-cost tree path
     * to the link's term node. The two paths from that junction node to the term node
     * form the PAS segments. See Bar-Gera (2010) Fig. 6.
     *
     * @param origin_index Index into origins_ array.
     * @param link_index The non-tree link that triggers PAS construction.
     * @return (least_cost_segment, non_least_cost_segment), or empty if construction fails.
     */
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
        int parent = this->parent_origin_least_cost_routes_tree_[origin_index][now];
        if (parent < 0) {
          return {{}, {}};
        }
        now = parent;
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
        if (q.empty()) {
          return {{}, {}};
        }
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

    /**
     * @brief Checks whether a link qualifies for PAS construction.
     *
     * A link qualifies if it: (1) is not in the least-cost tree (cost-effective check),
     * (2) carries positive bush flow (flow-effective check), and (3) doesn't already
     * have a PAS for the same terminal link pair. See Bar-Gera (2010) Section 6.1-6.2.
     */
    bool LinkPasConstructionCondition(int origin_index, int link_index) {
      if (this->links_origin_least_cost_routes_tree_[origin_index][link_index]) {
        return false;
      }
      if (this->bush_links_flows_[origin_index][link_index] < this->computation_threshold_) {
        return false;
      }
      int term = this->network_.links()[link_index].term;
      int parent_node = parent_origin_least_cost_routes_tree_[origin_index][term];
      if (parent_node < 0) {
        return false;
      }
      int least_cost_routes_link = -1;
      for (int adj_link : this->network_.adjacency()[parent_node]) {
        if (this->network_.links()[adj_link].term == term) {
          least_cost_routes_link = adj_link;
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
        if (LinkPasConstructionCondition(origin_index, link_index)) {
          auto pas = OriginLinkPasConstruction(origin_index, link_index);
          if (!pas.first.empty() && !pas.second.empty()) {
            OriginLinkPasInitialization(origin_index, pas);
          }
        }
      }
    }

    struct PasFlowInfo {
      std::pair<T, T> totals;
      std::pair<std::unordered_map<int, T>, std::unordered_map<int, T>> per_origin;
    };

    PasFlowInfo ComputePasFlowInfo(int pas_hash) {
      PasFlowInfo info;
      info.totals = {0, 0};
      const auto& pas = this->reserved_pas_hash_.at(pas_hash);
      for (auto origin_index : this->pas_users_.at(pas_hash)) {
        T min_flow = -1;
        for (auto link_index : pas.first) {
          if (min_flow == -1) {
            min_flow = this->bush_links_flows_[origin_index][link_index];
          }
          else {
            min_flow = std::min(min_flow, this->bush_links_flows_[origin_index][link_index]);
          }
        }
        info.per_origin.first[origin_index] = min_flow;
        info.totals.first += min_flow;
      }
      for (auto origin_index : this->pas_users_.at(pas_hash)) {
        T min_flow = -1;
        for (auto link_index : pas.second) {
          if (min_flow == -1) {
            min_flow = this->bush_links_flows_[origin_index][link_index];
          }
          else {
            min_flow = std::min(min_flow, this->bush_links_flows_[origin_index][link_index]);
          }
        }
        info.per_origin.second[origin_index] = min_flow;
        info.totals.second += min_flow;
      }
      return info;
    }

    std::pair <T, T> PasTotalAmountOfFlow(int pas_hash) {
      std::pair <T, T> total_amount = { 0, 0 };
      const auto& pas = this->reserved_pas_hash_.at(pas_hash);
      for (auto origin_index : this->pas_users_.at(pas_hash)) {
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
      for (auto origin_index : this->pas_users_.at(pas_hash)) {
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
      const auto& pas = this->reserved_pas_hash_.at(pas_hash);
      for (auto link_index : pas.first) {
        pas_delay.first += this->network_.links()[link_index].Delay(this->network_.links()[link_index].flow + flow_shift.first);
      }
      for (auto link_index : pas.second) {
        pas_delay.second += this->network_.links()[link_index].Delay(this->network_.links()[link_index].flow + flow_shift.second);
      }
      return (pas_delay.first > pas_delay.second);
    }

    /**
     * @brief Computes per-origin proportional share for a PAS flow shift.
     *
     * Each origin's share is proportional to its flow on the segment being reduced.
     * See Bar-Gera (2010) Section 7, Eq. 19-31.
     */
    std::unordered_map <int, T> PasProportion(const int pas_hash, const std::pair <T, T> flow_shift, const std::pair <T, T> total_pas_flow,
        const std::pair <std::unordered_map <int, T>, std::unordered_map <int, T>>& pas_origins_flows) {
      std::unordered_map <int, T> proportion;
      if (flow_shift.first < 0) {
        for (auto origin_index : this->pas_users_.at(pas_hash)) {
          proportion[origin_index] = pas_origins_flows.first.at(origin_index) / total_pas_flow.first;
        }
      }
      else {
        for (auto origin_index : this->pas_users_.at(pas_hash)) {
          proportion[origin_index] = pas_origins_flows.second.at(origin_index) / total_pas_flow.second;
        }
      }
      return proportion;
    }

    /**
     * @brief Applies a flow shift to a PAS, distributing proportionally across origins.
     *
     * The shift is split among participating origins proportional to their current
     * flow on the losing segment (the segment flow is being taken from).
     * Updates both per-origin bush flows and aggregate network link flows.
     * See Bar-Gera (2010) Section 7, Eq. 19-31.
     *
     * @param flow_shift (delta_1, delta_2) where delta_1 + delta_2 = 0.
     *   Negative value = flow removed from that segment; positive = flow added.
     */
    void ImplementPasFlowShift(const int pas_hash, const std::pair <T, T> flow_shift, std::pair <T, T> total_pas_flow,
        const std::pair <std::unordered_map <int, T>, std::unordered_map <int, T>>& pas_origins_flows) {
      if (std::abs(flow_shift.first) < this->computation_threshold_) {
        return;
      }
      std::unordered_map <int, T> proportion = PasProportion(pas_hash, flow_shift, total_pas_flow, pas_origins_flows);
      const auto& cur_pas = this->reserved_pas_hash_.at(pas_hash);
      for (auto origin_index : this->pas_users_.at(pas_hash)) {
        for (auto link_index : cur_pas.first) {
          this->bush_links_flows_[origin_index][link_index] += proportion[origin_index] * flow_shift.first;
        }
      }
      for (auto origin_index : this->pas_users_.at(pas_hash)) {
        for (auto link_index : cur_pas.second) {
          this->bush_links_flows_[origin_index][link_index] += proportion[origin_index] * flow_shift.second;
        }
      }
      for (auto link_id : cur_pas.first) {
        this->network_.SetLinkFlow(link_id, flow_shift.first);
      }
      for (auto link_id : cur_pas.second) {
        this->network_.SetLinkFlow(link_id, flow_shift.second);
      }
    }

    /// @brief Removes a PAS from all data structures (hash table, user sets, origin maps).
    void EliminatePas(const int pas_hash) {
      auto pas_it = this->reserved_pas_hash_.find(pas_hash);
      if (pas_it == this->reserved_pas_hash_.end()) {
        pas_users_.erase(pas_hash);
        return;
      }
      const auto& cur_pas = pas_it->second;
      if (cur_pas.first.empty() || cur_pas.second.empty()) {
        pas_users_.erase(pas_hash);
        reserved_pas_hash_.erase(pas_it);
        return;
      }
      std::pair <int, int> term_links = { cur_pas.first[cur_pas.first.size() - 1], cur_pas.second[cur_pas.second.size() - 1] };
      if (term_links.first > term_links.second) {
        std::swap(term_links.first, term_links.second);
      }
      auto users_it = pas_users_.find(pas_hash);
      if (users_it == pas_users_.end()) {
        reserved_pas_hash_.erase(pas_it);
        return;
      }
      for (auto origin_index : users_it->second) {
        for (auto link_index : cur_pas.first) {
          origin_link_corresponding_pas_[origin_index][link_index] = -1;
        }
        origin_corresponding_pas_set_[origin_index].erase(pas_hash);
        origin_link_pair_pas_[origin_index].erase(term_links);
      }
      for (auto origin_index : users_it->second) {
        for (auto link_index : cur_pas.second) {
          origin_link_corresponding_pas_[origin_index][link_index] = -1;
        }
      }
      pas_users_.erase(pas_hash);
      // REQUIRES CONSIDERATION
      reserved_pas_hash_.erase(pas_hash);
    }

    // If direction is true than shift is performed from the first to the second part of PAS
    void ImplementFullPasFlowShift(const int pas_hash, const bool direction, const std::pair <T, T> total_pas_flow,
        const std::pair <std::unordered_map <int, T>, std::unordered_map <int, T>>& pas_origins_flows, bool elimination_flag) {
      const auto& cur_pas = this->reserved_pas_hash_.at(pas_hash);
      std::pair <T, T> flow_shift;
      if (direction) {
        flow_shift = { -total_pas_flow.first, total_pas_flow.first };
      }
      else {
        flow_shift = { total_pas_flow.second, -total_pas_flow.second };
      }
      if (std::abs(flow_shift.first)) {
        ImplementPasFlowShift(pas_hash, flow_shift, total_pas_flow, pas_origins_flows);
      }
      if (elimination_flag) {
        EliminatePas(pas_hash);
      }
    }

    bool TryFullPasFlowShift(const int pas_hash, std::pair <T, T> total_flow,
        const std::pair <std::unordered_map <int, T>, std::unordered_map <int, T>>& pas_origins_flows, bool elimination_flag) {
      if (PasFlowShiftResult(pas_hash, { -total_flow.first, total_flow.first })) {
        ImplementFullPasFlowShift(pas_hash, true, total_flow, pas_origins_flows, elimination_flag);
        return true;
      }
      if (!PasFlowShiftResult(pas_hash, { total_flow.second, -total_flow.second })) {
        ImplementFullPasFlowShift(pas_hash, false, total_flow, pas_origins_flows, elimination_flag);
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
      const auto& pas = this->reserved_pas_hash_.at(pas_hash);
      const std::vector <int>& pas_part = PasFlowShiftResult(pas_hash, { 0, 0 }) ? pas.first : pas.second;
      for (auto origin_index : this->pas_users_.at(pas_hash)) {
        for (auto link_index : pas_part) {
          if (this->bush_links_flows_[origin_index][link_index] < computation_threshold_) {
            users_to_ignore.insert(origin_index);
            break;
          }
        }
      }
      return users_to_ignore;
    }

    /**
     * @brief Performs flow shift on a PAS using the configured shift method.
     *
     * First attempts a full flow transfer (moving all flow from the costlier segment).
     * If that overshoots, delegates to the shift method (Newton step or line search)
     * for a partial shift. See Bar-Gera (2010) Table 3.
     */
    void PasFlowShift(const int pas_hash, bool elimination_flag) {
      auto flow_info = this->ComputePasFlowInfo(pas_hash);
      if (TryFullPasFlowShift(pas_hash, flow_info.totals, flow_info.per_origin, elimination_flag)) {
        return;
      }
      if (shift_method_ == nullptr) {
        throw std::runtime_error("shift_method_ is not initialized");
      }
      std::pair <T, T> flow_shift = shift_method_->FlowShift(this->pas_flow_shift_starting_point_[pas_hash], this->reserved_pas_hash_.at(pas_hash), flow_info.totals);
      this->pas_flow_shift_starting_point_[pas_hash] = abs(flow_shift.first);
      ImplementPasFlowShift(pas_hash, flow_shift, flow_info.totals, flow_info.per_origin);
    }



    std::vector <int> GetPasSubset() {
      std::vector <int> pas_subset;
      for (const auto& [pas_hash, pas] : this->reserved_pas_hash_) {
        pas_subset.push_back(pas_hash);
      }
      std::shuffle(pas_subset.begin(), pas_subset.end(), rng_);
      return pas_subset;
    }

    void DeltaPasQueuePreparation() {
      while (!delta_pas_queue.empty()) {
        delta_pas_queue.pop();
      }
      pas_processed_last_queue_update = 0;
    }

    /// @brief Returns the absolute cost difference between the two PAS segments.
    T PasDelta(int pas) {
      const auto& p = reserved_pas_hash_.at(pas);
      return std::abs(Link<T>::GetLinksDelay(this->network_.links(), p.first) - Link<T>::GetLinksDelay(this->network_.links(), p.second));
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

    /// @brief Processes the highest-priority PAS from the queue (shift flow, optionally eliminate).
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
      // Skip stale entries (PAS was eliminated but entry remained in queue)
      if (!reserved_pas_hash_.count(pas_hash)) {
        return;
      }
      PasFlowShift(pas_hash, elimination_flag);
      if (DeltaPasQueuePushRequirement(pas_hash)) {
        delta_pas_queue.push({ PasDelta(pas_hash), pas_hash });
      } else if (elimination_flag && reserved_pas_hash_.count(pas_hash)) {
        EliminatePas(pas_hash);
      }
    }

    /**
     * @brief Stall recovery: replaces all bush flows with AON on current shortest paths.
     *
     * Removes old per-origin bush flows from the network, then recomputes shortest
     * paths and assigns full demand. Maintains flow conservation per origin.
     * Triggered when RGAP stalls for 3 consecutive checks.
     */
    void RestartBushFlows() {
      for (int origin_index = 0; origin_index < this->number_of_origins_; origin_index++) {
        // Remove old bush flows from network
        for (int link_id = 0; link_id < this->network_.number_of_links(); link_id++) {
          if (bush_links_flows_[origin_index][link_id] > 0) {
            this->network_.SetLinkFlow(link_id, -bush_links_flows_[origin_index][link_id]);
            bush_links_flows_[origin_index][link_id] = 0;
          }
        }
        // Add new AON flows on current shortest paths
        auto shortest_routes = this->network_.ComputeSingleOriginBestRoutes(origins_[origin_index]);
        for (const auto& route : shortest_routes) {
          T demand = this->network_.od_pairs()[route.first].GetDemand();
          for (int link_id : route.second) {
            bush_links_flows_[origin_index][link_id] += demand;
            this->network_.SetLinkFlow(link_id, demand);
          }
        }
      }
    }

    /// @brief Eliminates all PAS (stall recovery companion to RestartBushFlows).
    void EliminateAllPas() {
      std::vector<int> all_hashes;
      for (const auto& [hash, _] : reserved_pas_hash_) {
        all_hashes.push_back(hash);
      }
      for (int hash : all_hashes) {
        EliminatePas(hash);
      }
      pas_flow_shift_starting_point_.clear();
      DeltaPasQueuePreparation();
    }

    /**
     * @brief One full equilibration iteration. See Bar-Gera (2010) Fig. 7.
     *
     * Phase 1 (per-origin, random order): build least-cost tree, construct new PAS,
     * perform pas_per_origin local PAS flow shifts.
     * Phase 2 (global): process all PAS multiple times with elimination enabled,
     * removing converged PAS from the active set.
     */
    void EquilibrationIteration(int iteration_number, bool statistics_recording = false) {
      UpdateDeltaPasQueue();
      std::vector <int> origin_order(number_of_origins_);
      for (int origin_index = 0; origin_index < this->number_of_origins_; origin_index++) {
        origin_order[origin_index] = origin_index;
      }
      std::shuffle(origin_order.begin(), origin_order.end(), rng_);
      for (int origin_index : origin_order) {
        BuildSingleOriginLeastCostRoutesTree(origin_index);
        SingleOriginLinksPasConstruction(origin_index);
        for (int cnt = 0; cnt < pas_per_origin_; cnt++) {
          PasProcessing(false);
        }
      }
      for (int cnt = 0; cnt < pas_multiplier_ * 4; cnt++) {
        for (std::size_t i = 0; i < reserved_pas_hash_.size(); i++) {
          PasProcessing(true);
        }
      }
      if (statistics_recording) {
        this->statistics_recorder_.RecordStatistics();
      }
    }

    const std::vector <Link<T>>& GetLinksRef() {
      return this->network_.links();
    }
  };
}

#endif

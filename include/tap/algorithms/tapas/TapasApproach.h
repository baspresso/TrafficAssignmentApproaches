#ifndef TAPAS_APPROACH_H
#define TAPAS_APPROACH_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <string>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include "../common/TrafficAssignmentApproach.h"
#include "../../core/Network.h"

namespace TrafficAssignment {

  /**
   * @brief TAPAS implementation based on the TAsK framework (Perederieieva et al.).
   *
   * Bush-based TAPAS algorithm with Paired Alternative Segments (PAS).
   * Features:
   * - Halley step with Armijo backtracking for flow shifts
   * - Stall detection with bush restart
   * - Random origin shuffle per iteration
   * - Sequential PAS processing (no priority queue)
   * - Post-iteration deleteUnusedPASAndMoveFlow pass
   * - nb_flow_moves lifecycle counter for aggressive PAS pruning
   *
   * @tparam T Numeric type for flow computations.
   */
  template <typename T>
  class TapasApproach : public TrafficAssignmentApproach<T> {
  public:
    TapasApproach(Network<T>& network,
                  T alpha = 1e-6,
                  int max_iterations = 200,
                  T mu = T(0.5),
                  T v = T(0.25))
      : TrafficAssignmentApproach<T>::TrafficAssignmentApproach(network, alpha),
        max_iterations_(max_iterations),
        mu_(mu),
        v_(v) {

      for (int i = 0; i < this->network_.number_of_zones(); i++) {
        if (this->network_.origin_info()[i].size() > 0) {
          origins_.push_back(i);
        }
      }
      number_of_origins_ = static_cast<int>(origins_.size());

      bush_flows_.resize(number_of_origins_);
      sp_parent_.resize(number_of_origins_);
      sp_distance_.resize(number_of_origins_);
      sp_tree_links_.resize(number_of_origins_,
          std::vector<bool>(this->network_.number_of_links(), false));

      for (int o = 0; o < number_of_origins_; o++) {
        bush_flows_[o].resize(this->network_.number_of_links(), T(0));
        sp_parent_[o].resize(this->network_.number_of_nodes(), -1);
        sp_distance_[o].resize(this->network_.number_of_nodes(), T(0));
      }
    }

    ~TapasApproach() = default;

    void Reset() override {
      TrafficAssignmentApproach<T>::Reset();
      ClearInternalState();
      has_valid_state_ = false;
    }

    void ComputeTrafficFlows(bool statistics_recording = false) override {
      if (!has_valid_state_) {
        this->network_.Reset();
        ClearInternalState();
        AllOrNothingAssignment();
      } else {
        ClearPASState();
      }

      if (statistics_recording) {
        this->statistics_recorder_.StartRecording(this->GetApproachName());
        this->statistics_recorder_.RecordStatistics();
      }

      T best_rgap = T(1);
      int stall_count = 0;
      const int kRgapCheckInterval = 5;

      for (int i = 0; i < max_iterations_; i++) {
        EquilibrationIteration(i + 1);

        bool should_check = (i == 0) || ((i + 1) % kRgapCheckInterval == 0);
        if (!should_check && !statistics_recording) continue;

        T rgap = this->network_.RelativeGap();
        if (statistics_recording) this->statistics_recorder_.RecordStatistics();
        if (std::abs(rgap) <= this->alpha_) break;

        if (should_check) {
          T abs_rgap = std::abs(rgap);
          if (abs_rgap < best_rgap * T(0.999)) {
            best_rgap = abs_rgap;
            stall_count = 0;
          } else {
            ++stall_count;
          }
          if (stall_count >= 3) {
            if (best_rgap < this->alpha_ * T(1000)) break;  // precision floor
            RestartBushFlows();
            ClearPASState();
            stall_count = 0;
          }
        }
      }

      has_valid_state_ = true;
    }

    std::string GetApproachName() override {
      return "TapasNewtonStep";
    }

  private:
    static constexpr T kZeroFlow = T(1e-15);
    static constexpr T kDirTol = T(1e-15);

    int max_iterations_;
    T mu_;
    T v_;
    int number_of_origins_;
    std::vector<int> origins_;
    bool has_valid_state_ = false;
    std::mt19937 rng_{42};

    // Per-origin bush data
    std::vector<std::vector<T>> bush_flows_;          // [origin][link]
    std::vector<std::vector<int>> sp_parent_;         // [origin][node] = parent node
    std::vector<std::vector<T>> sp_distance_;         // [origin][node] = SP distance
    std::vector<std::vector<bool>> sp_tree_links_;    // [origin][link] = in SP tree?

    // PAS data
    struct PAS {
      std::vector<int> cheap_segment;
      std::vector<int> expensive_segment;
      std::unordered_set<int> origins;  // origin indices using this PAS
      int nb_flow_moves = 0;
    };
    std::vector<PAS> pas_list_;

    // Adaptive threshold counter (TAsK: nbIter_)
    int pas_creation_iteration_ = 1;

    /// @brief Clears all internal state for a cold start.
    void ClearInternalState() {
      for (int o = 0; o < number_of_origins_; o++) {
        std::fill(bush_flows_[o].begin(), bush_flows_[o].end(), T(0));
        std::fill(sp_parent_[o].begin(), sp_parent_[o].end(), -1);
        std::fill(sp_distance_[o].begin(), sp_distance_[o].end(), T(0));
        std::fill(sp_tree_links_[o].begin(), sp_tree_links_[o].end(), false);
      }
      ClearPASState();
    }

    /// @brief Clears PAS data structures (for warm start or cold start).
    void ClearPASState() {
      pas_list_.clear();
      pas_creation_iteration_ = 1;
    }

    /// @brief All-or-nothing initialization: assigns full demand to shortest paths.
    void AllOrNothingAssignment() {
      for (int origin_index = 0; origin_index < number_of_origins_; origin_index++) {
        auto shortest_routes = this->network_.ComputeSingleOriginBestRoutes(origins_[origin_index]);
        this->network_.AddRoutes(shortest_routes);
        for (std::size_t r = 0; r < shortest_routes.size(); r++) {
          T demand = this->network_.od_pairs()[shortest_routes[r].first].GetDemand();
          for (int link_id : shortest_routes[r].second) {
            bush_flows_[origin_index][link_id] += demand;
          }
        }
      }
    }

    // =========================================================================
    // Cycle removal (DFS)
    // =========================================================================

    void RemoveCycleFlow(int origin_index, const std::vector<int>& cycle) {
      T min_flow = bush_flows_[origin_index][cycle[0]];
      for (int link_id : cycle) {
        min_flow = std::min(min_flow, bush_flows_[origin_index][link_id]);
      }
      for (int link_id : cycle) {
        if (bush_flows_[origin_index][link_id] - min_flow < kZeroFlow)
          bush_flows_[origin_index][link_id] = T(0);
        else
          bush_flows_[origin_index][link_id] -= min_flow;
        T current_flow = this->network_.links()[link_id].flow;
        if (current_flow - min_flow < kZeroFlow)
          this->network_.SetLinkFlow(link_id, -current_flow);
        else
          this->network_.SetLinkFlow(link_id, -min_flow);
      }
    }

    void DfsForCycleIdentification(int origin_index, int cur_node,
                                   std::vector<int>& link_used,
                                   std::vector<int>& state,
                                   bool& cycle_found) {
      if (cycle_found || state[cur_node] != 0) return;
      state[cur_node] = 1;
      for (int now : this->network_.adjacency()[cur_node]) {
        if (cycle_found) return;
        int next_node = this->network_.links()[now].term;
        if (bush_flows_[origin_index][now] < kZeroFlow || state[next_node] == 2) continue;
        if (state[next_node] == 1) {
          cycle_found = true;
          std::vector<int> cycle;
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

    /// @brief Collect all nodes that belong to the bush (have positive bush flow on an adjacent link).
    std::vector<int> GetBushNodes(int origin_index) const {
      std::vector<bool> is_bush_node(this->network_.number_of_nodes(), false);
      for (int link_id = 0; link_id < this->network_.number_of_links(); link_id++) {
        if (bush_flows_[origin_index][link_id] >= kZeroFlow) {
          is_bush_node[this->network_.links()[link_id].init] = true;
          is_bush_node[this->network_.links()[link_id].term] = true;
        }
      }
      std::vector<int> nodes;
      for (int i = 0; i < this->network_.number_of_nodes(); i++)
        if (is_bush_node[i]) nodes.push_back(i);
      return nodes;
    }

    void RemoveCyclicFlows(int origin_index) {
      bool cycle_found = true;
      while (cycle_found) {
        cycle_found = false;
        auto bush_nodes = GetBushNodes(origin_index);
        std::vector<int> link_used(this->network_.number_of_nodes(), -1);
        std::vector<int> state(this->network_.number_of_nodes(), 0);
        for (int node : bush_nodes) {
          if (state[node] == 0) {
            DfsForCycleIdentification(origin_index, node, link_used, state, cycle_found);
            if (cycle_found) break;
          }
        }
      }
    }

    // =========================================================================
    // Dijkstra shortest path tree
    // =========================================================================

    void BuildDijkstraTree(int origin_index) {
      std::fill(sp_tree_links_[origin_index].begin(), sp_tree_links_[origin_index].end(), false);
      std::fill(sp_parent_[origin_index].begin(), sp_parent_[origin_index].end(), -1);

      int origin = origins_[origin_index];
      std::priority_queue<std::pair<T,int>, std::vector<std::pair<T,int>>,
                          std::greater<std::pair<T,int>>> pq;
      pq.push({T(0), -1});  // (cost, link_index) — -1 means origin node
      std::vector<bool> processed(this->network_.number_of_nodes(), false);

      while (!pq.empty()) {
        auto [cur_delay, link_index] = pq.top();
        pq.pop();
        int u = (link_index == -1) ? origin : this->network_.links()[link_index].term;
        if (processed[u]) continue;
        processed[u] = true;
        sp_distance_[origin_index][u] = cur_delay;
        if (u != origin) {
          sp_parent_[origin_index][u] = this->network_.links()[link_index].init;
          sp_tree_links_[origin_index][link_index] = true;
        }
        for (int adj : this->network_.adjacency()[u]) {
          int v = this->network_.links()[adj].term;
          if (!processed[v]) {
            pq.push({cur_delay + this->network_.links()[adj].Delay(), adj});
          }
        }
      }
    }

    // =========================================================================
    // PAS cost recalculation
    // =========================================================================

    /// @brief Recalculate cheap/expensive labels for all PAS, swapping if reversed.
    void RecalculateAllPASCosts() {
      for (auto& pas : pas_list_) {
        T cheap_cost = Link<T>::GetLinksDelay(this->network_.links(), pas.cheap_segment);
        T exp_cost = Link<T>::GetLinksDelay(this->network_.links(), pas.expensive_segment);
        if (cheap_cost > exp_cost) {
          std::swap(pas.cheap_segment, pas.expensive_segment);
        }
      }
    }

    // =========================================================================
    // PAS creation (TAsK: PASManager::createNewPAS)
    // =========================================================================

    /// @brief Linear scan for existing PAS matching cheap/expensive last links (TAsK: PASManager::pasExist).
    int FindExistingPAS(int sp_link, int exp_link) const {
      for (std::size_t i = 0; i < pas_list_.size(); i++) {
        if (!pas_list_[i].cheap_segment.empty() && !pas_list_[i].expensive_segment.empty()
            && pas_list_[i].cheap_segment.back() == sp_link
            && pas_list_[i].expensive_segment.back() == exp_link) {
          return static_cast<int>(i);
        }
      }
      return -1;
    }

    /// @brief Adaptive PAS creation threshold: 10 * 10^(-iteration).
    T CalcPasCreationThreshold() const {
      return T(10.0) * std::pow(T(10.0), -T(pas_creation_iteration_));
    }

    /// @brief Backward BFS to find PAS: follows links with bush_flow > flow_threshold.
    /// Returns (cheap_segment, expensive_segment) or empty if not found.
    std::pair<std::vector<int>, std::vector<int>>
    ConstructPASByBFS(int origin_index, int link_index, T flow_threshold) {
      int term_node = this->network_.links()[link_index].term;
      int init_node = this->network_.links()[link_index].init;

      // Build set of nodes on SP path from origin to term_node
      std::unordered_set<int> sp_path_nodes;
      {
        int node = term_node;
        while (node != origins_[origin_index]) {
          int parent = sp_parent_[origin_index][node];
          if (parent < 0) return {{}, {}};
          node = parent;
          sp_path_nodes.insert(node);
        }
      }

      // BFS backward from init_node along links with sufficient bush flow
      std::queue<int> q;
      std::unordered_set<int> visited;
      std::unordered_map<int, int> used_link;  // node -> link that brought us here

      used_link[init_node] = link_index;
      q.push(init_node);
      visited.insert(init_node);

      int junction_node = -1;
      while (!q.empty()) {
        int node = q.front();
        q.pop();

        if (sp_path_nodes.count(node)) {
          junction_node = node;
          break;
        }

        for (int rev_link : this->network_.reverse_adjacency()[node]) {
          int prev_node = this->network_.links()[rev_link].init;
          if (!visited.count(prev_node) && bush_flows_[origin_index][rev_link] >= flow_threshold) {
            visited.insert(prev_node);
            q.push(prev_node);
            used_link[prev_node] = rev_link;
            // Early termination: check immediately when adding neighbor (TAsK behavior)
            if (sp_path_nodes.count(prev_node)) {
              junction_node = prev_node;
              break;
            }
          }
        }
        if (junction_node >= 0) break;
      }

      if (junction_node < 0) return {{}, {}};

      // Reconstruct expensive segment: junction -> ... -> term_node (via BFS links)
      std::vector<int> expensive_segment;
      {
        int node = junction_node;
        while (node != term_node) {
          int link = used_link.at(node);
          expensive_segment.push_back(link);
          node = this->network_.links()[link].term;
        }
      }

      // Reconstruct cheap segment: junction -> ... -> term_node (via SP tree)
      std::vector<int> cheap_segment;
      {
        int node = term_node;
        while (node != junction_node) {
          int parent = sp_parent_[origin_index][node];
          if (parent < 0) return {{}, {}};
          // Find the link from parent to node
          int found_link = -1;
          for (int adj : this->network_.adjacency()[parent]) {
            if (this->network_.links()[adj].term == node) {
              found_link = adj;
              break;
            }
          }
          if (found_link < 0) return {{}, {}};
          cheap_segment.push_back(found_link);
          node = parent;
        }
        std::reverse(cheap_segment.begin(), cheap_segment.end());
      }

      return {cheap_segment, expensive_segment};
    }

    /// @brief Check if PAS is cost-effective: cost_diff >= mu * reduced_cost.
    bool CheckCostEffective(std::size_t pas_idx, T mu_reduced_cost) const {
      const auto& pas = pas_list_[pas_idx];
      T cost_diff = std::abs(
          Link<T>::GetLinksDelay(this->network_.links(), pas.expensive_segment)
          - Link<T>::GetLinksDelay(this->network_.links(), pas.cheap_segment));
      return cost_diff >= mu_reduced_cost;
    }

    /// @brief Check if PAS is flow-effective: min origin flow on expensive segment >= v * exp_link_flow.
    bool CheckFlowEffective(std::size_t pas_idx, int origin_index, int exp_link_index) const {
      const auto& pas = pas_list_[pas_idx];
      T min_flow = T(-1);
      for (int link_id : pas.expensive_segment) {
        T f = bush_flows_[origin_index][link_id];
        if (min_flow < T(0) || f < min_flow) min_flow = f;
      }
      if (min_flow < T(0)) return false;
      return min_flow >= v_ * bush_flows_[origin_index][exp_link_index];
    }

    bool CheckEffective(std::size_t pas_idx, int origin_index, int exp_link_index, T mu_rc) const {
      return CheckCostEffective(pas_idx, mu_rc) && CheckFlowEffective(pas_idx, origin_index, exp_link_index);
    }

    /// @brief Add origin to existing PAS if it has flow on the expensive segment.
    bool TryAddOriginToExistingPAS(std::size_t pas_idx, int origin_index) {
      auto& pas = pas_list_[pas_idx];
      if (pas.origins.count(origin_index)) return true;

      // Verify origin has flow on expensive segment
      for (int link_id : pas.expensive_segment) {
        if (bush_flows_[origin_index][link_id] <= T(0)) return false;
      }
      pas.origins.insert(origin_index);
      return true;
    }

    /// @brief Creates new PAS for the given origin and non-tree link.
    /// Implements TAsK PASManager::createNewPAS two-pass logic.
    void CreateNewPAS(int origin_index, int link_index) {
      // Compute reduced cost
      int from_node = this->network_.links()[link_index].init;
      int to_node = this->network_.links()[link_index].term;
      T reduced_cost = sp_distance_[origin_index][from_node]
                      + this->network_.links()[link_index].Delay()
                      - sp_distance_[origin_index][to_node];
      if (reduced_cost < kDirTol) return;

      // Find SP link at terminal node
      int parent_node = sp_parent_[origin_index][to_node];
      if (parent_node < 0) return;
      int sp_link = -1;
      for (int adj : this->network_.adjacency()[parent_node]) {
        if (this->network_.links()[adj].term == to_node) {
          sp_link = adj;
          break;
        }
      }
      if (sp_link < 0) return;

      T mu_reduced_cost = mu_ * reduced_cost;

      // Linear scan for existing PAS with matching last links (TAsK: PASManager::pasExist)
      int existing_idx = FindExistingPAS(sp_link, link_index);
      if (existing_idx >= 0) {
        if (TryAddOriginToExistingPAS(static_cast<std::size_t>(existing_idx), origin_index)) {
          if (CheckEffective(static_cast<std::size_t>(existing_idx), origin_index, link_index, mu_reduced_cost)) {
            return;  // Reused existing effective PAS
          }
        }
      }

      // Pass 1: Relaxed BFS (any positive flow)
      bool is_effective = false;
      {
        auto [cheap, expensive] = ConstructPASByBFS(origin_index, link_index, kZeroFlow);
        if (!cheap.empty() && !expensive.empty()) {
          T cost_diff = std::abs(
              Link<T>::GetLinksDelay(this->network_.links(), expensive)
              - Link<T>::GetLinksDelay(this->network_.links(), cheap));
          if (cost_diff >= kDirTol) {
            std::size_t idx = pas_list_.size();
            PAS new_pas;
            new_pas.cheap_segment = std::move(cheap);
            new_pas.expensive_segment = std::move(expensive);
            new_pas.origins.insert(origin_index);
            new_pas.nb_flow_moves = 0;
            pas_list_.push_back(std::move(new_pas));
            is_effective = CheckEffective(idx, origin_index, link_index, mu_reduced_cost);
          }
        }
      }

      // Pass 2: Strict BFS if not effective and reduced cost > adaptive threshold
      if (!is_effective && reduced_cost > CalcPasCreationThreshold()) {
        T flow_threshold = v_ * bush_flows_[origin_index][link_index];
        auto [cheap, expensive] = ConstructPASByBFS(origin_index, link_index, flow_threshold);
        if (!cheap.empty() && !expensive.empty()) {
          T cost_diff = std::abs(
              Link<T>::GetLinksDelay(this->network_.links(), expensive)
              - Link<T>::GetLinksDelay(this->network_.links(), cheap));
          if (cost_diff >= kDirTol) {
            PAS new_pas;
            new_pas.cheap_segment = std::move(cheap);
            new_pas.expensive_segment = std::move(expensive);
            new_pas.origins.insert(origin_index);
            new_pas.nb_flow_moves = 0;
            pas_list_.push_back(std::move(new_pas));
          }
        }
      }
    }

    // =========================================================================
    // MoveFlow — Halley step with Armijo backtracking
    // =========================================================================

    /// @brief Performs Halley flow shift on a PAS with Armijo backtracking.
    /// Returns true if flow was actually moved.
    bool MoveFlow(std::size_t pas_idx) {
      auto& pas = pas_list_[pas_idx];

      // Decrement flow move counter
      --pas.nb_flow_moves;

      // Recalculate segment costs, swap if direction reversed
      T cheap_cost = Link<T>::GetLinksDelay(this->network_.links(), pas.cheap_segment);
      T exp_cost = Link<T>::GetLinksDelay(this->network_.links(), pas.expensive_segment);
      if (cheap_cost > exp_cost) {
        std::swap(pas.cheap_segment, pas.expensive_segment);
        std::swap(cheap_cost, exp_cost);
      }

      T cost_diff = exp_cost - cheap_cost;
      if (cost_diff < kDirTol) return false;

      // Newton step: dFlow = cost_diff / sum(DelayDer on both segments)
      T der_sum = Link<T>::GetLinksDelayDer(this->network_.links(), pas.cheap_segment)
                + Link<T>::GetLinksDelayDer(this->network_.links(), pas.expensive_segment);
      // When der_sum == 0 (constant-delay links with b=0), full shift is correct:
      // the only equilibrium puts all flow on the cheapest path.
      T dFlow = (der_sum != T(0)) ? (cost_diff / der_sum) : std::numeric_limits<T>::infinity();

      // Halley correction using second derivatives
      if (der_sum != T(0)) {
        T from_2nd = Link<T>::GetLinksDelaySecondDer(this->network_.links(), pas.expensive_segment);
        T to_2nd = Link<T>::GetLinksDelaySecondDer(this->network_.links(), pas.cheap_segment);
        T halley = T(1) - dFlow * (from_2nd - to_2nd) / (T(2) * der_sum);
        if (halley > T(0.1)) dFlow /= halley;
      }

      // Compute per-origin available flow on expensive segment (totalShift)
      T totalShift = T(0);
      std::unordered_map<int, T> origin_min_flow;
      for (int origin_index : pas.origins) {
        T min_flow = std::numeric_limits<T>::max();
        for (int link_id : pas.expensive_segment) {
          min_flow = std::min(min_flow, bush_flows_[origin_index][link_id]);
        }
        origin_min_flow[origin_index] = min_flow;
        totalShift += min_flow;
      }

      if (totalShift <= kZeroFlow) return false;

      // Clamp Newton step
      dFlow = std::min(dFlow, totalShift);
      if (dFlow <= kZeroFlow) return false;

      // Armijo backtracking: halve step if it overshoots past equilibrium
      for (int i = 0; i < 4 && dFlow > kZeroFlow; ++i) {
        T cf = T(0), ct = T(0);
        for (int l : pas.expensive_segment) cf += this->network_.links()[l].Delay(this->network_.links()[l].flow - dFlow);
        for (int l : pas.cheap_segment)     ct += this->network_.links()[l].Delay(this->network_.links()[l].flow + dFlow);
        if (ct <= cf) break;
        dFlow /= T(2);
      }
      if (dFlow <= kZeroFlow) return false;

      // Distribute proportionally across origins
      for (int origin_index : pas.origins) {
        T origin_shift = (origin_min_flow[origin_index] / totalShift) * dFlow;
        if (origin_shift <= kZeroFlow) continue;

        for (int link_id : pas.expensive_segment) {
          bush_flows_[origin_index][link_id] -= origin_shift;
          if (bush_flows_[origin_index][link_id] < kZeroFlow) {
            bush_flows_[origin_index][link_id] = T(0);
          }
        }
        for (int link_id : pas.cheap_segment) {
          bush_flows_[origin_index][link_id] += origin_shift;
        }
      }

      // Update aggregate link flows
      for (int link_id : pas.expensive_segment) {
        T current = this->network_.links()[link_id].flow;
        if (current - dFlow < kZeroFlow) {
          this->network_.SetLinkFlow(link_id, -current);
        } else {
          this->network_.SetLinkFlow(link_id, -dFlow);
        }
      }
      for (int link_id : pas.cheap_segment) {
        this->network_.SetLinkFlow(link_id, dFlow);
      }

      // Flow was moved — increment counter
      ++pas.nb_flow_moves;
      return true;
    }

    // =========================================================================
    // MoveFlowOnAllPAS — shift ALL active PAS
    // =========================================================================

    void MoveFlowOnAllPAS() {
      for (std::size_t i = 0; i < pas_list_.size(); i++) {
        MoveFlow(i);
      }
    }

    // =========================================================================
    // DeleteUnusedPASAndMoveFlow (TAsK: PASManager::deleteUnusedPASAndMoveFlow)
    // =========================================================================

    void DeleteUnusedPASAndMoveFlow() {
      ++pas_creation_iteration_;  // TAsK increments at start of this method
      std::size_t write = 0;
      for (std::size_t i = 0; i < pas_list_.size(); i++) {
        MoveFlow(i);
        // TAsK isUnused(): read counter, reset to 0, return retVal < 0
        int retVal = pas_list_[i].nb_flow_moves;
        pas_list_[i].nb_flow_moves = 0;
        if (retVal < 0) continue;  // unused, skip
        if (write != i) pas_list_[write] = std::move(pas_list_[i]);
        write++;
      }
      pas_list_.resize(write);
    }

    // =========================================================================
    // Bush restart (stall recovery)
    // =========================================================================

    /// @brief Restarts bush flows: subtract from network, recompute shortest paths, reassign on shortest paths.
    void RestartBushFlows() {
      for (int origin_index = 0; origin_index < number_of_origins_; origin_index++) {
        // Subtract this origin's bush flows from the network
        for (int link_id = 0; link_id < this->network_.number_of_links(); link_id++) {
          if (bush_flows_[origin_index][link_id] > T(0)) {
            this->network_.SetLinkFlow(link_id, -bush_flows_[origin_index][link_id]);
            bush_flows_[origin_index][link_id] = T(0);
          }
        }

        // Recompute shortest paths with current costs (without this origin's flow)
        auto shortest_routes = this->network_.ComputeSingleOriginBestRoutes(origins_[origin_index]);

        // Reassign demand on shortest paths
        for (std::size_t r = 0; r < shortest_routes.size(); r++) {
          T demand = this->network_.od_pairs()[shortest_routes[r].first].GetDemand();
          for (int link_id : shortest_routes[r].second) {
            bush_flows_[origin_index][link_id] += demand;
            this->network_.SetLinkFlow(link_id, demand);
          }
        }
      }
    }

    // =========================================================================
    // Equilibration iteration
    // =========================================================================

    void EquilibrationIteration(int iteration_number) {
      // Random origin shuffle
      std::vector<int> order(number_of_origins_);
      std::iota(order.begin(), order.end(), 0);
      std::shuffle(order.begin(), order.end(), rng_);

      for (int origin_index : order) {
        // 1. Remove cyclic flows
        RemoveCyclicFlows(origin_index);

        // 2. Build Dijkstra SP tree
        BuildDijkstraTree(origin_index);

        // 3. Recalculate all PAS costs (swap cheap/expensive if needed)
        RecalculateAllPASCosts();

        // 4. Create new PAS — node-first iteration (matching TAsK)
        auto bush_nodes = GetBushNodes(origin_index);
        for (int node : bush_nodes) {
          int sp_par = sp_parent_[origin_index][node];
          if (sp_par < 0) continue;  // origin or unreachable — no SP incoming link
          // Find SP tree link incoming to this node
          int sp_link_at_node = -1;
          for (int adj : this->network_.adjacency()[sp_par]) {
            if (this->network_.links()[adj].term == node) {
              sp_link_at_node = adj;
              break;
            }
          }
          // Process non-SP incoming links with positive bush flow
          for (int rev_link : this->network_.reverse_adjacency()[node]) {
            if (rev_link == sp_link_at_node) continue;
            if (bush_flows_[origin_index][rev_link] < kZeroFlow) continue;
            CreateNewPAS(origin_index, rev_link);
          }
        }

        // 5. Move flow on ALL active PAS (TAsK: moveFlow on all global PAS)
        MoveFlowOnAllPAS();
      }

      // 6. Post-iteration: deleteUnusedPASAndMoveFlow (also increments pas_creation_iteration_)
      DeleteUnusedPASAndMoveFlow();
    }
  };

}  // namespace TrafficAssignment

#endif

#ifndef TAPAS_APPROACH_H
#define TAPAS_APPROACH_H

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <queue>
#include "../TrafficAssignmentApproach.h"
#include "../../core/Network.h"
#include "BushData.h"
#include "CycleRemoval.h"
#include "PasManager.h"

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
        v_(v),
        pas_manager_(network, mu, v, kZeroFlow, kDirTol, bushes_, origins_) {

      for (int i = 0; i < this->network_.number_of_zones(); i++) {
        if (this->network_.origin_info()[i].size() > 0) {
          origins_.push_back(i);
        }
      }
      number_of_origins_ = static_cast<int>(origins_.size());

      bushes_.resize(number_of_origins_);
      for (int o = 0; o < number_of_origins_; o++) {
        bushes_[o].Initialize(this->network_.number_of_links(), this->network_.number_of_nodes());
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
        pas_manager_.Clear();
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
            pas_manager_.Clear();
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
    int number_of_origins_ = 0;
    std::vector<int> origins_;
    bool has_valid_state_ = false;
    std::mt19937 rng_{42};

    std::vector<OriginBush<T>> bushes_;
    PasManager<T> pas_manager_;

    void ClearInternalState() {
      for (int o = 0; o < number_of_origins_; o++) {
        bushes_[o].Clear();
      }
      pas_manager_.Clear();
    }

    void AllOrNothingAssignment() {
      for (int origin_index = 0; origin_index < number_of_origins_; origin_index++) {
        auto shortest_routes = this->network_.ComputeSingleOriginBestRoutes(origins_[origin_index]);
        this->network_.AddRoutes(shortest_routes);
        for (std::size_t r = 0; r < shortest_routes.size(); r++) {
          T demand = this->network_.od_pairs()[shortest_routes[r].first].GetDemand();
          for (int link_id : shortest_routes[r].second) {
            bushes_[origin_index].flows[link_id] += demand;
          }
        }
      }
    }

    // =========================================================================
    // Dijkstra shortest path tree
    // =========================================================================

    void BuildDijkstraTree(int origin_index) {
      auto& bush = bushes_[origin_index];
      std::fill(bush.sp_tree_links.begin(), bush.sp_tree_links.end(), false);
      std::fill(bush.sp_parent.begin(), bush.sp_parent.end(), -1);

      int origin = origins_[origin_index];
      std::priority_queue<std::pair<T,int>, std::vector<std::pair<T,int>>,
                          std::greater<std::pair<T,int>>> pq;
      pq.push({T(0), -1});
      std::vector<bool> processed(this->network_.number_of_nodes(), false);

      while (!pq.empty()) {
        auto [cur_delay, link_index] = pq.top();
        pq.pop();
        int u = (link_index == -1) ? origin : this->network_.links()[link_index].term;
        if (processed[u]) continue;
        processed[u] = true;
        bush.sp_distance[u] = cur_delay;
        if (u != origin) {
          bush.sp_parent[u] = this->network_.links()[link_index].init;
          bush.sp_tree_links[link_index] = true;
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
    // Bush restart (stall recovery)
    // =========================================================================

    void RestartBushFlows() {
      for (int origin_index = 0; origin_index < number_of_origins_; origin_index++) {
        auto& bush = bushes_[origin_index];
        for (int link_id = 0; link_id < this->network_.number_of_links(); link_id++) {
          if (bush.flows[link_id] > T(0)) {
            this->network_.SetLinkFlow(link_id, -bush.flows[link_id]);
            bush.flows[link_id] = T(0);
          }
        }

        auto shortest_routes = this->network_.ComputeSingleOriginBestRoutes(origins_[origin_index]);

        for (std::size_t r = 0; r < shortest_routes.size(); r++) {
          T demand = this->network_.od_pairs()[shortest_routes[r].first].GetDemand();
          for (int link_id : shortest_routes[r].second) {
            bush.flows[link_id] += demand;
            this->network_.SetLinkFlow(link_id, demand);
          }
        }
      }
    }

    // =========================================================================
    // Equilibration iteration
    // =========================================================================

    void EquilibrationIteration(int iteration_number) {
      std::vector<int> order(number_of_origins_);
      std::iota(order.begin(), order.end(), 0);
      std::shuffle(order.begin(), order.end(), rng_);

      for (int origin_index : order) {
        // 1. Remove cyclic flows
        RemoveCyclicFlows(bushes_[origin_index], this->network_, kZeroFlow);

        // 2. Build Dijkstra SP tree
        BuildDijkstraTree(origin_index);

        // 3. Recalculate all PAS costs (swap cheap/expensive if needed)
        pas_manager_.RecalculateAllPASCosts();

        // 4. Create new PAS — node-first iteration (matching TAsK)
        auto bush_nodes = GetBushNodes(bushes_[origin_index], this->network_, kZeroFlow);
        for (int node : bush_nodes) {
          int sp_par = bushes_[origin_index].sp_parent[node];
          if (sp_par < 0) continue;
          int sp_link_at_node = -1;
          for (int adj : this->network_.adjacency()[sp_par]) {
            if (this->network_.links()[adj].term == node) {
              sp_link_at_node = adj;
              break;
            }
          }
          for (int rev_link : this->network_.reverse_adjacency()[node]) {
            if (rev_link == sp_link_at_node) continue;
            if (bushes_[origin_index].flows[rev_link] < kZeroFlow) continue;
            pas_manager_.CreateNewPAS(origin_index, rev_link);
          }
        }

        // 5. Move flow on ALL active PAS
        pas_manager_.MoveFlowOnAllPAS();
      }

      // 6. Post-iteration cleanup
      pas_manager_.DeleteUnusedPASAndMoveFlow();
    }
  };

}  // namespace TrafficAssignment

#endif

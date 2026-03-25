#ifndef PAS_MANAGER_H
#define PAS_MANAGER_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <algorithm>
#include <cmath>
#include <limits>
#include "../../core/Network.h"
#include "BushData.h"

namespace TrafficAssignment {

  template <typename T>
  class PasManager {
  public:
    struct PAS {
      std::vector<int> cheap_segment;
      std::vector<int> expensive_segment;
      std::unordered_set<int> origins;
      int nb_flow_moves = 0;
    };

    PasManager(Network<T>& network, T mu, T v, T zero_flow, T dir_tol,
               std::vector<OriginBush<T>>& bushes,
               const std::vector<int>& origins)
      : network_(network), mu_(mu), v_(v), kZeroFlow_(zero_flow), kDirTol_(dir_tol),
        bushes_(bushes), origins_(origins) {}

    void Clear() {
      pas_list_.clear();
      pas_creation_iteration_ = 1;
      pas_lookup_.clear();
    }

    void RecalculateAllPASCosts() {
      for (auto& pas : pas_list_) {
        T cheap_cost = Link<T>::GetLinksDelay(network_.links(), pas.cheap_segment);
        T exp_cost = Link<T>::GetLinksDelay(network_.links(), pas.expensive_segment);
        if (cheap_cost > exp_cost) {
          std::swap(pas.cheap_segment, pas.expensive_segment);
        }
      }
    }

    int FindExistingPAS(int sp_link, int exp_link) const {
      auto it = pas_lookup_.find({sp_link, exp_link});
      if (it != pas_lookup_.end()) {
        return it->second;
      }
      return -1;
    }

    T CalcPasCreationThreshold() const {
      return T(10.0) * std::pow(T(10.0), -T(pas_creation_iteration_));
    }

    std::pair<std::vector<int>, std::vector<int>>
    ConstructPASByBFS(int origin_index, int link_index, T flow_threshold) {
      int term_node = network_.links()[link_index].term;
      int init_node = network_.links()[link_index].init;

      std::unordered_set<int> sp_path_nodes;
      {
        int node = term_node;
        while (node != origins_[origin_index]) {
          int parent = bushes_[origin_index].sp_parent[node];
          if (parent < 0) return {{}, {}};
          node = parent;
          sp_path_nodes.insert(node);
        }
      }

      std::queue<int> q;
      std::unordered_set<int> visited;
      std::unordered_map<int, int> used_link;

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

        for (int rev_link : network_.reverse_adjacency()[node]) {
          int prev_node = network_.links()[rev_link].init;
          if (!visited.count(prev_node) && bushes_[origin_index].flows[rev_link] >= flow_threshold) {
            visited.insert(prev_node);
            q.push(prev_node);
            used_link[prev_node] = rev_link;
            if (sp_path_nodes.count(prev_node)) {
              junction_node = prev_node;
              break;
            }
          }
        }
        if (junction_node >= 0) break;
      }

      if (junction_node < 0) return {{}, {}};

      std::vector<int> expensive_segment;
      {
        int node = junction_node;
        while (node != term_node) {
          int link = used_link.at(node);
          expensive_segment.push_back(link);
          node = network_.links()[link].term;
        }
      }

      std::vector<int> cheap_segment;
      {
        int node = term_node;
        while (node != junction_node) {
          int parent = bushes_[origin_index].sp_parent[node];
          if (parent < 0) return {{}, {}};
          int found_link = -1;
          for (int adj : network_.adjacency()[parent]) {
            if (network_.links()[adj].term == node) {
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

    bool CheckCostEffective(std::size_t pas_idx, T mu_reduced_cost) const {
      const auto& pas = pas_list_[pas_idx];
      T cost_diff = std::abs(
          Link<T>::GetLinksDelay(network_.links(), pas.expensive_segment)
          - Link<T>::GetLinksDelay(network_.links(), pas.cheap_segment));
      return cost_diff >= mu_reduced_cost;
    }

    bool CheckFlowEffective(std::size_t pas_idx, int origin_index, int exp_link_index) const {
      const auto& pas = pas_list_[pas_idx];
      T min_flow = T(-1);
      for (int link_id : pas.expensive_segment) {
        T f = bushes_[origin_index].flows[link_id];
        if (min_flow < T(0) || f < min_flow) min_flow = f;
      }
      if (min_flow < T(0)) return false;
      return min_flow >= v_ * bushes_[origin_index].flows[exp_link_index];
    }

    bool CheckEffective(std::size_t pas_idx, int origin_index, int exp_link_index, T mu_rc) const {
      return CheckCostEffective(pas_idx, mu_rc) && CheckFlowEffective(pas_idx, origin_index, exp_link_index);
    }

    bool TryAddOriginToExistingPAS(std::size_t pas_idx, int origin_index) {
      auto& pas = pas_list_[pas_idx];
      if (pas.origins.count(origin_index)) return true;

      for (int link_id : pas.expensive_segment) {
        if (bushes_[origin_index].flows[link_id] <= T(0)) return false;
      }
      pas.origins.insert(origin_index);
      return true;
    }

    void CreateNewPAS(int origin_index, int link_index) {
      int from_node = network_.links()[link_index].init;
      int to_node = network_.links()[link_index].term;
      T reduced_cost = bushes_[origin_index].sp_distance[from_node]
                      + network_.links()[link_index].Delay()
                      - bushes_[origin_index].sp_distance[to_node];
      if (reduced_cost < kDirTol_) return;

      int parent_node = bushes_[origin_index].sp_parent[to_node];
      if (parent_node < 0) return;
      int sp_link = -1;
      for (int adj : network_.adjacency()[parent_node]) {
        if (network_.links()[adj].term == to_node) {
          sp_link = adj;
          break;
        }
      }
      if (sp_link < 0) return;

      T mu_reduced_cost = mu_ * reduced_cost;

      int existing_idx = FindExistingPAS(sp_link, link_index);
      if (existing_idx >= 0) {
        if (TryAddOriginToExistingPAS(static_cast<std::size_t>(existing_idx), origin_index)) {
          if (CheckEffective(static_cast<std::size_t>(existing_idx), origin_index, link_index, mu_reduced_cost)) {
            return;
          }
        }
      }

      // Pass 1: Relaxed BFS (any positive flow)
      bool is_effective = false;
      {
        auto [cheap, expensive] = ConstructPASByBFS(origin_index, link_index, kZeroFlow_);
        if (!cheap.empty() && !expensive.empty()) {
          T cost_diff = std::abs(
              Link<T>::GetLinksDelay(network_.links(), expensive)
              - Link<T>::GetLinksDelay(network_.links(), cheap));
          if (cost_diff >= kDirTol_) {
            std::size_t idx = pas_list_.size();
            PAS new_pas;
            new_pas.cheap_segment = std::move(cheap);
            new_pas.expensive_segment = std::move(expensive);
            new_pas.origins.insert(origin_index);
            new_pas.nb_flow_moves = 0;
            pas_list_.push_back(std::move(new_pas));
            UpdateLookup(idx);
            is_effective = CheckEffective(idx, origin_index, link_index, mu_reduced_cost);
          }
        }
      }

      // Pass 2: Strict BFS if not effective and reduced cost > adaptive threshold
      if (!is_effective && reduced_cost > CalcPasCreationThreshold()) {
        T flow_threshold = v_ * bushes_[origin_index].flows[link_index];
        auto [cheap, expensive] = ConstructPASByBFS(origin_index, link_index, flow_threshold);
        if (!cheap.empty() && !expensive.empty()) {
          T cost_diff = std::abs(
              Link<T>::GetLinksDelay(network_.links(), expensive)
              - Link<T>::GetLinksDelay(network_.links(), cheap));
          if (cost_diff >= kDirTol_) {
            std::size_t idx = pas_list_.size();
            PAS new_pas;
            new_pas.cheap_segment = std::move(cheap);
            new_pas.expensive_segment = std::move(expensive);
            new_pas.origins.insert(origin_index);
            new_pas.nb_flow_moves = 0;
            pas_list_.push_back(std::move(new_pas));
            UpdateLookup(idx);
          }
        }
      }
    }

    bool MoveFlow(std::size_t pas_idx) {
      auto& pas = pas_list_[pas_idx];

      --pas.nb_flow_moves;

      T cheap_cost = Link<T>::GetLinksDelay(network_.links(), pas.cheap_segment);
      T exp_cost = Link<T>::GetLinksDelay(network_.links(), pas.expensive_segment);
      if (cheap_cost > exp_cost) {
        std::swap(pas.cheap_segment, pas.expensive_segment);
        std::swap(cheap_cost, exp_cost);
      }

      T cost_diff = exp_cost - cheap_cost;
      if (cost_diff < kDirTol_) return false;

      T der_sum = Link<T>::GetLinksDelayDer(network_.links(), pas.cheap_segment)
                + Link<T>::GetLinksDelayDer(network_.links(), pas.expensive_segment);
      T dFlow = (der_sum != T(0)) ? (cost_diff / der_sum) : std::numeric_limits<T>::infinity();

      if (der_sum != T(0)) {
        T from_2nd = Link<T>::GetLinksDelaySecondDer(network_.links(), pas.expensive_segment);
        T to_2nd = Link<T>::GetLinksDelaySecondDer(network_.links(), pas.cheap_segment);
        T halley = T(1) - dFlow * (from_2nd - to_2nd) / (T(2) * der_sum);
        if (halley > T(0.1)) dFlow /= halley;
      }

      T totalShift = T(0);
      std::unordered_map<int, T> origin_min_flow;
      for (int origin_index : pas.origins) {
        T min_flow = std::numeric_limits<T>::max();
        for (int link_id : pas.expensive_segment) {
          min_flow = std::min(min_flow, bushes_[origin_index].flows[link_id]);
        }
        origin_min_flow[origin_index] = min_flow;
        totalShift += min_flow;
      }

      if (totalShift <= kZeroFlow_) return false;

      dFlow = std::min(dFlow, totalShift);
      if (dFlow <= kZeroFlow_) return false;

      for (int i = 0; i < 4 && dFlow > kZeroFlow_; ++i) {
        T cf = T(0), ct = T(0);
        for (int l : pas.expensive_segment) cf += network_.links()[l].Delay(network_.links()[l].flow - dFlow);
        for (int l : pas.cheap_segment)     ct += network_.links()[l].Delay(network_.links()[l].flow + dFlow);
        if (ct <= cf) break;
        dFlow /= T(2);
      }
      if (dFlow <= kZeroFlow_) return false;

      for (int origin_index : pas.origins) {
        T origin_shift = (origin_min_flow[origin_index] / totalShift) * dFlow;
        if (origin_shift <= kZeroFlow_) continue;

        for (int link_id : pas.expensive_segment) {
          bushes_[origin_index].flows[link_id] -= origin_shift;
          if (bushes_[origin_index].flows[link_id] < kZeroFlow_) {
            bushes_[origin_index].flows[link_id] = T(0);
          }
        }
        for (int link_id : pas.cheap_segment) {
          bushes_[origin_index].flows[link_id] += origin_shift;
        }
      }

      for (int link_id : pas.expensive_segment) {
        T current = network_.links()[link_id].flow;
        if (current - dFlow < kZeroFlow_) {
          network_.SetLinkFlow(link_id, -current);
        } else {
          network_.SetLinkFlow(link_id, -dFlow);
        }
      }
      for (int link_id : pas.cheap_segment) {
        network_.SetLinkFlow(link_id, dFlow);
      }

      ++pas.nb_flow_moves;
      return true;
    }

    void MoveFlowOnAllPAS() {
      for (std::size_t i = 0; i < pas_list_.size(); i++) {
        MoveFlow(i);
      }
    }

    void DeleteUnusedPASAndMoveFlow() {
      ++pas_creation_iteration_;
      std::size_t write = 0;
      for (std::size_t i = 0; i < pas_list_.size(); i++) {
        MoveFlow(i);
        int retVal = pas_list_[i].nb_flow_moves;
        pas_list_[i].nb_flow_moves = 0;
        if (retVal < 0) continue;  // unused, skip
        if (write != i) pas_list_[write] = std::move(pas_list_[i]);
        write++;
      }
      pas_list_.resize(write);
      RebuildLookup();
    }

  private:
    Network<T>& network_;
    T mu_;
    T v_;
    T kZeroFlow_;
    T kDirTol_;
    std::vector<OriginBush<T>>& bushes_;
    const std::vector<int>& origins_;

    std::vector<PAS> pas_list_;
    int pas_creation_iteration_ = 1;

    // O(1) PAS lookup by (cheap_last_link, expensive_last_link)
    struct PairHash {
      std::size_t operator()(const std::pair<int,int>& p) const {
        return std::hash<long long>()(static_cast<long long>(p.first) * 1000003LL + p.second);
      }
    };
    std::unordered_map<std::pair<int,int>, int, PairHash> pas_lookup_;

    void UpdateLookup(std::size_t idx) {
      const auto& pas = pas_list_[idx];
      if (!pas.cheap_segment.empty() && !pas.expensive_segment.empty()) {
        pas_lookup_[{pas.cheap_segment.back(), pas.expensive_segment.back()}] = static_cast<int>(idx);
      }
    }

    void RebuildLookup() {
      pas_lookup_.clear();
      for (std::size_t i = 0; i < pas_list_.size(); i++) {
        UpdateLookup(i);
      }
    }
  };

}  // namespace TrafficAssignment

#endif  // PAS_MANAGER_H

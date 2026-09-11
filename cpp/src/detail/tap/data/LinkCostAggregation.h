#ifndef LINK_COST_AGGREGATION_H
#define LINK_COST_AGGREGATION_H

#include <vector>
#include <unordered_set>
#include <Eigen/Dense>
#include "Link.h"

namespace TrafficAssignment {

  template <typename T>
  T GetLinksDelay(const std::vector<Link<T>>& links, const std::vector<int>& links_list) {
    T ans = 0;
    for (auto now : links_list) {
      ans += links[now].Delay();
    }
    return ans;
  }

  template <typename T>
  T GetLinksDelayDer(const std::vector<Link<T>>& links, const std::vector<int>& links_list) {
    T ans = 0;
    for (auto now : links_list) {
      ans += links[now].DelayDer();
    }
    return ans;
  }

  template <typename T>
  T GetLinksDelaySecondDer(const std::vector<Link<T>>& links, const std::vector<int>& links_list) {
    T ans = 0;
    for (auto now : links_list) {
      ans += links[now].DelaySecondDer();
    }
    return ans;
  }

  template <typename T>
  bool CheckNonZeroLinksCapacity(const std::vector<Link<T>>& links, const std::vector<int>& links_list) {
    for (auto now : links_list) {
      if (links[now].capacity != 0) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  T CalculateJacobiElement(const std::vector<int>& route_i,
                           const std::vector<int>& route_j,
                           const std::vector<Link<T>>& links) {
    T sum = 0;
    std::unordered_set<int> links_i(route_i.begin(), route_i.end());
    for (int link_id : route_j) {
      if (links_i.count(link_id)) {
        sum += links[link_id].DelayDer();
      }
    }
    return sum;
  }

  template <typename T>
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> RoutesJacobiMatrix(
      const std::vector<std::vector<int>>& routes,
      const std::vector<Link<T>>& links) {
    const size_t R = routes.size();
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> J(R, R);

    std::vector<std::unordered_set<int>> route_sets(R);
    for (size_t i = 0; i < R; ++i) {
      route_sets[i].insert(routes[i].begin(), routes[i].end());
    }

    for (size_t i = 0; i < R; ++i) {
      T sum = 0;
      for (int link_id : routes[i]) {
        sum += links[link_id].DelayDer();
      }
      J(i, i) = sum;

      for (size_t j = i + 1; j < R; ++j) {
        T val = 0;
        for (int link_id : routes[j]) {
          if (route_sets[i].count(link_id)) {
            val += links[link_id].DelayDer();
          }
        }
        J(i, j) = val;
        J(j, i) = val;
      }
    }
    return J;
  }

}  // namespace TrafficAssignment

#endif  // LINK_COST_AGGREGATION_H

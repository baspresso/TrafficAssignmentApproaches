#ifndef NETWORK_METRICS_H
#define NETWORK_METRICS_H

#include <vector>
#include <limits>
#include "Network.h"

namespace TrafficAssignment {

  template <typename T>
  T ObjectiveFunction(const Network<T>& network) {
    T total = 0;
    for (const auto& link : network.links()) {
      total += link.DelayInteg();
    }
    return total;
  }

  template <typename T>
  T TotalTravelTime(const Network<T>& network) {
    T total = 0;
    for (const auto& link : network.links()) {
      total += link.flow * link.Delay();
    }
    return total;
  }

  template <typename T>
  T CalculatePathDelay(const Network<T>& network, const std::vector<int>& path) {
    T delay = 0;
    for (int link_id : path) {
      delay += network.links()[link_id].Delay();
    }
    return delay;
  }

  template <typename T>
  T RelativeGap(const Network<T>& network) {
    T numerator = 0;
    T denominator = 0;

    for (int origin = 0; origin < network.number_of_zones(); ++origin) {
      auto best_routes = network.ComputeSingleOriginBestRoutes(origin);
      for (const auto& [od_index, path] : best_routes) {
        numerator += network.od_pairs()[od_index].GetDemand() *
                     CalculatePathDelay(network, path);
      }
    }

    for (const auto& link : network.links()) {
      denominator += link.flow * link.Delay();
    }

    if (denominator < std::numeric_limits<T>::epsilon()) {
      return std::numeric_limits<T>::quiet_NaN();
    }

    return 1 - (numerator / denominator);
  }

}  // namespace TrafficAssignment

#endif  // NETWORK_METRICS_H

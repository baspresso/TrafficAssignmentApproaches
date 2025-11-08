#ifndef ROUTE_BASED_NEWTON_STEP_SHIFT_METHOD_H
#define ROUTE_BASED_NEWTON_STEP_SHIFT_METHOD_H

#include "./RouteBasedShiftMethod.h"
#include <set>
#include <chrono>
#include <thread>

namespace TrafficAssignment {
  template <typename T>
  class RouteBasedNewtonStepShiftMethod : public RouteBasedShiftMethod <T> {    
  public:
    explicit RouteBasedNewtonStepShiftMethod(Network<T>& network)
        : RouteBasedShiftMethod<T>(network) {}

    std::vector <T> FlowShift(int od_index) {
      auto& od_pair = this->od_pairs()[od_index];
      const int routes_count = od_pair.GetRoutesCount();
      auto [min_max_diff, min_index, max_index] = MinMaxRoutes(od_pair); 
      if (min_max_diff < computational_treshold_) {
        return od_pair.GetRoutesFlow();
      }
      auto [min_route_non_intersecting, max_route_non_intersecting] =
        MinMaxRoutesNonIntersectingLinks(od_pair, min_index, max_index);
      T delta = FlowAmount(min_route_non_intersecting, max_route_non_intersecting);
      std::vector <T> new_flows = od_pair.GetRoutesFlow();
      new_flows[min_index] += delta;
      new_flows[max_index] -= delta;
      return new_flows;
    }
  private:
    const T computational_treshold_ = 1e-10;
    std::tuple <T, int, int> MinMaxRoutes(const OriginDestinationPair <T>& od_pair) {
        std::vector <std::vector <int>> routes = od_pair.GetRoutes();
        const int routes_count = od_pair.GetRoutesCount();
        int min_index = 0, max_index = 0;
        T min_delay = Link<T>::GetLinksDelay(this->network_.links(), routes[0]);
        T max_delay = min_delay;
        for (int i = 1; i < routes_count; i++) {
            T delay = Link<T>::GetLinksDelay(this->network_.links(), routes[i]);
            if (delay < min_delay) {
              min_index = i;
              min_delay = delay;
            }
            if (delay > max_delay) {
              max_index = i;
              max_delay = delay;
            }
        }
        return {max_delay - min_delay, min_index, max_index};
    }

    std::tuple<std::vector<int>, std::vector<int>> 
    MinMaxRoutesNonIntersectingLinks(const OriginDestinationPair<T>& od_pair,
                                     int min_index, 
                                     int max_index) {
      std::vector <std::vector <int>> routes = od_pair.GetRoutes();
      std::vector <int> min_route = routes[min_index];
      std::vector <int> max_route = routes[max_index];
      std::set <int> intersecting_links; 
      for (auto min_route_link : min_route) {
        if (std::find(max_route.begin(), max_route.end(), min_route_link) != max_route.end()) {
          intersecting_links.insert(min_route_link);
        }
      }
      std::vector <int> min_route_non_intersecting;
      std::vector <int> max_route_non_intersecting;
      for (auto min_route_link : min_route) {
        if (!intersecting_links.count(min_route_link)) {
          min_route_non_intersecting.push_back(min_route_link);
        }
      }
      for (auto max_route_link : max_route) {
        if (!intersecting_links.count(max_route_link)) {
          max_route_non_intersecting.push_back(max_route_link);
        }
      }
      return {min_route_non_intersecting, max_route_non_intersecting};
    }

    T FlowAmount(std::vector <int> min_route_non_intersecting, std::vector <int> max_route_non_intersecting) {
      std::pair <T, T> routes_delays = { 0, 0 };
      T sum_routes_delays_der = 0;
      for (auto link_index : min_route_non_intersecting) {
        routes_delays.first += this->network_.links()[link_index].Delay(this->network_.links()[link_index].flow);
        sum_routes_delays_der += this->network_.links()[link_index].DelayDer(this->network_.links()[link_index].flow);
      }
      for (auto link_index : max_route_non_intersecting) {
        routes_delays.second += this->network_.links()[link_index].Delay(this->network_.links()[link_index].flow);
        sum_routes_delays_der += this->network_.links()[link_index].DelayDer(this->network_.links()[link_index].flow);
      }
      T delta = std::abs(routes_delays.first - routes_delays.second) / sum_routes_delays_der;
      return delta;
    }
  };
}

#endif // ROUTE_BASED_NEWTON_STEP_SHIFT_METHOD_H
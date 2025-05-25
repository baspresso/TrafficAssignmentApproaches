#ifndef ROUTE_BASED_KRYLATOV_2023_APPROACH_H
#define ROUTE_BASED_KRYLATOV_2023_APPROACH_H

#include "RouteBasedApproach.h"
#include "./components/RouteBasedKrylatov2023ShiftMethod.h"

namespace TrafficAssignment {
  template <typename T>
  class RouteBasedKrylatov2023Approach : public RouteBasedApproach <T> {
  public:
    RouteBasedKrylatov2023Approach(Network<T>& network, T alpha = 1e-14)
      : RouteBasedApproach<T>(network, alpha) {
        this->shift_method_ = new RouteBasedKrylatov2023ShiftMethod<T>(
            network
        );
    }

    ~RouteBasedKrylatov2023Approach() = default;

    std::string GetApproachName() override {
      return "RouteBasedKrylatov2023";
    }
  };
}

#endif 
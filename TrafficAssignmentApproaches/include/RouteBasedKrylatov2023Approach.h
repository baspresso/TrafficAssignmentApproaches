#ifndef ROUTE_BASED_KRYLATOV_2023_APPROACH_H
#define ROUTE_BASED_KRYLATOV_2023_APPROACH_H

#include "RouteBasedApproach.h"
#include "RouteBasedKrylatov2023ShiftMethod.h"

namespace TrafficAssignment {
  template <typename T>
  class RouteBasedKrylatov2023Approach : public RouteBasedApproach <T> {
  public:
    RouteBasedKrylatov2023Approach(std::string dataset_name, T alpha = 1e-6) :
      RouteBasedApproach<T>(dataset_name, alpha) {
      this->shift_method_ = new RouteBasedKrylatov2023ShiftMethod<T>(this->GetLinksRef(), this->GetOriginDestinationPairsRef());
    }

    ~RouteBasedKrylatov2023Approach() = default;
  };
}

#endif ROUTE_BASED_KRYLATOV_2023_APPROACH_H
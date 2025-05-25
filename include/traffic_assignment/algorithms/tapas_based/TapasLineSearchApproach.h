#ifndef TAPAS_LINE_SEARCH_APPROACH_H
#define TAPAS_LINE_SEARCH_APPROACH_H

#include "TapasApproach.h"
#include "./components/TapasLineSearchShiftMethod.h"

namespace TrafficAssignment {
  template <typename T>
  class TapasLineSearchApproach : public TapasApproach <T> {
  public:
    TapasLineSearchApproach(Network<T>& network, T alpha = 1e-6) :
      TapasApproach<T>(network, alpha) {
      this->shift_method_ = new TapasLineSearchShiftMethod<T>(this->GetLinksRef());
    }

    ~TapasLineSearchApproach() = default;

    std::string GetApproachName() override {
      return "TapasLineSearch";
    }
  };
}

#endif
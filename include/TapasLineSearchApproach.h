#ifndef TAPAS_LINE_SEARCH_APPROACH_H
#define TAPAS_LINE_SEARCH_APPROACH_H

#include "TapasApproach.h"
#include "TapasLineSearchShiftMethod.h"

namespace TrafficAssignment {
  template <typename T>
  class TapasLineSearchApproach : public TapasApproach <T> {
  public:
    TapasLineSearchApproach(std::string dataset_name, T alpha = 1e-14) :
      TapasApproach<T>(dataset_name, alpha) {
      this->shift_method_ = new TapasLineSearchShiftMethod<T>(this->GetLinksRef());
    }

    ~TapasLineSearchApproach() = default;

    std::string GetApproachName() override {
      return "TapasLineSearch";
    }
  };
}

#endif
#ifndef TAPAS_NEWTON_STEP_APPROACH_H
#define TAPAS_NEWTON_STEP_APPROACH_H

#include "TapasApproach.h"
#include "./components/TapasNewtonStepShiftMethod.h"

namespace TrafficAssignment {
  template <typename T>
  class TapasNewtonStepApproach : public TapasApproach <T> {
  public:
    TapasNewtonStepApproach(Network<T>& network, T alpha = 1e-6) :
      TapasApproach<T>(network, alpha) {
      this->shift_method_ = new TapasNewtonStepShiftMethod<T>(this->GetLinksRef());
    }

    ~TapasNewtonStepApproach() = default;

    std::string GetApproachName() override {
      return "TapasNewtonStep";
    }
  };
}

#endif
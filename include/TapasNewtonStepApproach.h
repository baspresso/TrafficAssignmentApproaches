#ifndef TAPAS_NEWTON_STEP_APPROACH_H
#define TAPAS_NEWTON_STEP_APPROACH_H

#include "TapasApproach.h"
#include "TapasNewtonStepShiftMethod.h"

namespace TrafficAssignment {
  template <typename T>
  class TapasNewtonStepApproach : public TapasApproach <T> {
  public:
    TapasNewtonStepApproach(std::string dataset_name, T alpha = 1e-14) :
      TapasApproach<T>(dataset_name, alpha) {
      this->shift_method_ = new TapasNewtonStepShiftMethod<T>(this->GetLinksRef());
    }

    ~TapasNewtonStepApproach() = default;

    std::string GetApproachName() override {
      return "TapasNewtonStep";
    }
  };
}

#endif
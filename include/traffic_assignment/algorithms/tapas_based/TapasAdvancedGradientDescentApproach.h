#ifndef TAPAS_ADVANCED_GRADIENT_DESCENT_APPROACH_H
#define TAPAS_ADVANCED_GRADIENT_DESCENT_APPROACH_H

#include "TapasApproach.h"
#include "TapasAdvancedGradientDescentShiftMethod.h"

namespace TrafficAssignment {
  template <typename T>
  class TapasAdvancedGradientDescentAppoach : public TapasApproach <T> {
  public:
    TapasAdvancedGradientDescentAppoach(Network<T>& network, T alpha = 1e-6) :
      TapasApproach<T>(network, alpha) {
      this->shift_method_ = new TapasAdvancedGradientDescentShiftMethod<T>(this->GetLinksRef());
    }

    ~TapasAdvancedGradientDescentAppoach() = default;

    std::string GetApproachName() override {
      return "TapasAdvancedGradientDescent";
    }
  };
}

#endif
#ifndef TAPAS_NEWTON_STEP_SHIFT_METHOD_H
#define TAPAS_NEWTON_STEP_SHIFT_METHOD_H

#include "TapasShiftMethod.h"

namespace TrafficAssignment {
  /**
   * @brief Newton step shift method for PAS-based flow redistribution.
   *
   * Computes the shift amount using the first-order Newton approximation:
   *   delta = |C(s1) - C(s2)| / sum(dc_a/df_a)
   * where the sum runs over all links in both PAS segments.
   * Adapts Perederieieva et al. (2015) Eq. 14 from route pairs to PAS segments.
   *
   * @tparam T Numeric type for flow computations.
   */
  template <typename T>
  class TapasNewtonStepShiftMethod : public TapasShiftMethod <T> {
  public:
    TapasNewtonStepShiftMethod(const std::vector <Link <T>>& links, T computation_threshold = T(1e-10))
      : TapasShiftMethod<T>(links, computation_threshold) {}

    ~TapasNewtonStepShiftMethod() {}

    /// @brief Computes Newton-step PAS flow shift, clamped to available flow.
    std::pair <T, T> FlowShift(T starting_point, const std::pair <std::vector <int>, std::vector <int>>& pas, std::pair <T, T> total_flow) override {
      // Flow is being tranfered from / to the first part of the PAS if direction is true / false
      std::pair <T, T> pas_delay = { 0, 0 };
      for (auto link_index : pas.first) {
        pas_delay.first += this->links_[link_index].Delay(this->links_[link_index].flow);
      }
      for (auto link_index : pas.second) {
        pas_delay.second += this->links_[link_index].Delay(this->links_[link_index].flow);
      }
      bool direction = (pas_delay.first > pas_delay.second); 
      T sum_pas_delays_der = 0;
      for (auto link_index : pas.first) {
        sum_pas_delays_der += this->links_[link_index].DelayDer(this->links_[link_index].flow);
      }
      for (auto link_index : pas.second) {
        sum_pas_delays_der += this->links_[link_index].DelayDer(this->links_[link_index].flow);
      }
      if (sum_pas_delays_der == 0) {
        return {0, 0};
      }
      T delta = std::abs(pas_delay.first - pas_delay.second) / sum_pas_delays_der;
      if (direction) {
        delta = std::min(total_flow.first, delta);
        return {-delta, delta};
      }
      else {
        delta = std::min(total_flow.second, delta);
        return {delta, -delta};
      }
    }
  };
}

#endif
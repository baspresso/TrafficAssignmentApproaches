#ifndef TAPAS_NEWTON_STEP_SHIFT_METHOD_H
#define TAPAS_NEWTON_STEP_SHIFT_METHOD_H

#include "TapasShiftMethod.h"

namespace TrafficAssignment {
  template <typename T>
  class TapasNewtonStepShiftMethod : public TapasShiftMethod <T> {
  public:
    TapasNewtonStepShiftMethod(std::vector <Link <T>>& links) : TapasShiftMethod<T>(links) {}

    ~TapasNewtonStepShiftMethod() {}

    std::pair <T, T> FlowShift(T starting_point, std::pair <std::vector <int>, std::vector <int>> pas, std::pair <T, T> total_flow) override {
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
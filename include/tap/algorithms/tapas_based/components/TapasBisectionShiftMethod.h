#ifndef TAPAS_BISECTION_SHIFT_METHOD_H
#define TAPAS_BISECTION_SHIFT_METHOD_H

#include "TapasShiftMethod.h"

namespace TrafficAssignment {
  /**
   * @brief Bisection line search shift method for PAS flow redistribution.
   *
   * Binary search on [0, max_shift] to find the equilibrating flow shift.
   * Halves interval until width <= 2 * precision. Based on TAsK reference
   * implementation's Bisection class.
   *
   * @tparam T Numeric type for flow computations.
   */
  template <typename T>
  class TapasBisectionShiftMethod : public TapasShiftMethod<T> {
  public:
    TapasBisectionShiftMethod(const std::vector<Link<T>>& links, T computation_threshold = T(1e-10))
      : TapasShiftMethod<T>(links, computation_threshold) {}

    ~TapasBisectionShiftMethod() {}

    std::pair<T, T> FlowShift(T starting_point, const std::pair<std::vector<int>, std::vector<int>>& pas, std::pair<T, T> total_flow) override {
      bool direction = this->FlowShiftResult(pas, {0, 0});
      T max_shift = direction ? total_flow.first : total_flow.second;

      if (max_shift < this->computation_threshold_) {
        return {0, 0};
      }

      T a = T(0);
      T b = max_shift;

      // Verify that shifting max_shift actually reverses the direction
      // If not, just shift the full amount
      if (this->DirectionFlowShiftResult(pas, b, direction) == direction) {
        if (direction) {
          return {-b, b};
        } else {
          return {b, -b};
        }
      }

      // Bisection: find the zero-crossing of the cost difference
      while ((b - a) > T(2) * this->computation_threshold_) {
        T mid = (a + b) / T(2);
        if (this->DirectionFlowShiftResult(pas, mid, direction) == direction) {
          a = mid;
        } else {
          b = mid;
        }
      }

      T flow_shift = (a + b) / T(2);
      if (flow_shift < this->computation_threshold_) {
        return {0, 0};
      }

      if (direction) {
        return {-flow_shift, flow_shift};
      } else {
        return {flow_shift, -flow_shift};
      }
    }
  };
}

#endif

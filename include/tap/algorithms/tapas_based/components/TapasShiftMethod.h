#ifndef TAPAS_SHIFT_METHOD_H
#define TAPAS_SHIFT_METHOD_H

#include <vector>
#include "../../../data/Link.h"

namespace TrafficAssignment {
  /**
   * @brief Abstract interface for PAS (Paired Alternative Segment) flow shift methods.
   *
   * Computes how much flow to transfer between the two segments of a PAS to reduce
   * disequilibrium. The shift procedure is described in Bar-Gera (2010) Table 3.
   *
   * @tparam T Numeric type for flow computations.
   */
  template <typename T>
  class TapasShiftMethod {
  public:
    TapasShiftMethod(const std::vector <Link <T>>& links, T computation_threshold = T(1e-10))
      : links_(links), computation_threshold_(computation_threshold) {}

    virtual ~TapasShiftMethod() = default;

    /**
     * @brief Computes flow shift amounts for both PAS segments.
     * @param starting_point Initial step size hint (from previous shift on this PAS).
     * @param pas The paired alternative segment: (segment_1 links, segment_2 links).
     * @param total_flow Total transferable flow on (segment_1, segment_2).
     * @return (delta_1, delta_2) where delta_1 + delta_2 = 0 (flow conservation).
     */
    virtual std::pair <T, T> FlowShift(T starting_point, const std::pair <std::vector <int>, std::vector <int>>& pas, std::pair <T, T> total_flow) = 0;

  protected:
    const std::vector <Link <T>>& links_;

    const T computation_threshold_;

    /// @brief Checks if segment 1 would still be costlier than segment 2 after applying flow_shift.
    /// @return True if cost(segment_1 + shift_1) > cost(segment_2 + shift_2).
    bool FlowShiftResult(const std::pair <std::vector <int>, std::vector <int>>& pas, const std::pair <T, T>& flow_shift) {
      std::pair <T, T> pas_delay = { 0, 0 };
      for (auto link_index : pas.first) {
        pas_delay.first += this->links_[link_index].Delay(this->links_[link_index].flow + flow_shift.first);
      }
      for (auto link_index : pas.second) {
        pas_delay.second += this->links_[link_index].Delay(this->links_[link_index].flow + flow_shift.second);
      }
      return (pas_delay.first > pas_delay.second);
    }

    /// @brief Directional variant of FlowShiftResult: tests shift in the given direction.
    bool DirectionFlowShiftResult(const std::pair <std::vector <int>, std::vector <int>>& pas, const T& flow_shift, const bool direction) {
      if (direction) {
        return FlowShiftResult(pas, { -flow_shift, flow_shift });
      }
      else {
        return FlowShiftResult(pas, { flow_shift, -flow_shift });
      }
    }
  };
}

#endif
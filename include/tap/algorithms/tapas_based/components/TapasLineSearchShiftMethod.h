#ifndef TAPAS_LINE_SEARCH_SHIFT_METHOD_H
#define TAPAS_LINE_SEARCH_SHIFT_METHOD_H

#include "TapasShiftMethod.h"

namespace TrafficAssignment {
  template <typename T>
  class TapasLineSearchShiftMethod : public TapasShiftMethod <T> {
  public:
    TapasLineSearchShiftMethod(const std::vector <Link <T>>& links) : TapasShiftMethod<T>(links) {}

    ~TapasLineSearchShiftMethod() {}

    std::pair <T, T> FlowShift(T starting_point, std::pair <std::vector <int>, std::vector <int>> pas, std::pair <T, T> total_flow) override {
      T flow_shift = starting_point;
      // Flow is being tranfered from / to the first part of the PAS if direction is true / false
      bool direction;
      if (this->FlowShiftResult(pas, { 0, 0 })) {
        direction = true;
      }
      else {
        direction = false;
      }
      T total_flow_transfering_from;
      if (direction) {
        total_flow_transfering_from = total_flow.first;
      }
      else {
        total_flow_transfering_from = total_flow.second;
      }
      flow_shift = std::min(flow_shift, total_flow_transfering_from);
      while (this->DirectionFlowShiftResult(pas, std::min(flow_shift * 2, total_flow_transfering_from), direction) == direction) {
        flow_shift *= 2;
      }
      while (this->DirectionFlowShiftResult(pas, flow_shift, direction) != direction) {
        flow_shift /= 2;
      }
      if (direction) {
        return { -flow_shift, flow_shift };
      }
      else {
        return { flow_shift, -flow_shift };
      }
    }
  };
}

#endif
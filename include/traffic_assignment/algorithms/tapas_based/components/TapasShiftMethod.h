#ifndef TAPAS_SHIFT_METHOD_H
#define TAPAS_SHIFT_METHOD_H

#include <vector>
#include "../../../data/Link.h"

namespace TrafficAssignment {
  template <typename T>
  class TapasShiftMethod {
  public:
    TapasShiftMethod(std::vector <Link <T>>& links) : links_(links) {}
    
    virtual ~TapasShiftMethod() = default;

    virtual std::pair <T, T> FlowShift(T starting_point, std::pair <std::vector <int>, std::vector <int>> pas, std::pair <T, T> total_flow) = 0;
  
  protected:
    std::vector <Link <T>>& links_;

    const double computational_treshold_ = 1e-7;

    bool FlowShiftResult(const std::pair <std::vector <int>, std::vector <int>> pas, const std::pair <T, T> flow_shift) {
      std::pair <T, T> pas_delay = { 0, 0 };
      for (auto link_index : pas.first) {
        pas_delay.first += this->links_[link_index].Delay(this->links_[link_index].flow + flow_shift.first);
      }
      for (auto link_index : pas.second) {
        pas_delay.second += this->links_[link_index].Delay(this->links_[link_index].flow + flow_shift.second);
      }
      return (pas_delay.first > pas_delay.second);
    }

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
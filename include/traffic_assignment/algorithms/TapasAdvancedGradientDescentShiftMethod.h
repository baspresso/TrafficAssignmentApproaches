#ifndef TAPAS_ADVANCED_GRADIENT_DESCENT_SHIFT_METHOD_H
#define TAPAS_ADVANCED_GRADIENT_DESCENT_SHIFT_METHOD_H

#include "TapasShiftMethod.h"

namespace TrafficAssignment {
  template <typename T>
  class TapasAdvancedGradientDescentShiftMethod : public TapasShiftMethod <T> {
  public:
    TapasAdvancedGradientDescentShiftMethod(std::vector <Link <T>>& links) : TapasShiftMethod<T>(links) {}

    ~TapasAdvancedGradientDescentShiftMethod() {}

    std::pair <T, T> FlowShift(T starting_point, std::pair <std::vector <int>, std::vector <int>> pas, std::pair <T, T> total_flow) override {

      T flow_shift = 0;
      // Flow is being tranfered from / to the first part of the PAS if direction is true / false
      bool direction;
      if (this->FlowShiftResult(pas, { 0, 0 })) {
        direction = true;
      }
      else {
        direction = false;
      }
      double min_flow = std::min(total_flow.first, total_flow.second);
      double max_flow = std::max(total_flow.first, total_flow.second);
      if (min_flow < this->computational_treshold_) {
        if (direction) {
          return { -max_flow / 2, max_flow / 2 };
        }
        else {
          return { max_flow / 2, -max_flow / 2 };
        }
      }
      std::pair <T, T> pas_der = { Link<T>::GetLinksDelayDer(this->links_, pas.first), Link<T>::GetLinksDelayDer(this->links_, pas.second) };
      std::pair <T, T> pas_second_der = { Link<T>::GetLinksDelaySecondDer(this->links_, pas.first), Link<T>::GetLinksDelaySecondDer(this->links_, pas.second) };
      std::pair <T, T> ders_ratio = { pas_der.first / (pas_second_der.first), pas_der.second / (pas_second_der.second) };
      flow_shift = -ders_ratio.first + (ders_ratio.first + ders_ratio.second) / (pas_second_der.first) / (1 / (pas_second_der.first) + 1 / (pas_second_der.second));
      //cout << flow_shift << '\n';
      //if (flow_shift < 0) {
      //  flow_shift = -std::min(total_flow.first, -flow_shift);
      //}
      //else {
      //  flow_shift = std::min(total_flow.second, flow_shift);
      //}

      //while (ObjectiveFunctionIncreaseResult(flow_shift, pas)) {
      //  flow_shift /= 2;
      //}

      std::cout << flow_shift << '\n';
      return { flow_shift, -flow_shift };
    }
  private:
    bool ObjectiveFunctionIncreaseResult(T flow_shift, const std::pair <std::vector <int>, std::vector <int>>& pas) {
      T current_res = 0, next_res = 0;
      for (auto link_index : pas.first) {
        current_res += this->links_[link_index].DelayInteg();
        next_res += this->links_[link_index].DelayInteg(this->links_[link_index].flow + flow_shift);
      }
      for (auto link_index : pas.second) {
        current_res += this->links_[link_index].DelayInteg();
        next_res += this->links_[link_index].DelayInteg(this->links_[link_index].flow - flow_shift);
      }
      return current_res < next_res;
    }
  };
}

#endif
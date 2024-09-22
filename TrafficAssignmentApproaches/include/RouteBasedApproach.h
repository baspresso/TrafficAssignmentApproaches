#ifndef ROUTE_BASED_APPROACH_H
#define ROUTE_BASED_APPROACH_H

#include <queue>
#include <string>
#include "TrafficAssignmentApproach.h"
#include "RouteBasedShiftMethod.h"

namespace TrafficAssignment {
  template <typename T>
  class RouteBasedApproach : public TrafficAssignmentApproach <T> {
  public:
    RouteBasedApproach(std::string dataset_name, T alpha = 1e-6) :
      TrafficAssignmentApproach<T>::TrafficAssignmentApproach(dataset_name, alpha) {
      shift_method_ = nullptr;
    }
    void ComputeTrafficFlows() override { 
      int iteration_count = 0;
      while (iteration_count++ < 10) {
        for (int origin_index = 0; origin_index < this->number_of_zones_; origin_index++) {
          this->SingleOriginBestRoutes(origin_index);
          for (int cnt = 0; cnt < origin_cycle_count; cnt++) {
            for (auto now : this->origin_info_[origin_index]) {
              PerformShift(now.second);
            }
          }
        }
        
        for (int count = 0; count < full_cycle_count; count++) {
          for (int od_pair_index = 0; od_pair_index < this->number_of_origin_destination_pairs_; od_pair_index++) {
            PerformShift(od_pair_index);
            //this->GetStatistics();
          }
          //this->GetStatistics();
        }
      }
    }
  protected:
    const int full_cycle_count = 20;
    const int origin_cycle_count = 5;
    RouteBasedShiftMethod <T>* shift_method_;

    void PerformShift(int od_pair_index) {
      shift_method_->PerformShift(od_pair_index);
    }
  };
}

#endif ROUTE_BASED_APPROACH_H
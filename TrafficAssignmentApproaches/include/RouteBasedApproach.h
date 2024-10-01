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

    ~RouteBasedApproach() {
      delete shift_method_;
    }

    void ComputeTrafficFlows() override { 
      this->statistics_recorder_.StartRecording(this, this->dataset_name_);
      int iteration_count = 0;
      for (int origin_index = 0; origin_index < this->number_of_zones_; origin_index++) {
        auto best_routes = this->SingleOriginBestRoutes(origin_index);
        this->AddNewOriginDestinationRoutes(best_routes);
      }
      while (iteration_count++ < 10) {
        for (int origin_index = 0; origin_index < this->number_of_zones_; origin_index++) {
          auto best_routes = this->SingleOriginBestRoutes(origin_index);
          this->AddNewOriginDestinationRoutes(best_routes);
          this->statistics_recorder_.RecordStatistics();

          for (int cnt = 0; cnt < origin_iteration_count; cnt++) {
            for (auto now : this->origin_info_[origin_index]) {
              if (this->origin_destination_pairs_[now.second].GetRoutesCount() > 1) {
                ExecuteFlowShift(now.second);
              }
            }
          }
        }
        
        for (int count = 0; count < full_iteration_count; count++) {
          for (int od_pair_index = 0; od_pair_index < this->number_of_origin_destination_pairs_; od_pair_index++) {
            ExecuteFlowShift(od_pair_index);
          }
          this->statistics_recorder_.RecordStatistics();
        }
      }
    }

  protected:
    const int full_iteration_count = 20;
    const int origin_iteration_count = 5;
    RouteBasedShiftMethod <T>* shift_method_;

    std::vector <T> FlowShift(int od_pair_index) {
      return shift_method_->FlowShift(od_pair_index);
    }

    std::vector <Link<T>>& GetLinksRef() {
      return this->links_;
    }

    std::vector <OriginDestinationPair<T>>& GetOriginDestinationPairsRef() {
      return this->origin_destination_pairs_;
    }

    void ExecuteFlowShift(int od_pair_index) {
      std::vector <T> flow_shift = FlowShift(od_pair_index);
      this->origin_destination_pairs_[od_pair_index].SetRoutesFlow(flow_shift);
    }

  };
}

#endif ROUTE_BASED_APPROACH_H
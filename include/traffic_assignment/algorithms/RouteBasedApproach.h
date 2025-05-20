#ifndef ROUTE_BASED_APPROACH_H
#define ROUTE_BASED_APPROACH_H

#include <queue>
#include <string>
#include "TrafficAssignmentApproach.h"
#include "RouteBasedShiftMethod.h"
#include <random>

namespace TrafficAssignment {
  template <typename T>
  class RouteBasedApproach : public TrafficAssignmentApproach <T> {
  public:
    RouteBasedApproach(std::string dataset_name, T alpha = 1e-14) :
      TrafficAssignmentApproach<T>::TrafficAssignmentApproach(dataset_name, alpha) {
      shift_method_ = nullptr;
      objective_function_expected_decrease_.resize(this->number_of_origin_destination_pairs_, 0);
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
      while (iteration_count++ < 100) {
        
        for (int od_index = 0; od_index < this->number_of_origin_destination_pairs_; od_index++) {
          objective_function_expected_decrease_[od_index] = 0;
        }
        for (int origin_index = 0; origin_index < this->number_of_zones_; origin_index++) {
          auto best_routes = this->SingleOriginBestRoutes(origin_index);
          this->AddNewOriginDestinationRoutes(best_routes);
          //this->statistics_recorder_.RecordStatistics();

          for (int cnt = 0; cnt < origin_iteration_count; cnt++) {
            for (auto now : this->origin_info_[origin_index]) {
              if (this->origin_destination_pairs_[now.second].GetRoutesCount() > 1) {
                ExecuteFlowShift(now.second);
              }
            }
          }
        }
        
        std::priority_queue <std::pair <T, int>> od_queue = OriginDestinationQueuePreparation();
        
        if (od_queue.empty()) {
          continue;
        }

        for (int count = 0; count < full_iteration_count; count++) {
          for (int od_pair_index = 0; od_pair_index < this->number_of_origin_destination_pairs_; od_pair_index++) {
            OriginDestinationPairProcessing(od_queue);
          }
          //this->statistics_recorder_.RecordStatistics();
        }
        this->statistics_recorder_.RecordStatistics();
      }
    }

  protected:
    const int full_iteration_count = 3;
    const int origin_iteration_count = 1;
    const T alpha = 0.7;
    RouteBasedShiftMethod <T>* shift_method_;
    std::vector <T> objective_function_expected_decrease_;

    std::vector <T> FlowShift(int od_pair_index) {
      return shift_method_->FlowShift(od_pair_index);
    }

    std::vector <Link<T>>& GetLinksRef() {
      return this->links_;
    }


    std::vector <OriginDestinationPair<T>>& GetOriginDestinationPairsRef() {
      return this->origin_destination_pairs_;
    }

    std::vector <int> GetRandomOrder(int number_of_elements) {
      std::vector <int> order(number_of_elements);
      for (int i = 0; i < number_of_elements; i++) {
        order[i] = i;
      }
      std::random_device rd;
      std::mt19937 g(rd());
      std::shuffle(order.begin(), order.end(), g);
      return order;
    }

    std::priority_queue <std::pair <T, int>> OriginDestinationQueuePreparation() {
      std::priority_queue <std::pair <T, int>> od_queue;
      for (int od_pair_index = 0; od_pair_index < this->number_of_origin_destination_pairs_; od_pair_index++) {
        if (this->origin_destination_pairs_[od_pair_index].GetRoutesCount() > 1) {
          od_queue.push({ objective_function_expected_decrease_[od_pair_index], od_pair_index });
        }
      }
      return od_queue;
    }

    void UpdateObjectiveFunctionExpectedDecrease(const T& objective_func_delta, int od_pair_index) {
      if (objective_function_expected_decrease_[od_pair_index] == 0) {
        objective_function_expected_decrease_[od_pair_index] = objective_func_delta;
      }
      else {
        objective_function_expected_decrease_[od_pair_index] = (1 - alpha) * objective_func_delta + alpha * objective_function_expected_decrease_[od_pair_index];
      }
    }

    void ExecuteFlowShift(int od_pair_index) {
      std::vector <std::vector <int>> od_routes = this->origin_destination_pairs_[od_pair_index].GetRoutes();
      T objective_func_before_shift = this->RoutesLinksObjectiveFunction(od_routes);
      std::vector <T> flow_shift = FlowShift(od_pair_index);
      this->origin_destination_pairs_[od_pair_index].SetRoutesFlow(flow_shift);
      T objective_func_after_shift = this->RoutesLinksObjectiveFunction(od_routes);
      UpdateObjectiveFunctionExpectedDecrease(std::abs(objective_func_before_shift - objective_func_after_shift), od_pair_index);
    }

    void OriginDestinationPairProcessing(std::priority_queue <std::pair <T, int>>& od_queue) {
      int od_pair_index = od_queue.top().second;
      od_queue.pop();
      ExecuteFlowShift(od_pair_index);
      od_queue.push({ objective_function_expected_decrease_[od_pair_index], od_pair_index });
    }

  };
}

#endif
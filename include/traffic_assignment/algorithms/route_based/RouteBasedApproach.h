#ifndef ROUTE_BASED_APPROACH_H
#define ROUTE_BASED_APPROACH_H

#include <queue>
#include <string>
#include "../common/TrafficAssignmentApproach.h"
#include "./components/RouteBasedShiftMethod.h"
#include "RouteBasedShiftMethodFactory.h"
#include <random>

namespace TrafficAssignment {
  template <typename T>
  class RouteBasedApproach : public TrafficAssignmentApproach <T> {
  public:
    RouteBasedApproach(Network<T>& network,
                       T alpha = 1e-14,
                       std::string shift_method_name = "Krylatov2023")
      : TrafficAssignmentApproach<T>::TrafficAssignmentApproach(network, alpha) 
    {
        shift_method_ = RouteBasedShiftMethodFactory<T>::GetInstance().Create(shift_method_name, this->network_);
        shift_method_name_ = shift_method_name;
        InitializeState();
    }

    ~RouteBasedApproach() = default;

    void ComputeTrafficFlows(bool statistics_recording = false) override { 
      if (statistics_recording) {
        this->statistics_recorder_.StartRecording(this->GetApproachName());
      }
      int iteration_count = 0;

      InitializeRoutes();
      
      if (statistics_recording) {
        this->statistics_recorder_.RecordStatistics();
      }

      while (iteration_count++ < 200) {
        ResetExpectedDecreases();
        ProcessOrigins();
        ProcessODPairs(statistics_recording);
        if (statistics_recording) {
          this->statistics_recorder_.RecordStatistics();
        }
      }
    }

  protected:
    const int full_iteration_count_ = 3;
    const int origin_iteration_count_ = 1;
    const T alpha_ = 0.7;
    std::unique_ptr <RouteBasedShiftMethod <T>> shift_method_;
    std::string shift_method_name_;
    std::vector <T> objective_function_expected_decrease_;

    void InitializeState() {
        objective_function_expected_decrease_.resize(
            this->network().number_of_od_pairs(), 0
        );
    }

    void InitializeRoutes() {
        for (int origin_index = 0; origin_index < this->network().number_of_zones(); origin_index++) {
            auto best_routes = this->network().ComputeSingleOriginBestRoutes(origin_index);
            this->network().AddRoutes(best_routes);
        }
    }

    void ResetExpectedDecreases() {
        std::fill(
            objective_function_expected_decrease_.begin(),
            objective_function_expected_decrease_.end(),
            T(0)
        );
    }

    void ProcessOrigins() {
        for (int origin_index = 0; origin_index < this->network().number_of_zones(); origin_index++) {
            UpdateOriginRoutes(origin_index);
            ProcessOriginFlows(origin_index);
        }
    }

    void UpdateOriginRoutes(int origin_index) {
        auto best_routes = this->network().ComputeSingleOriginBestRoutes(origin_index);
        this->network().AddRoutes(best_routes);
    }

    void ProcessOriginFlows(int origin_index) {
        for (int cnt = 0; cnt < origin_iteration_count_; cnt++) {
            ProcessOriginODPairs(origin_index);
        }
    }

    void ProcessOriginODPairs(int origin_index) {
        for (const auto& [dest, od_index] : this->network().origin_info()[origin_index]) {
            if (this->network().od_pairs()[od_index].GetRoutesCount() > 1) {
                ExecuteFlowShift(od_index);
            }
        }
    }

    void ProcessODPairs(bool statistics_recording) {
        auto od_queue = PrepareODQueue();
        if (od_queue.empty()) return;

        for (int count = 0; count < full_iteration_count_; count++) {
            ProcessODQueue(od_queue);
            if (count % 10 == 1) {
                if (statistics_recording) {
                  this->statistics_recorder_.RecordStatistics();
                }
            }
        }
    }

    std::priority_queue<std::pair<T, int>> PrepareODQueue() {
        std::priority_queue<std::pair<T, int>> od_queue;
        const auto& od_pairs = this->network().od_pairs();
        
        for (int od_index = 0; od_index < od_pairs.size(); od_index++) {
            if (od_pairs[od_index].GetRoutesCount() > 1) {
                od_queue.emplace(objective_function_expected_decrease_[od_index], od_index);
            }
        }
        return od_queue;
    }

    void ProcessODQueue(std::priority_queue<std::pair<T, int>>& od_queue) {
        auto temp_queue = od_queue;  // Work on a copy of the queue
        while (!temp_queue.empty()) {
            const auto [_, od_index] = temp_queue.top();
            temp_queue.pop();
            ExecuteFlowShift(od_index);
        }
    }


    void ExecuteFlowShift(int od_index) {
        auto& od_pair = this->network().mutable_od_pairs()[od_index];
        const auto routes = od_pair.GetRoutes();
        const T before = CalculateRoutesObjective(routes);
        if (this->network().od_pairs()[od_index].GetRoutesCount() == 1) {
            UpdateExpectedDecrease(od_index, before, CalculateRoutesObjective(routes));
            return;
        }
        const auto flow_shift = shift_method_->FlowShift(od_index);
        od_pair.SetRoutesFlow(flow_shift);
        
        UpdateExpectedDecrease(od_index, before, CalculateRoutesObjective(routes));
    }

    void UpdateExpectedDecrease(int od_index, T before, T after) {
        T delta = std::abs(before - after);
        auto& current = objective_function_expected_decrease_[od_index];
        
        current = (current == 0) ? delta : 
            (1 - alpha_) * delta + alpha_ * current;
    }
    
    std::vector <T> FlowShift(int od_pair_index) {
      return shift_method_->FlowShift(od_pair_index);
    }

    T CalculateRoutesObjective(const std::vector<std::vector<int>>& routes) {
        T total = 0;
        std::unordered_set<int> unique_links;
        
        for (const auto& route : routes) {
            unique_links.insert(route.begin(), route.end());
        }
        
        for (int link_id : unique_links) {
            total += this->network().links()[link_id].DelayInteg();
        }
        
        return total;
    }

    std::string GetApproachName() override {
      return "RouteBased" + shift_method_name_;
    }
  };
}

#endif
#ifndef DEMAND_BASED_APPROACH_H
#define DEMAND_BASED_APPROACH_H

#include "TapasApproach.h"
#include <limits>
#include <stack>

namespace TrafficAssignment {

  template <typename T>
  class DemandBasedApproach : public TrafficAssignmentApproach <T> {
  public:
    DemandBasedApproach(std::string dataset_name, T alpha = 1e-14) :
      TrafficAssignmentApproach<T>::TrafficAssignmentApproach(dataset_name, alpha) {
        node_dest_demand_.assign(this->number_of_nodes_, std::vector <T> (this->number_of_zones_, 0));
        node_dest_delay_estimation_.assign(this->number_of_nodes_, std::vector <T> (this->number_of_zones_, -1));
        link_dest_flow_.assign(this->number_of_links_, std::vector <T> (this->number_of_zones_, 0));
        zone_distances_.assign(this->number_of_zones_, std::vector <T> (this->number_of_nodes_, INF));
        zone_processed_.assign(this->number_of_zones_, std::vector <bool> (this->number_of_nodes_, false));
        zone_prev_node_.assign(this->number_of_zones_, std::vector <int> (this->number_of_nodes_, 0));
        zone_used_link_.assign(this->number_of_zones_, std::vector <int> (this->number_of_nodes_, 0));
        for (auto odp : this->origin_destination_pairs_) {
          auto origin_dest = odp.GetOriginDestination();
          node_dest_demand_[origin_dest.second][origin_dest.first] = odp.GetDemand();
        }
      }

    ~DemandBasedApproach() {
    }

    void ComputeTrafficFlows() override {
      this->statistics_recorder_.StartRecording(this, this->dataset_name_);
      int iteration_count = 0;
      while (iteration_count++ < 10000) {
        //std::cout << iteration_count << '\n';
        for (int zone_index = 0; zone_index < this->number_of_zones_; zone_index++) {
          ZoneProcessing(zone_index, k_);
        }
        RecordStatistics();
        FlowPumpOut(1.0 / ((iteration_count + 2) / std::log2(iteration_count + 2)));
      }
      FinishTransfer();
      //FullFlowTransfer();
      std::cout <<  std::setprecision(20) << this->ObjectiveFunction() << ' ' << this->RelativeGap() << '\n';
    }

    std::string GetApproachName() override {
      return "DemandBased";
    }

  protected:
    std::vector <std::vector <T>> node_dest_demand_;
    std::vector <std::vector <T>> node_dest_delay_estimation_;
    std::vector <std::vector <T>> link_dest_flow_;
    std::vector <std::vector <T>> zone_distances_;
    std::vector <std::vector <bool>> zone_processed_;
    std::vector <std::vector <int>> zone_prev_node_;
    std::vector <std::vector <int>> zone_used_link_;
    std::vector <std::vector <T>> original_zone_demand_;
    //const T INF = std::numeric_limits<T>::max();
    const T INF = 1e9;
    const T k_ = 0.5;

    void ZoneProcessing(int zone_index, T beta) {
      std::priority_queue <std::pair <T, int>, std::vector <std::pair <T, int>>, std::greater <std::pair <T, int>>> pq;
      std::fill(zone_distances_[zone_index].begin(), zone_distances_[zone_index].end(), INF);
      std::fill(zone_processed_[zone_index].begin(), zone_processed_[zone_index].end(), false);
      std::fill(zone_prev_node_[zone_index].begin(), zone_prev_node_[zone_index].end(), -1);
      std::fill(zone_used_link_[zone_index].begin(), zone_used_link_[zone_index].end(), -1);
      pq.push({0, zone_index});
      zone_distances_[zone_index][zone_index] = 0;
      while (!pq.empty()) {
        int cur_node = pq.top().second;
        pq.pop();
        if (zone_processed_[zone_index][cur_node]) {
          continue;
        }
        zone_processed_[zone_index][cur_node] = true;
        if (cur_node != zone_index) {
          T flow_amount = k_ * node_dest_demand_[cur_node][zone_index];
          node_dest_demand_[cur_node][zone_index] -= flow_amount;
          node_dest_demand_[zone_prev_node_[zone_index][cur_node]][zone_index] += flow_amount;
          this->links_[zone_used_link_[zone_index][cur_node]].flow += flow_amount;
          zone_distances_[zone_index][cur_node] = zone_distances_[zone_index][zone_prev_node_[zone_index][cur_node]] +
           this->links_[zone_used_link_[zone_index][cur_node]].Delay();
        }
        for (auto link_index : this->reverse_adjacency_list_[cur_node]) {
          int next_node = this->links_[link_index].init;
          if (zone_processed_[zone_index][next_node]) {
            continue;
          }  
          if (zone_distances_[zone_index][next_node] > zone_distances_[zone_index][cur_node] + this->links_[link_index].Delay()) {
            zone_distances_[zone_index][next_node] = zone_distances_[zone_index][cur_node] + this->links_[link_index].Delay();
            zone_prev_node_[zone_index][next_node] = cur_node;
            zone_used_link_[zone_index][next_node] = link_index;
            pq.push({zone_distances_[zone_index][next_node], next_node});
          }
        }
      }
    }

    void FlowPumpOut(T gamma) {
      for (int link_index = 0; link_index < this->number_of_links_; link_index++) {
        this->links_[link_index].flow *= (1 - gamma);
      }
      for (int zone_index = 0; zone_index < this->number_of_zones_; zone_index++) {
        for (int node_index = 0; node_index < this->number_of_nodes_; node_index++) {
          node_dest_demand_[zone_index][node_index] *= (1 - gamma);
        }
      }
      for (auto odp : this->origin_destination_pairs_) {
          auto origin_dest = odp.GetOriginDestination();
          node_dest_demand_[origin_dest.second][origin_dest.first] += gamma * odp.GetDemand();
        }
    }

    void RecordStatistics() {
      this->statistics_recorder_.PauseRecording();
      std::vector <std::vector <T>> node_dest_demand_copy(node_dest_demand_);
      std::vector <T> links_flow_copy(this->number_of_links_);
      for (int i = 0; i < this->number_of_links_; i++) {
        links_flow_copy[i] = this->links_[i].flow;
      }
      FullFlowTransfer();
      this->statistics_recorder_.ContinueRecording();
      this->statistics_recorder_.RecordStatistics();
    }

    void ZoneFullFlowTransfer(int zone_index) {
      std::priority_queue <std::pair <T, int>, std::vector <std::pair <T, int>>, std::greater <std::pair <T, int>>> pq;
      std::stack <int> processing_order;
      std::fill(zone_distances_[zone_index].begin(), zone_distances_[zone_index].end(), INF);
      std::fill(zone_processed_[zone_index].begin(), zone_processed_[zone_index].end(), false);
      std::fill(zone_prev_node_[zone_index].begin(), zone_prev_node_[zone_index].end(), -1);
      std::fill(zone_used_link_[zone_index].begin(), zone_used_link_[zone_index].end(), -1);
      pq.push({0, zone_index});
      zone_distances_[zone_index][zone_index] = 0;
      while (!pq.empty()) {
        int cur_node = pq.top().second;
        pq.pop();
        if (zone_processed_[zone_index][cur_node]) {
          continue;
        }
        zone_processed_[zone_index][cur_node] = true;
        if (cur_node != zone_index) {
          processing_order.push(cur_node);
        }
        for (auto link_index : this->reverse_adjacency_list_[cur_node]) {
          int next_node = this->links_[link_index].init;
          if (zone_processed_[zone_index][next_node]) {
            continue;
          } 
          //std::cout << zone_distances_[zone_index][next_node] << "!!!!" << zone_distances_[zone_index][cur_node] + this->links_[link_index].Delay() << '\n';  
          if (zone_distances_[zone_index][next_node] > zone_distances_[zone_index][cur_node] + this->links_[link_index].Delay()) {
            zone_distances_[zone_index][next_node] = zone_distances_[zone_index][cur_node] + this->links_[link_index].Delay();
            zone_prev_node_[zone_index][next_node] = cur_node;
            zone_used_link_[zone_index][next_node] = link_index;
            pq.push({zone_distances_[zone_index][next_node], next_node});
          }
        }
      }
      while (!processing_order.empty()) {
        int cur_node = processing_order.top();
        processing_order.pop();
        node_dest_demand_[zone_prev_node_[zone_index][cur_node]][zone_index] += node_dest_demand_[cur_node][zone_index];
        this->links_[zone_used_link_[zone_index][cur_node]].flow += node_dest_demand_[cur_node][zone_index];
        node_dest_demand_[cur_node][zone_index] = 0;
      }
    }

    void FullFlowTransfer() {
      for (int zone_index = 0; zone_index < this->number_of_zones_; zone_index++) {
        ZoneFullFlowTransfer(zone_index);
      }
    }

    void FinishTransfer() {
      int finish_iteration = 100;
      while (finish_iteration--) {
        for (int zone_index = 0; zone_index < this->number_of_zones_; zone_index++) {
            ZoneProcessing(zone_index, k_ / 3);
        }
      }
      FullFlowTransfer();
    }
    /*void DelayEstimation() {
      for (int node_index = 0; node_index < this->number_of_nodes_; node_index++) {
        for (auto link_index : this->adjacency_list_[node_index]) {
          for (int dest_index = 0; dest_index < this->number_of_zones_; dest_index++) {
            if (node_dest_delay_estimation_[node_index][dest_index] == -1 &&
                node_dest_delay_estimation_[this->links_[link_index].term][dest_index] != -1) {
              node_dest_delay_estimation_[node_index][dest_index] =
               node_dest_delay_estimation_[this->links_[link_index].term][dest_index] + this->links_[link_index].Delay();
            }
            else {
              if (node_dest_delay_estimation_[this->links_[link_index].term][dest_index] != -1) {
                node_dest_delay_estimation_[this->links_[link_index].term][dest_index] + this->links_[link_index].Delay();
              }
            }
          } 
        }
      }
    }*/

  };
}
#endif
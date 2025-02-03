#ifndef DEMAND_BASED_APPROACH_H
#define DEMAND_BASED_APPROACH_H

#include "TapasApproach.h"
#include <limits>

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
      int iteration_count = 0;
      while (iteration_count++ < 1000) {
        for (int zone_index = 0; zone_index < this->number_of_zones_; zone_index++) {
          ZoneProcessing(zone_index);
        }
        if (iteration_count % 5 == 0) {

        }
        FlowPumpOut(1.0 / (iteration_count + 3));
      }
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
    const int INF = std::numeric_limits<T>::max();
    const T k_ = 0.3;

    void ZoneProcessing(int zone_index) {
      std::priority_queue <std::pair <T, int>, std::vector <std::pair <T, int>>, std::greater <std::pair <T, int>>> pq;
      std::fill(zone_distances_[zone_index].begin(), zone_distances_[zone_index].end(), INF);
      std::fill(zone_processed_[zone_index].begin(), zone_processed_[zone_index].end(), false);
      std::fill(zone_prev_node_[zone_index].begin(), zone_prev_node_[zone_index].end(), -1);
      std::fill(zone_used_link_[zone_index].begin(), zone_used_link_[zone_index].end(), -1);
      pq.push({0, zone_index});
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
        }

        zone_distances_[zone_index][cur_node] = zone_distances_[zone_index][zone_prev_node_[zone_index][cur_node]] +
         this->links_[zone_used_link_[zone_index][cur_node]].Delay();
        

        for (auto link_index : this->reverse_adjacency_list_[cur_node]) {
          int next_node = this->links_[link_index].init;
          if (zone_processed_[zone_index][next_node]) {
            continue;
          }   
          if (zone_distances_[zone_index][next_node] > zone_distances_[zone_index][cur_node] + this->links_[link_index].Delay()) {
            zone_distances_[zone_index][next_node] = zone_distances_[zone_index][cur_node] + this->links_[link_index].Delay();
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
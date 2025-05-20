#ifndef DEMAND_BASED_APPROACH_H
#define DEMAND_BASED_APPROACH_H

#include "TapasApproach.h"
#include <limits>
#include <stack>

namespace TrafficAssignment {

  template <typename T>
  class PumpOutDemandBasedApproach : public TrafficAssignmentApproach <T> {
  public:
    PumpOutDemandBasedApproach(std::string dataset_name, T alpha = 1e-14) :
      TrafficAssignmentApproach<T>::TrafficAssignmentApproach(dataset_name, alpha) {
        zone_node_demand_.resize(this->number_of_zones_, std::vector <T> (this->number_of_nodes_, 0));
        //node_dest_demand_.resize(this->number_of_zones_, std::vector <T> (this->number_of_zones_, 0));
        //node_dest_delay_estimation_.resize(this->number_of_nodes_, std::vector <T> (this->number_of_zones_, -1));
        zone_link_flow_.resize(this->number_of_zones_, std::vector <T> (this->number_of_links_, 0));
        zone_distances_.resize(this->number_of_zones_, std::vector <T> (this->number_of_nodes_, INF));
        zone_processed_.resize(this->number_of_zones_, std::vector <bool> (this->number_of_nodes_, false));
        zone_prev_node_.resize(this->number_of_zones_, std::vector <int> (this->number_of_nodes_, 0));
        zone_used_link_.resize(this->number_of_zones_, std::vector <int> (this->number_of_nodes_, 0));
        for (auto odp : this->origin_destination_pairs_) {
          auto origin_dest = odp.GetOriginDestination();
          zone_node_demand_[origin_dest.second][origin_dest.first] = odp.GetDemand();
          //node_dest_demand_[origin_dest.second][origin_dest.first] = odp.GetDemand();
        }
        //std::cout << "Const_Finish\n";
      }

    ~PumpOutDemandBasedApproach() {
    }

    void ComputeTrafficFlows() override {
      this->statistics_recorder_.StartRecording(this, this->dataset_name_);
      int iteration_count = 0;
      while (iteration_count++ < number_of_iterations_) {
        //std::cout << iteration_count << '\n';
        for (int zone_index = 0; zone_index < this->number_of_zones_; zone_index++) {
          ZoneProcessing(zone_index);
        }
        //std::cout << "Statistics1\n";
        RecordStatistics();
        //std::cout << "Statistics2\n";
        FlowPumpOut(1.0 / (iteration_count + 2));
      }
      FinishTransfer();
      //FullFlowTransfer();
      //std::cout << 12345 << ' ' <<  std::setprecision(20) << this->ObjectiveFunction() << ' ' << this->RelativeGap() << '\n';
    }

    std::string GetApproachName() override {
      return "DemandBased";
    }

  protected:
  
    // Stores the demand for each node to each zone
    // The outer vector represents zones and the inner vector represents nodes.
    // Information is stored in a reversed order for uniform format
    std::vector <std::vector <T>> zone_node_demand_;

    // Stores the flow distribution on each link corresponding to each zone
    // The outer vector represents zones and the inner vector represents links.
    std::vector <std::vector <T>> zone_link_flow_;

    // Stores the links delay estimation from each node to each zone
    // The outer vector represents zones and the inner vector represents nodes.
    // Information is stored in a reversed order for uniform format
    std::vector <std::vector <T>> zone_distances_;

    // Stores information if node was processed during outer zone processing
    // The outer vector represents zones and the inner vector represents nodes.
    std::vector <std::vector <bool>> zone_processed_;

    // Stores the node that was the predecessor for a node during outer zone processing
    // The outer vector represents zones and the inner vector represents nodes.
    std::vector <std::vector <int>> zone_prev_node_;

    // Stores a link that was the used to get to a node during outer zone processing
    // The outer vector represents zones and the inner vector represents nodes.
    std::vector <std::vector <int>> zone_used_link_;

    const int number_of_iterations_ = 1e3;
    const T computation_threshold_ = 1e-14;
    //const T INF = std::numeric_limits<T>::max();
    const T INF = 1e9;
    const T pass_forward_processing_coefficient_ = 0.5;
    const int node_iteration_processing_count_ = 10;
    
    void NodeProcessing(int zone_index, int node_index, T pass_forward_coefficient, int node_iteration_count) {
      if (node_index == zone_index) {
        return;
      }
      int mn_node_index = zone_prev_node_[zone_index][node_index];
      int mn_link_index = zone_used_link_[zone_index][node_index];
      T flow_amount = zone_node_demand_[zone_index][node_index] * pass_forward_coefficient;
      //T flow_amount = node_dest_demand_[zone_index][node_index] * pass_forward_coefficient;
      zone_node_demand_[zone_index][node_index] -= flow_amount;
      //node_dest_demand_[zone_index][node_index] -= flow_amount;
      zone_node_demand_[zone_index][mn_node_index] += flow_amount;
      //node_dest_demand_[zone_index][mn_node_index] += flow_amount;
      this->links_[mn_link_index].flow += flow_amount;
      zone_link_flow_[zone_index][mn_link_index] += flow_amount;
      for (int cnt = 0; cnt < node_iteration_count; cnt++) {
        NodeFlowShift(zone_index, node_index);
      }
      zone_distances_[zone_index][node_index] = zone_distances_[zone_index][mn_node_index] + this->links_[mn_link_index].Delay();
    }

    void PerformNodeFlowShift(int zone_index, int node_index, T delta, int mn_link_index, int mx_link_index,
      int mn_node_index, int mx_node_index) {
      //std::cout << delta << ' ';
      delta = std::min(delta, zone_link_flow_[zone_index][mx_link_index]);
      delta = std::min(delta, zone_node_demand_[zone_index][mx_node_index]);
      //std::cout << delta << '\n';
      //delta = std::min(delta, node_dest_demand_[zone_index][mx_node_index]);
      this->links_[mx_link_index].flow -= delta;
      this->links_[mn_link_index].flow += delta;
      zone_link_flow_[zone_index][mx_link_index] -= delta;
      zone_link_flow_[zone_index][mn_link_index] += delta;
      zone_node_demand_[zone_index][mx_node_index] -= delta;
      //node_dest_demand_[zone_index][mx_node_index] -= delta;
      zone_node_demand_[zone_index][mn_node_index] += delta;
      //node_dest_demand_[zone_index][mn_node_index] += delta;
    }

    void NodeFlowShift(int zone_index, int node_index) {
      int mn_node_index = -1, mx_node_index = -1;
      int mn_link_index = 0, mx_link_index = 0;
      T mn_dist = 0, mx_dist = 0;
      for (int link_index : this->adjacency_list_[node_index]) {
        int transfer_node_index = this->links_[link_index].term;
        if (zone_processed_[zone_index][transfer_node_index]) { 
          T now_dist = zone_distances_[zone_index][transfer_node_index] + this->links_[link_index].Delay();
          if ((mn_node_index == -1) || (now_dist < mn_dist)) {
            mn_node_index = transfer_node_index;
            mn_link_index = link_index;
            mn_dist = now_dist;
          }
          if ((zone_link_flow_[zone_index][link_index] > computation_threshold_) && ((mx_node_index == -1) || (now_dist > mx_dist))) {
            mx_node_index = transfer_node_index;
            mx_link_index = link_index;
            mx_dist = now_dist;
          }
        }
      }
      if ((mx_node_index == -1) || (mn_node_index == mx_node_index)) {
        return;
      }
      T delta = (mx_dist - mx_dist) / (this->links_[mx_link_index].DelayDer() + this->links_[mn_link_index].DelayDer());
      //T delta = std::min(T(1.0), this->links_[mx_link_index].flow);
      //while ((delta < this->links_[mx_link_index].flow) &&
      //  (zone_distances_[zone_index][mx_node_index] + this->links_[mx_link_index].Delay(this->links_[mx_link_index].flow - delta) >
      //  zone_distances_[zone_index][mn_node_index] + this->links_[mn_link_index].Delay(this->links_[mn_link_index].flow + delta))) {
      //  delta *= 2;
      //}
      //while (zone_distances_[zone_index][mx_node_index] + this->links_[mx_link_index].Delay(this->links_[mx_link_index].flow - delta) <
      //zone_distances_[zone_index][mn_node_index] + this->links_[mn_link_index].Delay(this->links_[mn_link_index].flow + delta)) {
      //  delta /= 2;
      //}
      PerformNodeFlowShift(zone_index, node_index, delta, mn_link_index, mx_link_index, mn_node_index, mx_node_index);
    }

    void ZoneProcessing(int zone_index) {
      //std::cout << "Processing " << 1 << '\n';
      std::priority_queue <std::pair <T, int>, std::vector <std::pair <T, int>>, std::greater <std::pair <T, int>>> pq;
      for (int i = 0; i < this->number_of_nodes_; i++) {
        zone_distances_[zone_index][i] = INF;
        zone_processed_[zone_index][i] = false;
        zone_prev_node_[zone_index][i] = -1;
        zone_used_link_[zone_index][i] = -1;
      }

      //std::fill(zone_distances_[zone_index].begin(), zone_distances_[zone_index].end(), INF);
      //std::fill(zone_processed_[zone_index].begin(), zone_processed_[zone_index].end(), false);
      //std::fill(zone_prev_node_[zone_index].begin(), zone_prev_node_[zone_index].end(), -1);
      //std::fill(zone_used_link_[zone_index].begin(), zone_used_link_[zone_index].end(), -1);
      pq.push({0, zone_index});
      //std::cout << "Processing " << 2 << '\n';
      zone_distances_[zone_index][zone_index] = 0;
      while (!pq.empty()) {
        //std::cout << "Processing " << 21 << '\n';
        int cur_node = pq.top().second;
        pq.pop();
        if (zone_processed_[zone_index][cur_node]) {
          continue;
        }
        zone_processed_[zone_index][cur_node] = true;
        if (cur_node != zone_index) {
          /*T flow_amount = k_ * node_dest_demand_[cur_node][zone_index];
          node_dest_demand_[cur_node][zone_index] -= flow_amount;
          node_dest_demand_[zone_prev_node_[zone_index][cur_node]][zone_index] += flow_amount;
          this->links_[zone_used_link_[zone_index][cur_node]].flow += flow_amount;
          zone_distances_[zone_index][cur_node] = zone_distances_[zone_index][zone_prev_node_[zone_index][cur_node]] +
           this->links_[zone_used_link_[zone_index][cur_node]].Delay();
           */
          //std::cout << 1 << '\n';
          NodeProcessing(zone_index, cur_node, pass_forward_processing_coefficient_, node_iteration_processing_count_);
          //std::cout << 2 << '\n';
        }
        //std::cout << "Processing " << 22 << '\n';
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
        //std::cout << "Processing " << 23 << '\n';
      }
      //std::cout << "Processing " << 3 << '\n';
    }

    void FlowPumpOut(T gamma) {
      for (int link_index = 0; link_index < this->number_of_links_; link_index++) {
        this->links_[link_index].flow *= (1 - gamma);
      }
      for (int zone_index = 0; zone_index < this->number_of_zones_; zone_index++) {
        for (int node_index = 0; node_index < this->number_of_nodes_; node_index++) {
          zone_node_demand_[zone_index][node_index] *= (1 - gamma);
          //node_dest_demand_[zone_index][node_index] *= (1 - gamma);
        }
        for (int link_index = 0; link_index < this->number_of_links_; link_index++) {
          zone_link_flow_[zone_index][link_index] *= (1 - gamma);
        }
      }
      for (auto odp : this->origin_destination_pairs_) {
          auto origin_dest = odp.GetOriginDestination();
          zone_node_demand_[origin_dest.second][origin_dest.first] += gamma * odp.GetDemand();
          //node_dest_demand_[origin_dest.second][origin_dest.first] += gamma * odp.GetDemand();
        }
    }

    void RecordStatistics() {
      this->statistics_recorder_.PauseRecording();
      std::vector <std::vector <T>> zone_node_demand_copy(zone_node_demand_);
      std::vector <T> links_flow_copy(this->number_of_links_);
      for (int i = 0; i < this->number_of_links_; i++) {
        links_flow_copy[i] = this->links_[i].flow;
      }
      FullFlowTransfer();
      this->statistics_recorder_.ContinueRecording();
      //std::cout << 1 << '\n';
      this->statistics_recorder_.RecordStatistics();

      this->statistics_recorder_.PauseRecording();
      zone_node_demand_ = zone_node_demand_copy;
      for (int i = 0; i < this->number_of_links_; i++) {
        this->links_[i].flow = links_flow_copy[i];
      }
      this->statistics_recorder_.ContinueRecording();
    }

    void ZoneFullFlowTransfer(int zone_index) {
      std::priority_queue <std::pair <T, int>, std::vector <std::pair <T, int>>, std::greater <std::pair <T, int>>> pq;
      std::stack <int> processing_order;
      for (int i = 0; i < this->number_of_nodes_; i++) {
        zone_distances_[zone_index][i] = INF;
        zone_processed_[zone_index][i] = false;
        zone_prev_node_[zone_index][i] = -1;
        zone_used_link_[zone_index][i] = -1;
      }
      //std::fill(zone_distances_[zone_index].begin(), zone_distances_[zone_index].end(), INF);
      //std::fill(zone_processed_[zone_index].begin(), zone_processed_[zone_index].end(), false);
      //std::fill(zone_prev_node_[zone_index].begin(), zone_prev_node_[zone_index].end(), -1);
      //std::fill(zone_used_link_[zone_index].begin(), zone_used_link_[zone_index].end(), -1);
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
        zone_node_demand_[zone_index][zone_prev_node_[zone_index][cur_node]] += zone_node_demand_[zone_index][cur_node];
        //node_dest_demand_[zone_prev_node_[zone_index][cur_node]][zone_index] += node_dest_demand_[cur_node][zone_index];
        this->links_[zone_used_link_[zone_index][cur_node]].flow += zone_node_demand_[zone_index][cur_node];
        zone_node_demand_[zone_index][cur_node] = 0;
        //node_dest_demand_[cur_node][zone_index] = 0;
      }
    }

    void FullFlowTransfer() {
      for (int zone_index = 0; zone_index < this->number_of_zones_; zone_index++) {
        ZoneFullFlowTransfer(zone_index);
      }
    }

    void FinishTransfer() {
      //std::cout << 1 << '\n';
      int finish_iteration = 100;
      while (finish_iteration--) {
        for (int zone_index = 0; zone_index < this->number_of_zones_; zone_index++) {
            ZoneProcessing(zone_index);
        }
      }
      FullFlowTransfer();
      //std::cout << 2 << '\n';
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
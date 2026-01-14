#ifndef DEMAND_BASED_APPROACH_H
#define DEMAND_BASED_APPROACH_H

#include "../common/TrafficAssignmentApproach.h"
#include <limits>
#include <stack>

namespace TrafficAssignment {

  template <typename T>
  class PumpOutDemandBasedApproach : public TrafficAssignmentApproach <T> {
  public:
    PumpOutDemandBasedApproach(Network<T>& network, T alpha = 1e-14) :
      TrafficAssignmentApproach<T>(network, alpha)  {
        InitializeFromNetwork();
    }



    ~PumpOutDemandBasedApproach() {
    }

    void ComputeTrafficFlows() override {
      this->statistics_recorder_.StartRecording(this->GetApproachName());
      int iteration_count = 0;
      while (iteration_count++ < number_of_iterations_) {
        //std::cout << iteration_count << '\n';
        for (int zone_index = 0; zone_index < this->number_of_zones_; zone_index++) {
          ZoneProcessing(zone_index);
        }
        RecordStatistics();
        FlowPumpOut(1.0 / (iteration_count + 2));
      }
      FinishTransfer();
    }

    std::string GetApproachName() override {
      return "DemandBased";
    }

  protected:
    // Network dimensions
    int number_of_zones_;
    int number_of_nodes_;
    int number_of_links_;

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
    
    void InitializeFromNetwork() {
        number_of_zones_ = this->network_.number_of_zones();
        number_of_nodes_ = this->network_.number_of_nodes();
        number_of_links_ = this->network_.number_of_links();

        zone_node_demand_.resize(number_of_zones_, std::vector<T>(number_of_nodes_, 0));
        zone_link_flow_.resize(number_of_zones_, std::vector<T>(number_of_links_, 0));
        zone_distances_.resize(number_of_zones_, std::vector<T>(number_of_nodes_, INF));
        zone_processed_.resize(number_of_zones_, std::vector<bool>(number_of_nodes_, false));
        zone_prev_node_.resize(number_of_zones_, std::vector<int>(number_of_nodes_, 0));
        zone_used_link_.resize(number_of_zones_, std::vector<int>(number_of_nodes_, 0));

        for (const auto& odp : this->network_.od_pairs()) {
            int origin = odp.GetOriginDestination().first;
            int dest = odp.GetOriginDestination().second;
            zone_node_demand_[dest][origin] = odp.GetDemand();
        }
    }

    void NodeProcessing(int zone_index, int node_index, T pass_forward_coefficient, int node_iteration_count) {
      if (node_index == zone_index) {
        return;
      }
      int mn_node_index = zone_prev_node_[zone_index][node_index];
      int mn_link_index = zone_used_link_[zone_index][node_index];
      T flow_amount = zone_node_demand_[zone_index][node_index] * pass_forward_coefficient;
      zone_node_demand_[zone_index][node_index] -= flow_amount;
      zone_node_demand_[zone_index][mn_node_index] += flow_amount;
      this->network_.SetLinkFlow(mn_link_index, flow_amount);

      zone_link_flow_[zone_index][mn_link_index] += flow_amount;
      for (int cnt = 0; cnt < node_iteration_count; cnt++) {
        NodeFlowShift(zone_index, node_index);
      }
      zone_distances_[zone_index][node_index] = zone_distances_[zone_index][mn_node_index] + this->network_.links()[mn_link_index].Delay();
    }

    void PerformNodeFlowShift(int zone_index, int node_index, T delta, int mn_link_index, int mx_link_index,
      int mn_node_index, int mx_node_index) {
      delta = std::min(delta, zone_link_flow_[zone_index][mx_link_index]);
      delta = std::min(delta, zone_node_demand_[zone_index][mx_node_index]);
      this->network_.SetLinkFlow(mx_link_index, -delta);
      this->network_.SetLinkFlow(mn_link_index, delta);

      zone_link_flow_[zone_index][mx_link_index] -= delta;
      zone_link_flow_[zone_index][mn_link_index] += delta;
      zone_node_demand_[zone_index][mx_node_index] -= delta;
      zone_node_demand_[zone_index][mn_node_index] += delta;
    }

    void NodeFlowShift(int zone_index, int node_index) {
      int mn_node_index = -1, mx_node_index = -1;
      int mn_link_index = 0, mx_link_index = 0;
      T mn_dist = 0, mx_dist = 0;
      for (const int link_index : this->network_.adjacency()[node_index]) {
        int transfer_node_index = this->network_.links()[link_index].term;
        if (zone_processed_[zone_index][transfer_node_index]) { 
          T now_dist = zone_distances_[zone_index][transfer_node_index] + this->network_.links()[link_index].Delay();
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
      T delta = (mx_dist - mx_dist) / (this->network_.links()[mx_link_index].DelayDer() + this->network_.links()[mn_link_index].DelayDer());
      PerformNodeFlowShift(zone_index, node_index, delta, mn_link_index, mx_link_index, mn_node_index, mx_node_index);
    }

    void ZoneProcessing(int zone_index) {
      std::priority_queue <std::pair <T, int>, std::vector <std::pair <T, int>>, std::greater <std::pair <T, int>>> pq;
      for (int i = 0; i < this->number_of_nodes_; i++) {
        zone_distances_[zone_index][i] = INF;
        zone_processed_[zone_index][i] = false;
        zone_prev_node_[zone_index][i] = -1;
        zone_used_link_[zone_index][i] = -1;
      }

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
          NodeProcessing(zone_index, cur_node, pass_forward_processing_coefficient_, node_iteration_processing_count_);
        }
        for (auto link_index : this->network_.reverse_adjacency()[cur_node]) {
          int next_node = this->network_.links()[link_index].init;
          if (zone_processed_[zone_index][next_node]) {
            continue;
          }  
          if (zone_distances_[zone_index][next_node] > zone_distances_[zone_index][cur_node] + this->network_.links()[link_index].Delay()) {
            zone_distances_[zone_index][next_node] = zone_distances_[zone_index][cur_node] + this->network_.links()[link_index].Delay();
            zone_prev_node_[zone_index][next_node] = cur_node;
            zone_used_link_[zone_index][next_node] = link_index;
            pq.push({zone_distances_[zone_index][next_node], next_node});
          }
        }
      }
    }

    void FlowPumpOut(T gamma) {
      for (int link_index = 0; link_index < this->number_of_links_; link_index++) {
        auto cur_flow = this->network_.links()[link_index].flow;
        this->network_.SetLinkFlow(link_index, -cur_flow);
        this->network_.SetLinkFlow(link_index, cur_flow * (1 - gamma));
      }
      for (int zone_index = 0; zone_index < this->number_of_zones_; zone_index++) {
        for (int node_index = 0; node_index < this->number_of_nodes_; node_index++) {
          zone_node_demand_[zone_index][node_index] *= (1 - gamma);
        }
        for (int link_index = 0; link_index < this->number_of_links_; link_index++) {
          zone_link_flow_[zone_index][link_index] *= (1 - gamma);
        }
      }
      for (auto odp : this->network_.od_pairs()) {
          auto origin_dest = odp.GetOriginDestination();
          zone_node_demand_[origin_dest.second][origin_dest.first] += gamma * odp.GetDemand();
        }
    }

    void RecordStatistics() {
      this->statistics_recorder_.PauseRecording();
      std::vector <std::vector <T>> zone_node_demand_copy(zone_node_demand_);
      std::vector <T> links_flow_copy(this->number_of_links_);
      for (int i = 0; i < this->number_of_links_; i++) {
        
        links_flow_copy[i] = this->network_.links()[i].flow;;
      }
      FullFlowTransfer();
      this->statistics_recorder_.ResumeRecording();
      this->statistics_recorder_.RecordStatistics();

      this->statistics_recorder_.PauseRecording();
      zone_node_demand_ = zone_node_demand_copy;
      for (int link_index = 0; link_index < this->number_of_links_; link_index++) {
        auto cur_flow = this->network_.links()[link_index].flow;
        this->network_.SetLinkFlow(link_index, -cur_flow);
        this->network_.SetLinkFlow(link_index, links_flow_copy[link_index]);  
      }
      this->statistics_recorder_.ResumeRecording();
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
        for (auto link_index : this->network_.reverse_adjacency()[cur_node]) {
          int next_node = this->network_.links()[link_index].init;
          if (zone_processed_[zone_index][next_node]) {
            continue;
          }   
          if (zone_distances_[zone_index][next_node] > zone_distances_[zone_index][cur_node] +this->network_.links()[link_index].Delay()) {
            zone_distances_[zone_index][next_node] = zone_distances_[zone_index][cur_node] + this->network_.links()[link_index].Delay();
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
        this->network_.SetLinkFlow(zone_used_link_[zone_index][cur_node], zone_node_demand_[zone_index][cur_node]);
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
      int finish_iteration = 100;
      while (finish_iteration--) {
        for (int zone_index = 0; zone_index < this->number_of_zones_; zone_index++) {
            ZoneProcessing(zone_index);
        }
      }
      FullFlowTransfer();
    }

  };
}
#endif
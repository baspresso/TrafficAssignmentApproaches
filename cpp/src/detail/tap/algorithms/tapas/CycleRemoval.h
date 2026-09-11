#ifndef CYCLE_REMOVAL_H
#define CYCLE_REMOVAL_H

#include <vector>
#include <algorithm>
#include "../../core/Network.h"
#include "BushData.h"

namespace TrafficAssignment {

  template <typename T>
  void RemoveCycleFlow(OriginBush<T>& bush, Network<T>& network,
                       const std::vector<int>& cycle, T zero_flow) {
    T min_flow = bush.flows[cycle[0]];
    for (int link_id : cycle) {
      min_flow = std::min(min_flow, bush.flows[link_id]);
    }
    for (int link_id : cycle) {
      if (bush.flows[link_id] - min_flow < zero_flow)
        bush.flows[link_id] = T(0);
      else
        bush.flows[link_id] -= min_flow;
      T current_flow = network.links()[link_id].flow;
      if (current_flow - min_flow < zero_flow)
        network.SetLinkFlow(link_id, -current_flow);
      else
        network.SetLinkFlow(link_id, -min_flow);
    }
  }

  template <typename T>
  void DfsForCycleIdentification(OriginBush<T>& bush, Network<T>& network,
                                 int cur_node,
                                 std::vector<int>& link_used,
                                 std::vector<int>& state,
                                 bool& cycle_found, T zero_flow) {
    if (cycle_found || state[cur_node] != 0) return;
    state[cur_node] = 1;
    for (int now : network.adjacency()[cur_node]) {
      if (cycle_found) return;
      int next_node = network.links()[now].term;
      if (bush.flows[now] < zero_flow || state[next_node] == 2) continue;
      if (state[next_node] == 1) {
        cycle_found = true;
        std::vector<int> cycle;
        cycle.push_back(now);
        int temp = cur_node;
        while (temp != next_node) {
          cycle.push_back(link_used[temp]);
          temp = network.links()[link_used[temp]].init;
        }
        RemoveCycleFlow(bush, network, cycle, zero_flow);
        return;
      }
      if (state[next_node] == 0) {
        link_used[next_node] = now;
        DfsForCycleIdentification(bush, network, next_node, link_used, state, cycle_found, zero_flow);
      }
    }
    state[cur_node] = 2;
  }

  template <typename T>
  std::vector<int> GetBushNodes(const OriginBush<T>& bush, const Network<T>& network, T zero_flow) {
    std::vector<bool> is_bush_node(network.number_of_nodes(), false);
    for (int link_id = 0; link_id < network.number_of_links(); link_id++) {
      if (bush.flows[link_id] >= zero_flow) {
        is_bush_node[network.links()[link_id].init] = true;
        is_bush_node[network.links()[link_id].term] = true;
      }
    }
    std::vector<int> nodes;
    for (int i = 0; i < network.number_of_nodes(); i++)
      if (is_bush_node[i]) nodes.push_back(i);
    return nodes;
  }

  template <typename T>
  void RemoveCyclicFlows(OriginBush<T>& bush, Network<T>& network, T zero_flow) {
    bool cycle_found = true;
    while (cycle_found) {
      cycle_found = false;
      auto bush_nodes = GetBushNodes(bush, network, zero_flow);
      std::vector<int> link_used(network.number_of_nodes(), -1);
      std::vector<int> state(network.number_of_nodes(), 0);
      for (int node : bush_nodes) {
        if (state[node] == 0) {
          DfsForCycleIdentification(bush, network, node, link_used, state, cycle_found, zero_flow);
          if (cycle_found) break;
        }
      }
    }
  }

}  // namespace TrafficAssignment

#endif  // CYCLE_REMOVAL_H

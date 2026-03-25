#ifndef BUSH_DATA_H
#define BUSH_DATA_H

#include <vector>
#include <algorithm>

namespace TrafficAssignment {

  template <typename T>
  struct OriginBush {
    std::vector<T> flows;              // [link_id]
    std::vector<int> sp_parent;        // [node_id] = parent node
    std::vector<T> sp_distance;        // [node_id] = SP distance
    std::vector<bool> sp_tree_links;   // [link_id] = in SP tree?

    void Initialize(int num_links, int num_nodes) {
      flows.assign(num_links, T(0));
      sp_parent.assign(num_nodes, -1);
      sp_distance.assign(num_nodes, T(0));
      sp_tree_links.assign(num_links, false);
    }

    void Clear() {
      std::fill(flows.begin(), flows.end(), T(0));
      std::fill(sp_parent.begin(), sp_parent.end(), -1);
      std::fill(sp_distance.begin(), sp_distance.end(), T(0));
      std::fill(sp_tree_links.begin(), sp_tree_links.end(), false);
    }
  };

}  // namespace TrafficAssignment

#endif  // BUSH_DATA_H

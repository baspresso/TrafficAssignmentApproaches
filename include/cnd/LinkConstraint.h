#ifndef CND_LINK_CONSTRAINT_H
#define CND_LINK_CONSTRAINT_H

#include <cstddef>
#include <string>

namespace TrafficAssignment {

/**
 * @class DirectedLinkCapacityConstraint
 * @brief Represents a constraint for a directed link (arc) in the network
 *
 * This structure represents capacity constraints for directed arcs,
 * identified by their origin (init_node) and destination (term_node) nodes.
 *
 * Lives in its own header (rather than DirectedConstraintLoader.h) so that
 * lightweight consumers such as CndStatisticsRecorder can use the complete
 * type without pulling in the loader or creating an include cycle.
 */
struct DirectedLinkCapacityConstraint {
  std::size_t init_node;               // Origin node
  std::size_t term_node;               // Destination node
  double lower_bound;                  // Minimum capacity
  double upper_bound;                  // Maximum capacity
  double investment_cost_param;        // Cost per unit capacity

  DirectedLinkCapacityConstraint() = default;

  DirectedLinkCapacityConstraint(
    std::size_t origin,
    std::size_t destination,
    double lb,
    double ub,
    double cost)
    : init_node(origin),
      term_node(destination),
      lower_bound(lb),
      upper_bound(ub),
      investment_cost_param(cost) {}

  /**
   * @brief Get string representation for logging
   */
  std::string ToString() const {
    return std::string("(") + std::to_string(init_node) + "->" +
           std::to_string(term_node) + ")";
  }

  /**
   * @brief Get unique arc ID (for compatibility with undirected networks)
   * Format: "init_node_term_node"
   */
  std::string GetArcID() const {
    return std::to_string(init_node) + "_" + std::to_string(term_node);
  }
};

}  // namespace TrafficAssignment

#endif  // CND_LINK_CONSTRAINT_H

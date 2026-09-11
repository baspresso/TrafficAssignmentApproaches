#pragma once

#include <traffic_assignment/options.hpp>

#include <string>
#include <vector>

namespace traffic_assignment {

/// Independent values captured after a solve, without references to live state.
/// A returned result reports termination; it is not an optimality certificate.
struct TapResult {
  std::string approach;
  std::string network_name;
  Real relative_gap = 0;
  Real total_travel_time = 0;
  Real beckmann_objective = 0;
  std::vector<Real> flows;
  std::vector<Real> link_costs;
  double solve_seconds = 0;
};

struct CndpResult {
  double objective = 0;
  double total_travel_time = 0;
  double budget = 0;
  double budget_upper_bound = 0;
  std::vector<Real> capacities;
  std::vector<Real> flows;
  std::vector<double> lower_bounds;
  std::vector<double> upper_bounds;
  double elapsed_seconds = 0;
};

}  // namespace traffic_assignment

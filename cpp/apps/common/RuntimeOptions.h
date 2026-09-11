#pragma once

#include "TomlConfigLoader.h"

namespace TrafficAssignment::Config {

// This adapter only translates CLI/TOML defaults; algorithm selection is native.
inline traffic_assignment::TapOptions ToTapOptions(const SolverConfig& config) {
  traffic_assignment::TapOptions options;
  options.approach = config.approach;
  options.relative_gap_tolerance = config.approach_alpha;
  options.shift_method = config.route_based.shift_method;
  options.route_search_threads = config.route_based.route_search_threads;
  if (config.max_standard_iterations > 0) options.max_iterations = config.max_standard_iterations;
  if (config.route_based.full_iteration_count > 0)
    options.full_iteration_count = config.route_based.full_iteration_count;
  if (config.route_based.origin_iteration_count > 0)
    options.origin_iteration_count = config.route_based.origin_iteration_count;
  if (config.route_based.ema_alpha > 0) options.ema_alpha = config.route_based.ema_alpha;
  if (config.tapas.mu > 0) options.mu = config.tapas.mu;
  if (config.tapas.v > 0) options.v = config.tapas.v;
  return options;
}

inline traffic_assignment::CndpOptions ToCndpOptions(const CndpConfig& config) {
  traffic_assignment::CndpOptions options;
  options.link_capacity_selection_threshold = config.solver.link_capacity_selection_threshold;
  options.budget_threshold = config.solver.budget_threshold;
  options.budget_function_multiplier = config.solver.budget_function_multiplier;
  options.budget_upper_bound = config.solver.budget_upper_bound;
  options.route_search_threads = config.solver.route_based.route_search_threads;
  options.metrics = config.metrics;
  options.progress_format = config.output.progress_format;
  options.statistics_enabled = true;
  options.final_diagnostics = true;
  return options;
}

}  // namespace TrafficAssignment::Config

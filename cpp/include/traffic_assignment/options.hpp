#pragma once

#include <cstddef>
#include <string>

namespace traffic_assignment {

// The supported native precision; optimizer interfaces may use double.
using Real = long double;

struct StepConfig {
  std::string type;           ///< Step type: "nlopt", "optimality_condition", or "optimlib".
  std::string name;           ///< Optional display name; auto-generated from type+algorithm if empty.
  int max_iterations = 100;   ///< Maximum iterations/evaluations for this step.
  double tolerance = 1e-4;    ///< Convergence tolerance (relative objective change).

  // NLopt-specific
  std::string algorithm;          ///< NLopt algorithm name (e.g., "LN_COBYLA", "GN_ISRES").
  std::string local_algorithm;    ///< Local optimizer for composite algorithms (AUGLAG, global+local).
  int local_max_iterations = 0;   ///< Max evaluations for local optimizer (0 = max_iterations/10).
  double local_tolerance = 0.0;   ///< Convergence tolerance for local optimizer (0 = use main tolerance).

  // Population-based (OptimLib) specific
  int population_size = 0;  ///< Population size for DE/PSO (0 = auto: max(200, 2*n_vars)).

  // Gradient descent specific
  double step_size = 1.0;    ///< Initial step size for gradient descent (Armijo starting alpha).
  double fd_epsilon = 1e-4;  ///< Finite difference perturbation size for gradient estimation.
  std::string gradient_method;  ///< Gradient estimator: "finite_difference" (default), "spsa", "sensitivity".
  std::string stochastic_optimizer;  ///< Optimizer for stochastic gradient: "sgd" (default), "momentum", "adam".
};

struct MetricsConfig {
  bool enable_trace = true;             ///< Write per-iteration trace CSV (objective vs time).
  bool enable_relative_gap = true;      ///< Sample TAP relative gap during trace recording.
  int relative_gap_sample_period = 10;  ///< Sample RGAP every N trace points (expensive to compute).
  int flush_every_n_points = 0;         ///< Flush trace to disk every N points (0 = flush at end).
  bool write_metadata_json = true;      ///< Write JSON file with algorithm config and run parameters.
  bool write_summary_csv = true;        ///< Append one-line summary to shared CSV (append-only).
  bool append_dataset_subdir = true;    ///< Append dataset name as subdirectory to output_root.
  std::string output_root = "performance_results"; ///< Root directory for all metrics output.
  std::string run_id;                   ///< Unique run identifier (auto-generated if empty).
  std::string scenario_name;            ///< Scenario label for grouping runs in analysis.
};

struct LinkConstraint {
  std::size_t init_node;               // Origin node
  std::size_t term_node;               // Destination node
  double lower_bound;                  // Minimum capacity
  double upper_bound;                  // Maximum capacity
  double investment_cost_param;        // Cost per unit capacity

  LinkConstraint() = default;

  LinkConstraint(
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

/// Options shared by the Python API and standalone TAP applications.
struct TapOptions {
  std::string approach = "tapas";
  Real relative_gap_tolerance = 1e-14L;
  int max_iterations = 200;
  Real mu = 0.5L;
  Real v = 0.25L;
  std::string shift_method = "NewtonStep";
  std::size_t route_search_threads = 1;
  int full_iteration_count = 3;
  int origin_iteration_count = 1;
  Real ema_alpha = 0.7L;
};

/// Output is opt-in for library callers. CLI adapters enable it explicitly.
struct CndpOptions {
  Real link_capacity_selection_threshold = 1e-3L;
  Real budget_threshold = 1e-1L;
  Real budget_function_multiplier = 5.0L;
  double budget_upper_bound = 100000.0;
  std::size_t route_search_threads = 1;
  bool exclude_insensitive_links = true;
  bool statistics_enabled = false;
  bool final_diagnostics = false;
  bool verbose = false;
  std::string progress_format = "none";
  MetricsConfig metrics;
};

}  // namespace traffic_assignment

#ifndef CND_OPTIMIZATION_CONTEXT_H
#define CND_OPTIMIZATION_CONTEXT_H

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <sstream>
#include <algorithm>
#include <limits>
#include <nlopt.hpp>
#include "../tap/algorithms/common/TrafficAssignmentApproach.h"
#include "../tap/core/Network.h"
#include "DirectedConstraintLoader.h"
#include "CndStatisticsRecorder.h"

namespace TrafficAssignment {

/**
 * @brief Shared evaluation context for bilevel CNDP optimization steps.
 *
 * Holds references to the network, TAP solver, and capacity constraints, plus shared
 * parameters (budget multiplier, bounds) and runtime counters. All optimization steps
 * operate through this context to evaluate the upper-level objective:
 *
 *   F(y) = TSTT(x(y)) + theta * sum(y_a - lb_a)
 *
 * where x(y) is the user-equilibrium flow for capacity vector y, and theta is the
 * budget function multiplier. Chiou (2005) Eq. 1.
 *
 * Provides crash recovery for TAP solver failures (retry with reset) and
 * budget feasibility enforcement via proportional projection.
 *
 * @tparam T Numeric type for flow/capacity computations.
 */
template <typename T>
class CndOptimizationContext {
public:
  // --- Named constants ---
  static constexpr double kDenominatorZeroGuard    = 1e-40;  ///< Guard against division by zero in Jacobian computations.
  static constexpr double kNegativeFlowTolerance   = -1e-8;  ///< Tolerance for detecting invalid negative flows.
  static constexpr double kCapacityBoundTolerance   = 1e-5;  ///< Tolerance for capacity-at-bound checks.
  static constexpr double kFlowSignificanceThreshold = 1e-3; ///< Minimum flow for a link to be considered active.
  static constexpr double kInitialTransferAmount    = 1.0;   ///< Starting capacity transfer in optimality iterations.
  static constexpr int    kProgressBarWidth         = 30;
  static constexpr int    kMaxInvalidStateWarnings  = 10;    ///< Maximum invalid-state warnings before suppression.
  static constexpr double kBinarySearchDivisor      = 2.0;   ///< Halving factor for transfer validation retries.
  static constexpr double kInvalidObjectivePenalty   = 1e30;  ///< Penalty returned when objective is non-finite.
  static constexpr double kDefaultBudgetViolationPenaltyFactor = 100.0; ///< Quadratic penalty coefficient for soft budget constraint.
  static constexpr double kConstraintTolerance      = 1e-4;  ///< NLopt constraint feasibility tolerance.

  /// @brief Console progress bar state for long-running optimization loops.
  struct ProgressState {
    int count = 0;
    int total = 0;
    std::chrono::steady_clock::time_point start_time;
    bool active = false;

    void Start(int total_target) {
      count = 0;
      total = total_target;
      start_time = std::chrono::steady_clock::now();
      active = total_target > 0;
    }

    void Finish() {
      if (!active) return;
      std::cout << '\n';
      active = false;
    }

    void Advance() { ++count; }
  };

  /**
   * @brief Mutable counters tracking iteration counts, TA solver timing, and crash recovery
   *        across all pipeline steps. Persists for the lifetime of one ComputeNetworkDesign() call.
   */
  struct RuntimeCounters {
    int standard_trace_step = 0;           ///< Cumulative evaluation count for standard (NLopt/OptimLib) steps.
    int optcond_trace_step = 0;            ///< Cumulative iteration count for optimality condition steps.
    int total_ta_run_count = 0;            ///< Total number of lower-level TAP solver invocations.
    double ta_compute_time_sum_seconds = 0.0; ///< Cumulative TA solver wall-clock time.
    double ta_compute_time_max_seconds = 0.0; ///< Maximum single TA solver invocation time.
    double max_budget_violation = 0.0;     ///< Worst observed budget violation max(0, B(y) - B_max).
    int invalid_ta_state_count = 0;        ///< Number of TA solver runs producing invalid state.
    int ta_recovery_success_count = 0;     ///< Successful recoveries via reset+retry.
    int ta_recovery_failure_count = 0;     ///< Failed recoveries (penalty objective returned).
    int invalid_value_log_count = 0;       ///< Warning counter for suppression after kMaxInvalidStateWarnings.

    void Reset() { *this = RuntimeCounters{}; }
  };

  // --- Non-owning references ---
  Network<T>& network;                        ///< Transportation network with mutable link capacities and flows.
  std::shared_ptr<TrafficAssignmentApproach<T>> approach; ///< Lower-level TAP solver (Beckmann UE).
  const std::vector<DirectedLinkCapacityConstraint>& constraints; ///< Per-link capacity bounds [lb_a, ub_a].

  // --- Shared parameters ---
  T budget_function_multiplier;               ///< theta: multiplier in budget cost theta * sum(y_a - lb_a). Chiou (2005) Eq. 2.
  double budget_upper_bound;                  ///< B: maximum allowable investment budget.
  T budget_threshold;                         ///< Tolerance for budget-at-bound detection in optimality iterations.
  T link_capacity_selection_threshold;        ///< Minimum distance from bounds to consider a link for capacity transfer.
  std::size_t route_search_thread_count;      ///< Parallel thread count for route enumeration in RouteBased approach.
  bool verbose;                               ///< Enable console output (progress bars, warnings).

  // --- Mutable state ---
  RuntimeCounters& counters;                  ///< Shared counters across all pipeline steps.
  CndStatisticsRecorder<T>& statistics_recorder; ///< Convergence trace and summary recorder.

  CndOptimizationContext(
    Network<T>& network_ref,
    std::shared_ptr<TrafficAssignmentApproach<T>> approach_ref,
    const std::vector<DirectedLinkCapacityConstraint>& constraints_ref,
    T budget_function_multiplier_val,
    double budget_upper_bound_val,
    T budget_threshold_val,
    T link_capacity_selection_threshold_val,
    std::size_t route_search_thread_count_val,
    bool verbose_val,
    RuntimeCounters& counters_ref,
    CndStatisticsRecorder<T>& statistics_recorder_ref
  )
    : network(network_ref),
      approach(approach_ref),
      constraints(constraints_ref),
      budget_function_multiplier(budget_function_multiplier_val),
      budget_upper_bound(budget_upper_bound_val),
      budget_threshold(budget_threshold_val),
      link_capacity_selection_threshold(link_capacity_selection_threshold_val),
      route_search_thread_count(route_search_thread_count_val),
      verbose(verbose_val),
      counters(counters_ref),
      statistics_recorder(statistics_recorder_ref) {}

  // --- Utility methods ---

  /// @brief Checks if a double value is finite (not NaN or infinity).
  static bool IsFiniteNumber(double value) {
    return std::isfinite(value);
  }

  bool IsFiniteScalar(T value) const {
    return IsFiniteNumber(static_cast<double>(value));
  }

  static T NaNValue() {
    return static_cast<T>(std::numeric_limits<double>::quiet_NaN());
  }

  /// @brief Computes investment cost from current network capacities: theta * sum(y_a - lb_a).
  ///        Chiou (2005) Eq. 2.
  T Budget() {
    T result = 0;
    for (int link_index = 0; link_index < network.number_of_links(); link_index++) {
      result += network.mutable_links()[link_index].capacity - constraints[link_index].lower_bound;
    }
    return result * budget_function_multiplier;
  }

  /// @brief Computes investment cost from explicit capacity vector x (NLopt callback interface).
  double BudgetFunction(int n, const double* x) {
    double result = 0;
    for (int i = 0; i < n; i++) {
      result += x[i] - constraints[i].lower_bound;
    }
    return result * budget_function_multiplier;
  }

  /// @brief Checks if any link has non-finite flow/capacity/delay or negative values.
  bool HasInvalidLinkState() const {
    for (const auto& link : network.links()) {
      const double flow = static_cast<double>(link.flow);
      const double capacity = static_cast<double>(link.capacity);
      const double delay = static_cast<double>(link.Delay());
      if (!IsFiniteNumber(flow) || !IsFiniteNumber(capacity) || !IsFiniteNumber(delay)) {
        return true;
      }
      if (capacity < 0.0 || delay < 0.0) {
        return true;
      }
      if (flow < kNegativeFlowTolerance) {
        return true;
      }
    }
    return false;
  }

  bool HasInvalidTrafficAssignmentState() {
    if (HasInvalidLinkState()) return true;
    const double total_travel_time = network.TotalTravelTime();
    return !IsFiniteNumber(total_travel_time) || total_travel_time < 0.0;
  }

  void PrintInvalidStateWarning(const std::string& context,
                                const std::string& stage) {
    if (!verbose) return;
    if (counters.invalid_value_log_count >= kMaxInvalidStateWarnings) return;
    ++counters.invalid_value_log_count;
    std::cout << "\n[WARN] Invalid TA state detected"
              << " context=" << context
              << " stage=" << stage
              << " (NaN/inf/negative values)."
              << std::endl;
  }

  /**
   * @brief Solves the lower-level TAP (Beckmann user equilibrium) with crash recovery.
   *
   * On first failure (NaN/inf flows), resets the solver state and retries once.
   * Updates ta_compute_seconds with the wall-clock time spent in the TA solver.
   *
   * @param ta_compute_seconds Output: cumulative TA solver time for this call.
   * @param context Label for warning messages (e.g., "standard_eval", "optcond_iter").
   * @param force_reset_before_first_attempt If true, resets solver state before the first attempt.
   * @return true if valid equilibrium flows were obtained.
   */
  bool ComputeTrafficFlowsSafely(double& ta_compute_seconds,
                                 const std::string& context,
                                 bool force_reset_before_first_attempt = false) {
    ta_compute_seconds = 0.0;

    auto run_attempt = [&](bool use_reset) -> bool {
      if (use_reset) {
        approach->Reset();
      }
      const auto ta_start_time = std::chrono::steady_clock::now();
      approach->ComputeTrafficFlows();
      ta_compute_seconds +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - ta_start_time).count();
      return !HasInvalidTrafficAssignmentState();
    };

    if (run_attempt(force_reset_before_first_attempt)) return true;

    ++counters.invalid_ta_state_count;
    PrintInvalidStateWarning(context, "first_attempt");

    if (run_attempt(true)) {
      ++counters.ta_recovery_success_count;
      return true;
    }

    ++counters.ta_recovery_failure_count;
    PrintInvalidStateWarning(context, "after_reset_retry");
    return false;
  }

  /**
   * @brief Projects capacity vector onto the budget constraint set.
   *
   * If budget B(x) > B_max, proportionally reduces all capacity increments (x_a - lb_a)
   * to bring total investment within budget. Uses two passes: proportional scaling,
   * then greedy sweep for any residual excess.
   */
  void EnforceBudgetFeasibility(std::vector<double>& x,
                                const std::vector<double>& lower_bounds) {
    const int n = static_cast<int>(x.size());
    double budget = BudgetFunction(n, x.data());
    if (budget <= budget_upper_bound) return;

    const double multiplier = static_cast<double>(budget_function_multiplier);
    if (multiplier <= 0.0) return;

    double excess = (budget - budget_upper_bound) / multiplier;
    if (excess <= 0.0) return;

    double total_slack = 0.0;
    for (int i = 0; i < n; ++i) {
      total_slack += std::max(0.0, x[i] - lower_bounds[i]);
    }
    if (total_slack <= 0.0) return;

    const double keep_ratio = std::max(0.0, (total_slack - excess) / total_slack);
    for (int i = 0; i < n; ++i) {
      const double slack = std::max(0.0, x[i] - lower_bounds[i]);
      x[i] = lower_bounds[i] + slack * keep_ratio;
    }

    budget = BudgetFunction(n, x.data());
    if (budget <= budget_upper_bound) return;

    excess = (budget - budget_upper_bound) / multiplier;
    for (int i = 0; i < n && excess > 0.0; ++i) {
      const double reducible = std::max(0.0, x[i] - lower_bounds[i]);
      const double reduce_by = std::min(reducible, excess);
      x[i] -= reduce_by;
      excess -= reduce_by;
    }
  }

  /// @brief Updates shared runtime counters with TA timing and budget violation metrics.
  void UpdateRuntimeAndBudgetMetrics(double ta_compute_seconds, double budget) {
    ++counters.total_ta_run_count;
    counters.ta_compute_time_sum_seconds += ta_compute_seconds;
    counters.ta_compute_time_max_seconds = std::max(counters.ta_compute_time_max_seconds, ta_compute_seconds);
    counters.max_budget_violation = std::max(
      counters.max_budget_violation,
      std::max(0.0, budget - budget_upper_bound)
    );
  }

  /// @brief Records an objective evaluation to the convergence trace (delegates to CndStatisticsRecorder).
  void LogQualityTimePoint(const std::string& phase,
                           int step,
                           double objective,
                           double total_travel_time,
                           double budget,
                           double ta_compute_seconds,
                           bool force_relative_gap = false) {
    if (!statistics_recorder.IsRecording()) return;

    double relative_gap = std::numeric_limits<double>::quiet_NaN();
    const bool values_finite =
      IsFiniteNumber(objective) &&
      IsFiniteNumber(total_travel_time) &&
      IsFiniteNumber(budget);
    if (values_finite && statistics_recorder.ShouldSampleRelativeGap(step, force_relative_gap)) {
      relative_gap = static_cast<double>(network.RelativeGap());
    }

    statistics_recorder.LogPoint(
      phase, step, objective, total_travel_time, budget,
      budget_upper_bound, static_cast<double>(budget_threshold),
      ta_compute_seconds, relative_gap
    );
  }

  // --- Progress bar helpers ---

  static void PrintProgressBar(ProgressState& state, const std::string& suffix) {
    if (!state.active) return;

    state.Advance();
    int total = std::max(1, state.total);
    int current = std::min(state.count, total);
    int left = std::max(0, total - current);

    double ratio = static_cast<double>(current) / static_cast<double>(total);
    int filled = static_cast<int>(ratio * kProgressBarWidth);

    std::string bar;
    bar.reserve(kProgressBarWidth);
    for (int i = 0; i < kProgressBarWidth; ++i) {
      if (i < filled) {
        bar.push_back('=');
      } else if (i == filled && filled < kProgressBarWidth) {
        bar.push_back('>');
      } else {
        bar.push_back(' ');
      }
    }

    auto now = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(now - state.start_time).count();
    double rate = elapsed > 0.0 ? static_cast<double>(current) / elapsed : 0.0;

    std::ostringstream line;
    line << '\r'
         << std::setw(3) << static_cast<int>(ratio * 100.0) << "%|"
         << bar << "| "
         << current << "/" << total
         << " left:" << left
         << " [" << std::fixed << std::setprecision(2) << rate << suffix << "]";
    std::cout << line.str() << std::flush;
  }

  static void PrintProgressBarWithMetrics(ProgressState& state,
                                          const std::string& rate_unit,
                                          double objective_function,
                                          double total_travel_time,
                                          double budget_function) {
    if (!state.active) return;

    state.Advance();
    int total = std::max(1, state.total);
    int current = std::min(state.count, total);
    int left = std::max(0, total - current);

    double ratio = static_cast<double>(current) / static_cast<double>(total);
    int filled = static_cast<int>(ratio * kProgressBarWidth);

    std::string bar;
    bar.reserve(kProgressBarWidth);
    for (int i = 0; i < kProgressBarWidth; ++i) {
      if (i < filled) {
        bar.push_back('=');
      } else if (i == filled && filled < kProgressBarWidth) {
        bar.push_back('>');
      } else {
        bar.push_back(' ');
      }
    }

    auto now = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(now - state.start_time).count();
    double rate = elapsed > 0.0 ? static_cast<double>(current) / elapsed : 0.0;

    std::ostringstream line;
    line << '\r'
         << std::setw(3) << static_cast<int>(ratio * 100.0) << "%|"
         << bar << "| "
         << current << "/" << total
         << " left:" << left
         << " [" << std::fixed << std::setprecision(2) << rate << rate_unit << "] "
         << std::defaultfloat
         << "objective_function=" << objective_function
         << " total_travel_time=" << total_travel_time
         << " budget_function=" << budget_function;
    std::cout << line.str() << std::flush;
  }

  // --- NLopt algorithm metadata ---

  /// @brief Returns a human-readable name for the given NLopt algorithm enum.
  static std::string GetAlgorithmName(nlopt::algorithm algo) {
    switch (algo) {
      case nlopt::LD_SLSQP:
        return "LD_SLSQP (Sequential Least-Squares Programming)";
      case nlopt::LD_MMA:
        return "LD_MMA (Method of Moving Asymptotes)";
      case nlopt::LN_COBYLA:
        return "LN_COBYLA (Constrained Optimization BY Linear Approximations)";
      case nlopt::GN_ISRES:
        return "GN_ISRES (Improved Stochastic Ranking Evolution Strategy)";
      case nlopt::GN_AGS:
        return "GN_AGS (Adaptive Global Search)";
      default:
        return "Custom Algorithm (" + std::to_string(algo) + ")";
    }
  }

  /// @brief Returns true if the NLopt algorithm supports nonlinear inequality constraints.
  ///        Used to choose between hard budget constraint and soft penalty fallback.
  static bool SupportsHardBudgetConstraint(nlopt::algorithm algorithm) {
    switch (algorithm) {
      case nlopt::LN_COBYLA:
      case nlopt::GN_ISRES:
      case nlopt::AUGLAG:
      case nlopt::AUGLAG_EQ:
      case nlopt::LD_MMA:
      case nlopt::LD_SLSQP:
        return true;
      default:
        return false;
    }
  }
};

}  // namespace TrafficAssignment

#endif  // CND_OPTIMIZATION_CONTEXT_H

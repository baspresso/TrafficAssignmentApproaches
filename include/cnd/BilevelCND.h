#ifndef BILEVEL_CND_H
#define BILEVEL_CND_H

#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <sstream>
#include <limits>
#include <fstream>
#include <filesystem>
#include "../tap/algorithms/TrafficAssignmentApproach.h"
#include "../tap/core/Network.h"
#include "DirectedConstraintLoader.h"
#include "CndStatisticsRecorder.h"
#include "OptimizationStep.h"
#include "CndOptimizationContext.h"
#include "OptimizationPipeline.h"
#include "OptimalitySensitivity.h"
#include <boost/multiprecision/cpp_dec_float.hpp>

namespace TrafficAssignment {

/**
 * @brief Top-level bilevel Continuous Network Design Problem (CNDP) solver.
 *
 * Solves the bilevel optimization problem from Chiou (2005) Eq. 1-2:
 *
 * Upper level: min F(y) = TSTT(x(y)) + theta * sum(y_a - lb_a)
 *              s.t.  lb_a <= y_a <= ub_a   (capacity bounds)
 *                    theta * sum(y_a - lb_a) <= B   (budget constraint)
 *
 * Lower level: x(y) = argmin Beckmann UE objective (user equilibrium traffic assignment)
 *
 * The solver uses a configurable pipeline of optimization steps (NLopt, optimality condition,
 * OptimLib) executed sequentially. Each step operates on the shared CndOptimizationContext
 * and modifies network link capacities y. After all steps complete, final results are
 * computed and recorded to trace CSV, metadata JSON, and summary CSV.
 *
 * @tparam T Numeric type for flow/capacity computations.
 */
template <typename T>
class BilevelCND {
  using MatrixXd = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
  using MatrixHightPrecision = Eigen::Matrix<
    boost::multiprecision::cpp_dec_float_50, Eigen::Dynamic, Eigen::Dynamic>;

public:
  /**
   * @brief Constructs the bilevel CNDP solver with a pipeline of optimization steps.
   *
   * @param network Transportation network (links, nodes, OD pairs).
   * @param approach Lower-level TAP solver (RouteBased or Tapas).
   * @param constraints Per-link capacity bounds [lb_a, ub_a] (one per link).
   * @param step_configs Pipeline step configurations (from [step.N] INI sections).
   * @param link_capacity_selection_threshold Minimum capacity margin for optimality iterations.
   * @param budget_threshold Tolerance for budget-at-bound detection.
   * @param budget_function_multiplier theta: investment cost multiplier. Chiou (2005) Eq. 2.
   * @param budget_upper_bound B: maximum allowable investment budget.
   * @param metrics_config Output configuration for traces, metadata, and summary.
   * @param route_search_thread_count Parallel threads for route enumeration.
   */
  BilevelCND(
    Network<T>& network,
    std::shared_ptr<TrafficAssignmentApproach<T>> approach,
    const std::vector<DirectedLinkCapacityConstraint>& constraints,
    const std::vector<OptimizationStepConfig>& step_configs,
    T link_capacity_selection_threshold = 1e-3,
    T budget_threshold = 1e-1,
    T budget_function_multiplier = 5,
    double budget_upper_bound = 100000,
    const CndMetricsConfig& metrics_config = CndMetricsConfig(),
    std::size_t route_search_thread_count = 1,
    const std::string& progress_format = "bar"
  )
    : network_(network),
      approach_(approach),
      constraints_(constraints),
      step_configs_(step_configs),
      link_capacity_selection_threshold_(link_capacity_selection_threshold),
      budget_threshold_(budget_threshold),
      budget_function_multiplier_(budget_function_multiplier),
      budget_upper_bound_(budget_upper_bound),
      verbose_(true),
      statistics_recorder_(network),
      metrics_config_(metrics_config),
      route_search_thread_count_(route_search_thread_count),
      progress_format_(progress_format)
  {
    ValidateConstraints(constraints, network);
    if (step_configs.empty()) {
      throw std::invalid_argument("step_configs must not be empty");
    }
    if (budget_upper_bound <= 0) {
      throw std::invalid_argument("budget_upper_bound must be positive");
    }
  }

  ~BilevelCND() = default;

  void SetMetricsConfig(const CndMetricsConfig& metrics_config) {
    metrics_config_ = metrics_config;
  }

  const CndMetricsConfig& GetMetricsConfig() const {
    return metrics_config_;
  }

  void SetRouteSearchThreadCount(std::size_t thread_count) {
    route_search_thread_count_ = thread_count;
    ApplyApproachRuntimeOptions();
  }

  std::size_t GetRouteSearchThreadCount() const {
    return route_search_thread_count_;
  }

  void SetVerbose(bool verbose) {
    verbose_ = verbose;
  }

  /**
   * @brief Main entry point: builds pipeline, creates context, executes, and records results.
   *
   * Initializes capacities to lower bounds, runs the optimization pipeline, computes
   * final TA equilibrium, and outputs optimality condition statistics and solution CSV.
   * On failure, records partial results before re-throwing.
   */
  void ComputeNetworkDesign() {
    ApplyApproachRuntimeOptions();
    PrintConfiguration();
    const auto optimization_start_time = std::chrono::steady_clock::now();

    typename CndOptimizationContext<T>::RuntimeCounters counters;
    counters.Reset();
    StartStatisticsRecording(counters);

    CndOptimizationContext<T> ctx(
      network_, approach_, constraints_,
      budget_function_multiplier_, budget_upper_bound_,
      budget_threshold_, link_capacity_selection_threshold_,
      route_search_thread_count_, verbose_,
      counters, statistics_recorder_,
      progress_format_
    );

    try {
      InitializeCapacities();

      auto pipeline = OptimizationPipeline<T>::BuildFromConfigs(step_configs_);
      pipeline.Execute(ctx);

      ComputeFinalResults(optimization_start_time, counters);
    } catch (...) {
      const double failed_elapsed_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - optimization_start_time)
          .count();
      StopStatisticsRecording(
        counters,
        failed_elapsed_seconds,
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        false
      );
      throw;
    }
  }

private:
  // --- Core dependencies ---
  Network<T>& network_;                        ///< Transportation network with mutable link capacities.
  std::shared_ptr<TrafficAssignmentApproach<T>> approach_; ///< Lower-level TAP solver.
  std::vector<DirectedLinkCapacityConstraint> constraints_; ///< Per-link capacity bounds [lb_a, ub_a].

  // --- Pipeline config ---
  std::vector<OptimizationStepConfig> step_configs_; ///< Ordered step configurations from [step.N] sections.

  // --- Shared parameters ---
  T link_capacity_selection_threshold_;        ///< Min capacity margin for optimality condition link selection.
  T budget_threshold_;                         ///< Tolerance for budget-at-bound detection.
  T budget_function_multiplier_;               ///< theta: investment cost multiplier. Chiou (2005) Eq. 2.
  double budget_upper_bound_;                  ///< B: maximum allowable investment budget.
  bool verbose_;

  // --- Statistics ---
  CndStatisticsRecorder<T> statistics_recorder_;
  CndMetricsConfig metrics_config_;
  std::size_t route_search_thread_count_;
  std::string progress_format_;

  // --- Constants for statistics progress ---
  static constexpr int kProgressBarWidth = 30;

  static void ValidateConstraints(
      const std::vector<DirectedLinkCapacityConstraint>& constraints,
      const Network<T>& network) {
    if (constraints.empty()) {
      throw std::invalid_argument("constraints must not be empty");
    }
    if (static_cast<int>(constraints.size()) != network.number_of_links()) {
      throw std::invalid_argument(
        "constraints size (" + std::to_string(constraints.size()) +
        ") must match number of links (" + std::to_string(network.number_of_links()) + ")");
    }
  }

  void ApplyApproachRuntimeOptions() {
    if (!approach_) return;
    approach_->SetRouteSearchThreadCount(route_search_thread_count_);
    route_search_thread_count_ = approach_->GetRouteSearchThreadCount();
  }

  void PrintConfiguration() {
    if (!verbose_) return;
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "BILEVEL CONTINUOUS NETWORK DESIGN PROBLEM" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    std::cout << "\nProblem Configuration:" << std::endl;
    std::cout << "  Design variables: " << constraints_.size() << std::endl;
    std::cout << "  Pipeline steps: " << step_configs_.size() << std::endl;
    for (std::size_t i = 0; i < step_configs_.size(); ++i) {
      const auto& sc = step_configs_[i];
      std::cout << "    [" << i << "] type=" << sc.type;
      if (!sc.algorithm.empty()) std::cout << " algorithm=" << sc.algorithm;
      if (!sc.local_algorithm.empty()) std::cout << " local=" << sc.local_algorithm;
      std::cout << " max_iter=" << sc.max_iterations;
      std::cout << " tol=" << sc.tolerance;
      std::cout << std::endl;
    }
    std::cout << "  Route search threads: " << route_search_thread_count_ << std::endl;
  }

  /// @brief Sets all link capacities to their lower bounds (starting point for optimization).
  void InitializeCapacities() {
    for (std::size_t i = 0; i < constraints_.size(); ++i) {
      network_.mutable_links()[i].capacity = constraints_[i].lower_bound;
    }
  }

  /// @brief Runs final TA solve, logs results, writes summary and solution CSV.
  void ComputeFinalResults(
      std::chrono::steady_clock::time_point optimization_start_time,
      typename CndOptimizationContext<T>::RuntimeCounters& counters) {
    // Final TA computation
    CndOptimizationContext<T> ctx(
      network_, approach_, constraints_,
      budget_function_multiplier_, budget_upper_bound_,
      budget_threshold_, link_capacity_selection_threshold_,
      route_search_thread_count_, verbose_,
      counters, statistics_recorder_,
      progress_format_
    );

    double final_ta_compute_seconds = 0.0;
    if (!ctx.ComputeTrafficFlowsSafely(final_ta_compute_seconds, "final")) {
      throw std::runtime_error("TA computation failed on final recovery attempt.");
    }

    const auto optimization_end_time = std::chrono::steady_clock::now();
    const double optimization_elapsed_seconds =
      std::chrono::duration<double>(optimization_end_time - optimization_start_time).count();
    const double total_travel_time = network_.TotalTravelTime();
    const double budget_function = static_cast<double>(ctx.Budget());
    const double objective_function = total_travel_time + budget_function;
    if (!CndOptimizationContext<T>::IsFiniteNumber(total_travel_time) ||
        !CndOptimizationContext<T>::IsFiniteNumber(budget_function) ||
        !CndOptimizationContext<T>::IsFiniteNumber(objective_function)) {
      ctx.PrintInvalidStateWarning("final", "objective_invalid");
      throw std::runtime_error("final objective contains non-finite values.");
    }
    ctx.UpdateRuntimeAndBudgetMetrics(final_ta_compute_seconds, budget_function);
    ctx.LogQualityTimePoint("final",
                            counters.optcond_trace_step,
                            objective_function,
                            total_travel_time,
                            budget_function,
                            final_ta_compute_seconds,
                            true);
    // Always emit structured [RESULT] line for machine parsing
    std::cout << "[RESULT]"
              << " optimization_time=" << std::fixed << std::setprecision(2) << optimization_elapsed_seconds
              << " objective_function=" << std::fixed << std::setprecision(10) << objective_function
              << " total_travel_time=" << std::fixed << std::setprecision(10) << total_travel_time
              << " budget_function=" << std::fixed << std::setprecision(10) << budget_function
              << std::endl;

    if (verbose_) {
      std::cout << "optimization_time = "
                << std::fixed << std::setprecision(2)
                << optimization_elapsed_seconds << "s "
                << "objective_function=" << std::fixed << std::setprecision(10) << objective_function
                << " total_travel_time=" << std::fixed << std::setprecision(10) << total_travel_time
                << " budget_function=" << std::fixed << std::setprecision(10) << budget_function << std::endl;
    }

    StopStatisticsRecording(
      counters,
      optimization_elapsed_seconds,
      objective_function,
      total_travel_time,
      budget_function,
      true
    );

    statistics_recorder_.WriteSolutionCSV(constraints_);
    RecordStatistics(ctx);
  }

  // =====================================================================
  // Statistics recording
  // =====================================================================

  void StartStatisticsRecording(
      typename CndOptimizationContext<T>::RuntimeCounters& counters) {
    CndRunMetadata metadata;
    metadata.dataset_name = network_.name();
    metadata.approach_name = approach_ ? approach_->GetApproachName() : "UnknownApproach";
    // Build algorithm name from first nlopt step, or "N/A"
    metadata.nlopt_algorithm_name = "N/A";
    int total_standard = 0;
    int total_optcond = 0;
    for (const auto& sc : step_configs_) {
      if (sc.type == "nlopt") {
        if (metadata.nlopt_algorithm_name == "N/A") {
          metadata.nlopt_algorithm_name = sc.algorithm;
        }
        total_standard += sc.max_iterations;
      } else if (sc.type == "optimality_condition") {
        total_optcond += sc.max_iterations;
      }
    }
    metadata.scenario_name = metrics_config_.scenario_name;
    metadata.route_search_thread_count = route_search_thread_count_;
    metadata.max_standard_iterations = total_standard;
    metadata.max_optimality_condition_iterations = total_optcond;
    metadata.optimization_tolerance = step_configs_.empty() ? 0.0 : step_configs_[0].tolerance;
    metadata.budget_upper_bound = budget_upper_bound_;
    metadata.budget_threshold = static_cast<double>(budget_threshold_);
    metadata.budget_function_multiplier = static_cast<double>(budget_function_multiplier_);
    metadata.link_capacity_selection_threshold = static_cast<double>(link_capacity_selection_threshold_);
    metadata.design_variables = constraints_.size();
    metadata.number_of_links = static_cast<std::size_t>(network_.number_of_links());
    metadata.number_of_od_pairs = static_cast<std::size_t>(network_.number_of_od_pairs());
    statistics_recorder_.StartRun(metrics_config_, metadata);
  }

  void StopStatisticsRecording(
      const typename CndOptimizationContext<T>::RuntimeCounters& counters,
      double total_elapsed_seconds,
      double final_objective,
      double final_total_travel_time,
      double final_budget,
      bool success) {
    if (!statistics_recorder_.IsRecording()) return;

    CndRunSummary summary;
    summary.status = success ? "success" : "failure";
    summary.total_elapsed_seconds = total_elapsed_seconds;
    summary.standard_eval_count = counters.standard_trace_step;
    summary.optimality_condition_iteration_count = counters.optcond_trace_step;
    summary.total_ta_run_count = counters.total_ta_run_count;
    summary.invalid_ta_state_count = counters.invalid_ta_state_count;
    summary.ta_recovery_success_count = counters.ta_recovery_success_count;
    summary.ta_recovery_failure_count = counters.ta_recovery_failure_count;
    summary.avg_ta_compute_seconds =
      counters.total_ta_run_count > 0
        ? counters.ta_compute_time_sum_seconds / counters.total_ta_run_count
        : 0.0;
    summary.max_ta_compute_seconds = counters.ta_compute_time_max_seconds;
    summary.final_objective = final_objective;
    summary.final_total_travel_time = final_total_travel_time;
    summary.final_budget = final_budget;
    if (std::isfinite(final_budget)) {
      summary.final_budget_violation = std::max(0.0, final_budget - budget_upper_bound_);
    } else {
      summary.final_budget_violation = std::numeric_limits<double>::quiet_NaN();
    }
    summary.best_feasible_objective = statistics_recorder_.best_feasible_objective();
    summary.max_budget_violation = counters.max_budget_violation;
    statistics_recorder_.StopRun(summary);
  }

  // =====================================================================
  // Statistics CSV output
  // =====================================================================

  void RecordStatistics(CndOptimizationContext<T>& ctx) {
    std::vector<T> link_condition_result = OptimalityCondition(ctx);
    std::vector<T> flow(network_.number_of_links(), 0);
    std::vector<T> lower_bound(network_.number_of_links(), 0);
    std::vector<T> upper_bound(network_.number_of_links(), 0);
    std::vector<T> capacity(network_.number_of_links(), 0);

    for (int i = 0; i < network_.number_of_links(); i++) {
      flow[i] = network_.mutable_links()[i].flow;
      lower_bound[i] = constraints_[i].lower_bound;
      upper_bound[i] = constraints_[i].upper_bound;
      capacity[i] = network_.mutable_links()[i].capacity;
    }
    std::vector<std::vector<T>> vectors = {
      flow, lower_bound, upper_bound, capacity, link_condition_result
    };
    std::vector<std::string> headers = {
      "flow", "lower_bound", "upper_bound", "capacity", "condition_result"
    };
    writeVectorsToCSV(vectors, headers);
  }

  /**
   * @brief Computes the optimality function phi_a for every link (post-optimization diagnostic).
   *
   * At a first-order optimal solution, phi_a should be equal across all active design links.
   * The spread of phi values indicates solution quality. Chiou (2005) Eq. 9.
   *
   * @return Vector of phi_a values, one per link.
   */
  std::vector<T> OptimalityCondition(CndOptimizationContext<T>& ctx) {
    const int n_links = network_.number_of_links();
    std::vector<T> link_condition_result(n_links, T(0));
    double ta_compute_seconds = 0.0;
    ctx.ComputeTrafficFlowsSafely(ta_compute_seconds, "optimality_condition");

    // Build OD cache inline (same logic as OptimalityConditionStep)
    const int n_od = network_.number_of_od_pairs();
    struct OdCache {
      std::vector<std::vector<int>> routes;
      MatrixHightPrecision          jacobi_inverse;
      MatrixHightPrecision          e_col;
      boost::multiprecision::cpp_dec_float_50 denominator;
      T                             demand;
      bool                          valid = false;
    };

    std::vector<OdCache> od_cache(n_od);
    for (int od = 0; od < n_od; ++od) {
      const int rc = network_.mutable_od_pairs()[od].GetRoutesCount();
      if (rc <= 0) continue;
      OdCache& c = od_cache[od];
      c.routes  = network_.mutable_od_pairs()[od].GetRoutes();
      c.demand  = network_.od_pairs()[od].GetDemand();

      // Filter routes where ALL links have trivial BPR params (zero Jacobian rows)
      std::erase_if(c.routes, [&](const std::vector<int>& route) {
        return OptimalitySensitivity<T>::IsRouteTrivial(route, network_.mutable_links());
      });
      const int filtered_rc = static_cast<int>(c.routes.size());
      if (filtered_rc <= 0) continue;

      c.e_col   = MatrixHightPrecision::Constant(filtered_rc, 1,
        boost::multiprecision::cpp_dec_float_50(1));
      auto jacobi_hp = ConvertEigenMatrix(
          RoutesJacobiMatrix(c.routes, network_.mutable_links()));
      c.jacobi_inverse = jacobi_hp.inverse();
      c.denominator    = (c.e_col.transpose() * c.jacobi_inverse * c.e_col)(0, 0);
      if (boost::multiprecision::abs(c.denominator) <
          static_cast<boost::multiprecision::cpp_dec_float_50>(
            CndOptimizationContext<T>::kDenominatorZeroGuard)) {
        continue;
      }
      c.valid = true;
    }

    // Progress bar for statistics
    typename CndOptimizationContext<T>::ProgressState stats_progress;
    stats_progress.Start(n_links, "StatisticsRecording", ctx.progress_format);

    for (int link_index = 0; link_index < n_links; ++link_index) {
      PrintStatisticsProgressAt(stats_progress, link_index);
      // OptimalityFunctionFromCache inline
      T result = T(0);
      for (const OdCache& c : od_cache) {
        if (!c.valid) continue;
        auto cap_der = ConvertEigenMatrix(CapacityDerColumn(c.routes, link_index));
        const auto numerator =
            (c.e_col.transpose() * c.jacobi_inverse * cap_der)(0, 0);
        const T local_value = static_cast<T>(-numerator / c.denominator);
        if (!ctx.IsFiniteScalar(local_value)) {
          result = ctx.NaNValue();
          break;
        }
        result += c.demand * local_value;
      }
      if (ctx.IsFiniteScalar(result)) {
        result /= budget_function_multiplier_;
        if (!ctx.IsFiniteScalar(result)) result = ctx.NaNValue();
      }
      link_condition_result[link_index] = result;
    }
    PrintStatisticsProgressAt(stats_progress, n_links);
    stats_progress.Finish();
    return link_condition_result;
  }

  // --- Helper methods for statistics ---

  MatrixXd CapacityDerColumn(const std::vector<std::vector<int>>& routes, int link_index) {
    MatrixXd res = MatrixXd::Constant(routes.size(), 1, 0.0);
    for (std::size_t i = 0; i < routes.size(); i++) {
      if (std::find(routes[i].begin(), routes[i].end(), link_index) != routes[i].end()) {
        res(i, 0) = network_.links()[link_index].DelayCapacityDer();
      }
    }
    return res;
  }

  MatrixHightPrecision ConvertEigenMatrix(const MatrixXd& source) {
    MatrixHightPrecision target(source.rows(), source.cols());
    for (int i = 0; i < source.rows(); ++i) {
      for (int j = 0; j < source.cols(); ++j) {
        target(i, j) = static_cast<boost::multiprecision::cpp_dec_float_50>(source(i, j));
      }
    }
    return target;
  }

  void writeVectorsToCSV(const std::vector<std::vector<T>>& vectors,
                         const std::vector<std::string>& headers) {
    auto output_dir = statistics_recorder_.output_dir();
    if (output_dir.empty()) {
      output_dir = std::filesystem::current_path() / "performance_results" / network_.name();
    }
    std::filesystem::create_directories(output_dir);
    auto file_path = output_dir / (network_.name() + "_optimality.csv");
    std::ofstream file(file_path);
    for (size_t i = 0; i < headers.size(); ++i) {
      file << headers[i];
      if (i != headers.size() - 1) file << ",";
    }
    file << "\n";

    size_t rows = vectors[0].size();
    for (size_t row = 0; row < rows; ++row) {
      for (size_t col = 0; col < vectors.size(); ++col) {
        file << vectors[col][row];
        if (col != vectors.size() - 1) file << ",";
      }
      file << "\n";
    }

    file.close();
  }

  void PrintStatisticsProgressAt(
      typename CndOptimizationContext<T>::ProgressState& stats_progress,
      int current) {
    if (stats_progress.format == "none") return;

    int total = std::max(1, stats_progress.total);

    if (stats_progress.format == "line") {
      auto now = std::chrono::steady_clock::now();
      double elapsed = std::chrono::duration<double>(now - stats_progress.start_time).count();
      double links_per_sec = elapsed > 0.0 ? static_cast<double>(current) / elapsed : 0.0;
      std::cout << "[PROGRESS] step=" << stats_progress.step_name
                << " current=" << current
                << " total=" << total
                << " rate=" << std::fixed << std::setprecision(2) << links_per_sec
                << " unit=links/s\n" << std::flush;
      return;
    }

    // format == "bar": traditional \r-based progress bar
    int left  = std::max(0, total - current);

    double ratio  = static_cast<double>(current) / static_cast<double>(total);
    int    filled = static_cast<int>(ratio * kProgressBarWidth);

    std::string bar;
    bar.reserve(kProgressBarWidth);
    for (int i = 0; i < kProgressBarWidth; ++i) {
      if (i < filled)                                    bar.push_back('=');
      else if (i == filled && filled < kProgressBarWidth) bar.push_back('>');
      else                                               bar.push_back(' ');
    }

    auto   now          = std::chrono::steady_clock::now();
    double elapsed      = std::chrono::duration<double>(now - stats_progress.start_time).count();
    double links_per_sec = elapsed > 0.0 ? static_cast<double>(current) / elapsed : 0.0;

    std::ostringstream line;
    line << '\r'
         << std::setw(3) << static_cast<int>(ratio * 100.0) << "%|"
         << bar << "| "
         << current << "/" << total
         << " left:" << left
         << " [" << std::fixed << std::setprecision(2) << links_per_sec << "links/s]";
    std::cout << line.str() << std::flush;
  }

  BilevelCND(const BilevelCND&) = delete;
  BilevelCND& operator=(const BilevelCND&) = delete;
};

}  // namespace TrafficAssignment

#endif  // BILEVEL_CND_H

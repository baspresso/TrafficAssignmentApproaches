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
#include "../tap/algorithms/common/TrafficAssignmentApproach.h"
#include "../tap/core/Network.h"
#include "DirectedConstraintLoader.h"
#include "CndStatisticsRecorder.h"
#include "OptimizationStep.h"
#include "CndOptimizationContext.h"
#include "OptimizationPipeline.h"
#include <boost/multiprecision/cpp_dec_float.hpp>

namespace TrafficAssignment {

template <typename T>
class BilevelCND {
  using MatrixXd = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
  using MatrixHightPrecision = Eigen::Matrix<
    boost::multiprecision::cpp_dec_float_50, Eigen::Dynamic, Eigen::Dynamic>;

public:
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
    std::size_t route_search_thread_count = 1
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
      route_search_thread_count_(route_search_thread_count)
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
      counters, statistics_recorder_
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
  Network<T>& network_;
  std::shared_ptr<TrafficAssignmentApproach<T>> approach_;
  std::vector<DirectedLinkCapacityConstraint> constraints_;

  // --- Pipeline config ---
  std::vector<OptimizationStepConfig> step_configs_;

  // --- Shared parameters ---
  T link_capacity_selection_threshold_;
  T budget_threshold_;
  T budget_function_multiplier_;
  double budget_upper_bound_;
  bool verbose_;

  // --- Statistics ---
  CndStatisticsRecorder<T> statistics_recorder_;
  CndMetricsConfig metrics_config_;
  std::size_t route_search_thread_count_;

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

  void InitializeCapacities() {
    for (std::size_t i = 0; i < constraints_.size(); ++i) {
      network_.mutable_links()[i].capacity = constraints_[i].lower_bound;
    }
  }

  void ComputeFinalResults(
      std::chrono::steady_clock::time_point optimization_start_time,
      typename CndOptimizationContext<T>::RuntimeCounters& counters) {
    // Final TA computation
    CndOptimizationContext<T> ctx(
      network_, approach_, constraints_,
      budget_function_multiplier_, budget_upper_bound_,
      budget_threshold_, link_capacity_selection_threshold_,
      route_search_thread_count_, verbose_,
      counters, statistics_recorder_
    );

    double final_ta_compute_seconds = 0.0;
    if (!ctx.ComputeTrafficFlowsSafely(final_ta_compute_seconds, "final")) {
      throw std::runtime_error("TA computation failed on final recovery attempt.");
    }

    std::cout << ctx.Budget() << '\n';
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
    std::ostringstream timing_line;
    timing_line << "optimization_time = "
                << std::fixed << std::setprecision(2)
                << optimization_elapsed_seconds << "s " << std::setprecision(10)
                << "objective_function=" << objective_function
                << " total_travel_time=" << total_travel_time
                << " budget_function=" << budget_function << "\n";
    std::cout << timing_line.str();

    StopStatisticsRecording(
      counters,
      optimization_elapsed_seconds,
      objective_function,
      total_travel_time,
      budget_function,
      true
    );

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
      c.e_col   = MatrixHightPrecision::Constant(rc, 1,
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
    stats_progress.total = n_links;
    stats_progress.start_time = std::chrono::steady_clock::now();
    std::cout << "StatisticsRecording-start\n";

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
    std::cout << '\n';
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
    std::filesystem::path current_path = "C:/Projects/TrafficAssignmentApproaches";
    auto file_path = current_path / "data" / "TransportationNetworks" / network_.name() /
                     (network_.name() + "_optimality.csv");
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
    int total = std::max(1, stats_progress.total);
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

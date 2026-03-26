#ifndef CND_STATISTICS_RECORDER_H
#define CND_STATISTICS_RECORDER_H

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include "../tap/core/Network.h"

namespace TrafficAssignment {

struct DirectedLinkCapacityConstraint;  // forward declaration (defined in DirectedConstraintLoader.h)

/**
 * @brief Configuration for CNDP metrics collection and output.
 */
struct CndMetricsConfig {
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

/**
 * @brief Metadata describing a CNDP run (written to JSON and summary CSV).
 */
struct CndRunMetadata {
  std::string dataset_name;
  std::string approach_name;
  std::string nlopt_algorithm_name;
  std::string scenario_name;
  std::size_t route_search_thread_count = 1;
  int max_standard_iterations = 0;
  int max_optimality_condition_iterations = 0;
  double optimization_tolerance = 0.0;
  double budget_upper_bound = 0.0;
  double budget_threshold = 0.0;
  double budget_function_multiplier = 0.0;
  double link_capacity_selection_threshold = 0.0;
  std::size_t design_variables = 0;
  std::size_t number_of_links = 0;
  std::size_t number_of_od_pairs = 0;
};

/**
 * @brief Single data point in the convergence trace (one per objective evaluation).
 */
struct CndTracePoint {
  double elapsed_seconds = 0.0;         ///< Wall-clock time since run start.
  std::string phase;                    ///< Optimization phase ("standard_eval", "optcond_iter", etc.).
  int step = 0;                         ///< Step counter within the current phase.
  double objective = 0.0;               ///< Upper-level objective F(y) = TSTT + BudgetCost.
  double best_feasible_objective = std::numeric_limits<double>::quiet_NaN(); ///< Best objective among budget-feasible solutions.
  double total_travel_time = 0.0;       ///< TSTT at this evaluation.
  double budget = 0.0;                  ///< Investment cost at this evaluation.
  double budget_violation = 0.0;        ///< max(0, budget - budget_upper_bound).
  double relative_gap = std::numeric_limits<double>::quiet_NaN(); ///< TAP relative gap (sampled periodically).
  double ta_compute_seconds = 0.0;      ///< Time spent solving the lower-level TAP.
};

/**
 * @brief Aggregated summary of a completed CNDP run (appended to summary CSV).
 */
struct CndRunSummary {
  std::string status = "success";       ///< "success" or "failure".
  double total_elapsed_seconds = 0.0;
  int standard_eval_count = 0;
  int optimality_condition_iteration_count = 0;
  int total_ta_run_count = 0;
  int invalid_ta_state_count = 0;
  int ta_recovery_success_count = 0;
  int ta_recovery_failure_count = 0;
  double avg_ta_compute_seconds = 0.0;
  double max_ta_compute_seconds = 0.0;
  double final_objective = 0.0;
  double final_total_travel_time = 0.0;
  double final_budget = 0.0;
  double final_budget_violation = 0.0;
  double best_feasible_objective = std::numeric_limits<double>::quiet_NaN();
  double max_budget_violation = 0.0;
};

/**
 * @brief Records convergence traces, metadata, and summary statistics for bilevel CNDP experiments.
 *
 * Outputs three file types per run:
 * - Trace CSV: per-iteration objective, budget, RGAP, timing (for convergence plots)
 * - Metadata JSON: algorithm configuration and problem parameters
 * - Summary CSV: one-line append-only record per run (for comparison tables)
 *
 * The summary CSV includes a Scenario column for grouping multiple runs in analysis.
 *
 * @tparam T Numeric type for flow/capacity computations.
 */
template <typename T>
class CndStatisticsRecorder {
public:
  explicit CndStatisticsRecorder(Network<T>& network)
    : network_(network),
      is_recording_(false),
      best_feasible_objective_(std::numeric_limits<double>::infinity()),
      flushed_points_(0) {}

  /// @brief Initializes a new recording session, creating output files.
  void StartRun(const CndMetricsConfig& config, const CndRunMetadata& metadata) {
    config_ = config;
    metadata_ = metadata;
    metadata_.dataset_name = network_.name();
    if (metadata_.number_of_links == 0) {
      metadata_.number_of_links = static_cast<std::size_t>(network_.number_of_links());
    }
    if (metadata_.number_of_od_pairs == 0) {
      metadata_.number_of_od_pairs = static_cast<std::size_t>(network_.number_of_od_pairs());
    }

    run_id_ = config_.run_id.empty() ? GenerateRunId() : config_.run_id;
    output_dir_ = ResolveOutputDir();
    std::filesystem::create_directories(output_dir_);

    const std::string approach_name =
      metadata_.approach_name.empty() ? "UnknownApproach" : metadata_.approach_name;

    trace_file_path_ =
      output_dir_ / ("BilevelCND_" + approach_name + "_" + run_id_ + "_quality_time.csv");
    metadata_file_path_ =
      output_dir_ / ("BilevelCND_" + approach_name + "_" + run_id_ + "_metadata.json");
    summary_file_path_ = output_dir_ / "BilevelCND_run_summary.csv";

    trace_points_.clear();
    flushed_points_ = 0;
    best_feasible_objective_ = std::numeric_limits<double>::infinity();
    run_start_time_ = std::chrono::steady_clock::now();
    is_recording_ = true;

    if (config_.enable_trace) {
      InitializeTraceFile();
    }
    if (config_.write_metadata_json) {
      WriteMetadataFile();
    }
  }

  bool IsRecording() const {
    return is_recording_;
  }

  bool ShouldSampleRelativeGap(int step, bool force_sample) const {
    if (!config_.enable_relative_gap) {
      return false;
    }
    if (force_sample) {
      return true;
    }
    if (config_.relative_gap_sample_period <= 0) {
      return false;
    }
    return step % config_.relative_gap_sample_period == 0;
  }

  /// @brief Records a single trace point (objective evaluation) to the convergence trace.
  void LogPoint(const std::string& phase,
                int step,
                double objective,
                double total_travel_time,
                double budget,
                double budget_upper_bound,
                double budget_feasibility_tolerance,
                double ta_compute_seconds,
                double relative_gap) {
    if (!is_recording_) {
      return;
    }

    const double budget_violation = std::max(0.0, budget - budget_upper_bound);
    if (budget_violation <= budget_feasibility_tolerance) {
      best_feasible_objective_ = std::min(best_feasible_objective_, objective);
    }

    const double best_feasible =
      std::isfinite(best_feasible_objective_)
        ? best_feasible_objective_
        : std::numeric_limits<double>::quiet_NaN();

    if (!config_.enable_trace) {
      return;
    }

    trace_points_.push_back(
      {
        ElapsedSeconds(),
        phase,
        step,
        objective,
        best_feasible,
        total_travel_time,
        budget,
        budget_violation,
        relative_gap,
        ta_compute_seconds
      }
    );

    if (config_.enable_trace && config_.flush_every_n_points > 0) {
      if (trace_points_.size() - flushed_points_ >=
          static_cast<std::size_t>(config_.flush_every_n_points)) {
        AppendPendingTracePoints();
      }
    }
  }

  /// @brief Finalizes the recording session, flushing trace and writing summary CSV.
  void StopRun(const CndRunSummary& summary) {
    if (!is_recording_) {
      return;
    }
    if (config_.enable_trace) {
      AppendPendingTracePoints();
    }
    if (config_.write_summary_csv) {
      WriteSummary(summary);
    }
    is_recording_ = false;
  }

  const std::string& run_id() const {
    return run_id_;
  }

  const std::filesystem::path& output_dir() const {
    return output_dir_;
  }

  /// @brief Writes per-link solution CSV with optimized capacities, flows, and bounds.
  void WriteSolutionCSV(
      const std::vector<DirectedLinkCapacityConstraint>& constraints) {
    if (output_dir_.empty()) return;

    const std::string approach_name =
      metadata_.approach_name.empty() ? "UnknownApproach" : metadata_.approach_name;
    const auto solution_path =
      output_dir_ / ("BilevelCND_" + approach_name + "_" + run_id_ + "_solution.csv");

    std::ofstream file(solution_path, std::ios::out);
    file << "link_id,init_node,term_node,optimized_capacity,"
            "flow,lower_bound,upper_bound\n";
    file << std::setprecision(10);

    const int n_links = network_.number_of_links();
    for (int i = 0; i < n_links; ++i) {
      const auto& link = network_.links()[i];
      file << i << ","
           << link.init << ","
           << link.term << ","
           << link.capacity << ","
           << link.flow << ","
           << constraints[i].lower_bound << ","
           << constraints[i].upper_bound << "\n";
    }
  }

  double best_feasible_objective() const {
    if (std::isfinite(best_feasible_objective_)) {
      return best_feasible_objective_;
    }
    return std::numeric_limits<double>::quiet_NaN();
  }

private:
  Network<T>& network_;
  CndMetricsConfig config_;
  CndRunMetadata metadata_;
  bool is_recording_;
  std::string run_id_;
  std::filesystem::path output_dir_;
  std::filesystem::path trace_file_path_;
  std::filesystem::path metadata_file_path_;
  std::filesystem::path summary_file_path_;
  std::vector<CndTracePoint> trace_points_;
  std::chrono::steady_clock::time_point run_start_time_;
  double best_feasible_objective_;
  std::size_t flushed_points_;

  std::string GenerateRunId() const {
    const auto now = std::chrono::system_clock::now();
    const auto epoch_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count();
    return std::to_string(epoch_ms);
  }

  std::string CurrentTimestampIso8601() const {
    const auto now = std::chrono::system_clock::now();
    const std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tm_now {};
#ifdef _WIN32
    localtime_s(&tm_now, &t);
#else
    localtime_r(&t, &tm_now);
#endif
    std::ostringstream oss;
    oss << std::put_time(&tm_now, "%Y-%m-%d %H:%M:%S");
    return oss.str();
  }

  double ElapsedSeconds() const {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - run_start_time_).count();
  }

  std::filesystem::path ResolveOutputDir() const {
    const auto cwd = std::filesystem::current_path();
    std::filesystem::path root = config_.output_root.empty()
                                   ? std::filesystem::path("performance_results")
                                   : std::filesystem::path(config_.output_root);

    if (!root.is_absolute()) {
      if (std::filesystem::exists(cwd / root)) {
        root = cwd / root;
      } else if (std::filesystem::exists(cwd.parent_path() / root)) {
        root = cwd.parent_path() / root;
      } else {
        root = cwd / root;
      }
    }
    if (config_.append_dataset_subdir) {
      return root / network_.name();
    }
    return root;
  }

  void InitializeTraceFile() const {
    std::ofstream file(trace_file_path_, std::ios::out);
    file << "RunId,Scenario,Dataset,Approach,Algorithm,ElapsedTime(s),Phase,Step,Objective,"
            "BestFeasibleObjective,TotalTravelTime,Budget,BudgetViolation,RelativeGap,"
            "TAComputeTime(s)\n";
  }

  void AppendPendingTracePoints() {
    if (flushed_points_ >= trace_points_.size()) {
      return;
    }
    std::ofstream file(trace_file_path_, std::ios::app);
    file << std::setprecision(10);
    const std::string approach_name =
      metadata_.approach_name.empty() ? "UnknownApproach" : metadata_.approach_name;
    for (std::size_t i = flushed_points_; i < trace_points_.size(); ++i) {
      const auto& point = trace_points_[i];
      file << run_id_ << ","
           << metadata_.scenario_name << ","
           << metadata_.dataset_name << ","
           << approach_name << ","
           << metadata_.nlopt_algorithm_name << ","
           << point.elapsed_seconds << ","
           << point.phase << ","
           << point.step << ","
           << point.objective << ","
           << point.best_feasible_objective << ","
           << point.total_travel_time << ","
           << point.budget << ","
           << point.budget_violation << ","
           << point.relative_gap << ","
           << point.ta_compute_seconds << "\n";
    }
    flushed_points_ = trace_points_.size();
  }

  void WriteMetadataFile() const {
    std::ofstream file(metadata_file_path_, std::ios::out);
    const std::string approach_name =
      metadata_.approach_name.empty() ? "UnknownApproach" : metadata_.approach_name;
    file << "{\n";
    file << "  \"run_id\": \"" << run_id_ << "\",\n";
    file << "  \"timestamp\": \"" << CurrentTimestampIso8601() << "\",\n";
    file << "  \"scenario\": \"" << metadata_.scenario_name << "\",\n";
    file << "  \"dataset\": \"" << metadata_.dataset_name << "\",\n";
    file << "  \"approach\": \"" << approach_name << "\",\n";
    file << "  \"algorithm\": \"" << metadata_.nlopt_algorithm_name << "\",\n";
    file << "  \"route_search_thread_count\": " << metadata_.route_search_thread_count << ",\n";
    file << "  \"design_variables\": " << metadata_.design_variables << ",\n";
    file << "  \"number_of_links\": " << metadata_.number_of_links << ",\n";
    file << "  \"number_of_od_pairs\": " << metadata_.number_of_od_pairs << ",\n";
    file << "  \"max_standard_iterations\": " << metadata_.max_standard_iterations << ",\n";
    file << "  \"max_optimality_condition_iterations\": "
         << metadata_.max_optimality_condition_iterations << ",\n";
    file << "  \"optimization_tolerance\": " << metadata_.optimization_tolerance << ",\n";
    file << "  \"budget_upper_bound\": " << metadata_.budget_upper_bound << ",\n";
    file << "  \"budget_threshold\": " << metadata_.budget_threshold << ",\n";
    file << "  \"budget_function_multiplier\": " << metadata_.budget_function_multiplier << ",\n";
    file << "  \"link_capacity_selection_threshold\": "
         << metadata_.link_capacity_selection_threshold << ",\n";
    file << "  \"metrics_config\": {\n";
    file << "    \"enable_trace\": " << (config_.enable_trace ? "true" : "false") << ",\n";
    file << "    \"enable_relative_gap\": " << (config_.enable_relative_gap ? "true" : "false")
         << ",\n";
    file << "    \"relative_gap_sample_period\": " << config_.relative_gap_sample_period << ",\n";
    file << "    \"flush_every_n_points\": " << config_.flush_every_n_points << ",\n";
    file << "    \"write_metadata_json\": " << (config_.write_metadata_json ? "true" : "false")
         << ",\n";
    file << "    \"write_summary_csv\": " << (config_.write_summary_csv ? "true" : "false")
         << "\n";
    file << "  }\n";
    file << "}\n";
  }

  void WriteSummary(const CndRunSummary& summary) const {
    const bool file_exists = std::filesystem::exists(summary_file_path_);
    std::ofstream file(summary_file_path_, std::ios::app);
    if (!file_exists) {
      file << "RunId,Timestamp,Scenario,Status,Dataset,Approach,Algorithm,DesignVariables,Links,ODPairs,"
              "MaxStandardIterations,MaxOptimalityIterations,StandardEvalCount,"
              "OptimalityConditionIterationCount,TotalTARunCount,InvalidTAStateCount,"
              "TARecoverySuccessCount,TARecoveryFailureCount,TotalElapsedTime(s),"
              "AvgTAComputeTime(s),MaxTAComputeTime(s),FinalObjective,FinalTotalTravelTime,"
              "FinalBudget,FinalBudgetViolation,BestFeasibleObjective,MaxBudgetViolation,"
              "RouteSearchThreads\n";
    }
    const std::string approach_name =
      metadata_.approach_name.empty() ? "UnknownApproach" : metadata_.approach_name;
    file << std::setprecision(10)
         << run_id_ << ","
         << CurrentTimestampIso8601() << ","
         << metadata_.scenario_name << ","
         << summary.status << ","
         << metadata_.dataset_name << ","
         << approach_name << ","
         << metadata_.nlopt_algorithm_name << ","
         << metadata_.design_variables << ","
         << metadata_.number_of_links << ","
         << metadata_.number_of_od_pairs << ","
         << metadata_.max_standard_iterations << ","
         << metadata_.max_optimality_condition_iterations << ","
         << summary.standard_eval_count << ","
         << summary.optimality_condition_iteration_count << ","
         << summary.total_ta_run_count << ","
         << summary.invalid_ta_state_count << ","
         << summary.ta_recovery_success_count << ","
         << summary.ta_recovery_failure_count << ","
         << summary.total_elapsed_seconds << ","
         << summary.avg_ta_compute_seconds << ","
         << summary.max_ta_compute_seconds << ","
         << summary.final_objective << ","
         << summary.final_total_travel_time << ","
         << summary.final_budget << ","
         << summary.final_budget_violation << ","
         << summary.best_feasible_objective << ","
         << summary.max_budget_violation << ","
         << metadata_.route_search_thread_count << "\n";
  }
};

}  // namespace TrafficAssignment

#endif  // CND_STATISTICS_RECORDER_H

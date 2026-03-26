#ifndef TOML_CONFIG_LOADER_H
#define TOML_CONFIG_LOADER_H

#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

#include <toml++/toml.hpp>

#include "ConfigUtils.h"
#include "../cnd/OptimizationStep.h"
#include "../cnd/CndStatisticsRecorder.h"

namespace TrafficAssignment::Config {

namespace fs = std::filesystem;

// ---------------------------------------------------------------------------
// Structured config types
// ---------------------------------------------------------------------------

struct NetworkConfig {
  std::string dataset = "SiouxFalls";
  std::string constraints_file;  // empty = auto-resolve from dataset
};

struct TapasConfig {
  long double mu = 0.0L;   // 0 = use default 0.5
  long double v = 0.0L;    // 0 = use default 0.25
  int pas_per_origin = 0;
  int pas_multiplier = 0;
  int rgap_check_interval = 0;
};

struct RouteBasedConfig {
  std::string shift_method = "NewtonStep";
  std::size_t route_search_threads = 1;
  int full_iteration_count = 0;
  int origin_iteration_count = 0;
  long double ema_alpha = 0.0L;
};

struct SolverConfig {
  std::string approach = "Tapas";
  long double approach_alpha = 1e-14L;
  int max_standard_iterations = 0;
  long double budget_function_multiplier = 5.0L;
  double budget_upper_bound = 100000.0;
  long double link_capacity_selection_threshold = 1e-3L;
  long double budget_threshold = 1e-1L;
  TapasConfig tapas;
  RouteBasedConfig route_based;
};

struct OutputConfig {
  bool verbose = true;
  std::string progress_format = "auto";
  bool print_effective_config = true;
  std::string quiet = "auto";
};

struct CndpConfig {
  NetworkConfig network;
  SolverConfig solver;
  OutputConfig output;
  CndMetricsConfig metrics;
  std::vector<OptimizationStepConfig> pipeline;
};

struct TapConfig {
  NetworkConfig network;
  SolverConfig solver;
  OutputConfig output;
  std::string output_root = "performance_results";
};

// ---------------------------------------------------------------------------
// TOML loading helpers (internal)
// ---------------------------------------------------------------------------

namespace detail {

inline void LoadNetworkConfig(const toml::table& tbl, NetworkConfig& net) {
  if (auto node = tbl["network"]["dataset"].value<std::string>()) {
    net.dataset = *node;
  }
  if (auto node = tbl["network"]["constraints_file"].value<std::string>()) {
    net.constraints_file = *node;
  }
}

inline void LoadSolverConfig(const toml::table& tbl, SolverConfig& solver) {
  const auto* s = tbl["solver"].as_table();
  if (!s) return;

  if (auto v = (*s)["approach"].value<std::string>()) solver.approach = *v;
  if (auto v = (*s)["approach_alpha"].value<double>()) solver.approach_alpha = static_cast<long double>(*v);
  if (auto v = (*s)["max_standard_iterations"].value<int64_t>()) solver.max_standard_iterations = static_cast<int>(*v);
  if (auto v = (*s)["budget_function_multiplier"].value<double>()) solver.budget_function_multiplier = static_cast<long double>(*v);
  if (auto v = (*s)["budget_upper_bound"].value<double>()) solver.budget_upper_bound = *v;
  if (auto v = (*s)["link_capacity_selection_threshold"].value<double>()) solver.link_capacity_selection_threshold = static_cast<long double>(*v);
  if (auto v = (*s)["budget_threshold"].value<double>()) solver.budget_threshold = static_cast<long double>(*v);

  // [solver.tapas]
  if (const auto* t = (*s)["tapas"].as_table()) {
    if (auto v = (*t)["mu"].value<double>()) solver.tapas.mu = static_cast<long double>(*v);
    if (auto v = (*t)["v"].value<double>()) solver.tapas.v = static_cast<long double>(*v);
    if (auto v = (*t)["pas_per_origin"].value<int64_t>()) solver.tapas.pas_per_origin = static_cast<int>(*v);
    if (auto v = (*t)["pas_multiplier"].value<int64_t>()) solver.tapas.pas_multiplier = static_cast<int>(*v);
    if (auto v = (*t)["rgap_check_interval"].value<int64_t>()) solver.tapas.rgap_check_interval = static_cast<int>(*v);
  }

  // [solver.route_based]
  if (const auto* r = (*s)["route_based"].as_table()) {
    if (auto v = (*r)["shift_method"].value<std::string>()) solver.route_based.shift_method = *v;
    if (auto v = (*r)["route_search_threads"].value<int64_t>()) solver.route_based.route_search_threads = static_cast<std::size_t>(*v);
    if (auto v = (*r)["full_iteration_count"].value<int64_t>()) solver.route_based.full_iteration_count = static_cast<int>(*v);
    if (auto v = (*r)["origin_iteration_count"].value<int64_t>()) solver.route_based.origin_iteration_count = static_cast<int>(*v);
    if (auto v = (*r)["ema_alpha"].value<double>()) solver.route_based.ema_alpha = static_cast<long double>(*v);
  }
}

inline void LoadOutputConfig(const toml::table& tbl, OutputConfig& out) {
  const auto* o = tbl["output"].as_table();
  if (!o) return;

  if (auto v = (*o)["verbose"].value<bool>()) out.verbose = *v;
  if (auto v = (*o)["progress_format"].value<std::string>()) out.progress_format = *v;
  if (auto v = (*o)["print_effective_config"].value<bool>()) out.print_effective_config = *v;
  if (auto v = (*o)["quiet"].value<std::string>()) out.quiet = *v;
}

inline void LoadMetricsConfig(const toml::table& tbl, CndMetricsConfig& m) {
  const auto* mt = tbl["metrics"].as_table();
  if (!mt) return;

  if (auto v = (*mt)["enable_trace"].value<bool>()) m.enable_trace = *v;
  if (auto v = (*mt)["enable_relative_gap"].value<bool>()) m.enable_relative_gap = *v;
  if (auto v = (*mt)["relative_gap_sample_period"].value<int64_t>()) m.relative_gap_sample_period = static_cast<int>(*v);
  if (auto v = (*mt)["flush_every_n_points"].value<int64_t>()) m.flush_every_n_points = static_cast<int>(*v);
  if (auto v = (*mt)["write_metadata_json"].value<bool>()) m.write_metadata_json = *v;
  if (auto v = (*mt)["write_summary_csv"].value<bool>()) m.write_summary_csv = *v;
  if (auto v = (*mt)["append_dataset_subdir"].value<bool>()) m.append_dataset_subdir = *v;
  if (auto v = (*mt)["output_root"].value<std::string>()) m.output_root = *v;
  if (auto v = (*mt)["run_id"].value<std::string>()) m.run_id = *v;
  if (auto v = (*mt)["scenario_name"].value<std::string>()) m.scenario_name = *v;
}

inline OptimizationStepConfig LoadStepFromTable(const toml::table& step) {
  OptimizationStepConfig cfg;

  if (auto v = step["type"].value<std::string>()) cfg.type = *v;
  if (auto v = step["name"].value<std::string>()) cfg.name = *v;
  if (auto v = step["max_iterations"].value<int64_t>()) cfg.max_iterations = static_cast<int>(*v);
  if (auto v = step["tolerance"].value<double>()) cfg.tolerance = *v;
  if (auto v = step["algorithm"].value<std::string>()) cfg.algorithm = *v;
  if (auto v = step["local_algorithm"].value<std::string>()) cfg.local_algorithm = *v;
  if (auto v = step["local_max_iterations"].value<int64_t>()) cfg.local_max_iterations = static_cast<int>(*v);
  if (auto v = step["local_tolerance"].value<double>()) cfg.local_tolerance = *v;
  if (auto v = step["population_size"].value<int64_t>()) cfg.population_size = static_cast<int>(*v);
  if (auto v = step["step_size"].value<double>()) cfg.step_size = *v;
  if (auto v = step["fd_epsilon"].value<double>()) cfg.fd_epsilon = *v;
  if (auto v = step["gradient_method"].value<std::string>()) cfg.gradient_method = *v;
  if (auto v = step["stochastic_optimizer"].value<std::string>()) cfg.stochastic_optimizer = *v;

  return cfg;
}

inline void LoadPipeline(const toml::table& tbl, std::vector<OptimizationStepConfig>& pipeline) {
  const auto* arr = tbl["pipeline"].as_array();
  if (!arr) return;

  for (const auto& elem : *arr) {
    const auto* step_tbl = elem.as_table();
    if (!step_tbl) continue;
    pipeline.push_back(LoadStepFromTable(*step_tbl));
  }
}

inline toml::table ParseTomlFile(const fs::path& path) {
  if (!fs::exists(path)) {
    throw std::runtime_error("Config file does not exist: " + path.string());
  }
  try {
    return toml::parse_file(path.string());
  } catch (const toml::parse_error& err) {
    throw std::runtime_error(
      "TOML parse error in " + path.string() + ": " + std::string(err.description())
    );
  }
}

}  // namespace detail

// ---------------------------------------------------------------------------
// Public loader functions
// ---------------------------------------------------------------------------

inline CndpConfig LoadCndpConfig(const fs::path& path) {
  auto tbl = detail::ParseTomlFile(path);
  CndpConfig config;

  detail::LoadNetworkConfig(tbl, config.network);
  detail::LoadSolverConfig(tbl, config.solver);
  detail::LoadOutputConfig(tbl, config.output);
  detail::LoadMetricsConfig(tbl, config.metrics);
  detail::LoadPipeline(tbl, config.pipeline);

  return config;
}

inline TapConfig LoadTapConfig(const fs::path& path) {
  auto tbl = detail::ParseTomlFile(path);
  TapConfig config;

  detail::LoadNetworkConfig(tbl, config.network);
  detail::LoadSolverConfig(tbl, config.solver);
  detail::LoadOutputConfig(tbl, config.output);

  // TAP-specific: output_root can live at top-level or under [solver]
  if (auto v = tbl["output_root"].value<std::string>()) {
    config.output_root = *v;
  } else if (auto v = tbl["solver"]["output_root"].value<std::string>()) {
    config.output_root = *v;
  }

  return config;
}

// ---------------------------------------------------------------------------
// Env var + CLI override functions (CndpConfig)
// ---------------------------------------------------------------------------

inline void ApplyEnvironmentOverrides(CndpConfig& config) {
  using namespace ConfigUtils;

  auto apply_str = [](const char* env_name, std::string& field) {
    if (const auto v = GetEnvValue(env_name)) field = *v;
  };
  auto apply_ld = [](const char* env_name, long double& field) {
    if (const auto v = GetEnvValue(env_name)) field = ParseNumber<long double>(*v, env_name);
  };
  auto apply_dbl = [](const char* env_name, double& field) {
    if (const auto v = GetEnvValue(env_name)) field = ParseNumber<double>(*v, env_name);
  };
  auto apply_int = [](const char* env_name, int& field) {
    if (const auto v = GetEnvValue(env_name)) field = ParseNumber<int>(*v, env_name);
  };
  auto apply_sizet = [](const char* env_name, std::size_t& field) {
    if (const auto v = GetEnvValue(env_name)) field = ParseNumber<std::size_t>(*v, env_name);
  };
  auto apply_bool = [](const char* env_name, bool& field) {
    if (const auto v = GetEnvValue(env_name)) field = ParseBool(*v, env_name);
  };

  // Network
  apply_str("CND_DATASET", config.network.dataset);
  apply_str("CND_CONSTRAINTS_FILE", config.network.constraints_file);

  // Solver
  apply_str("CND_APPROACH", config.solver.approach);
  apply_ld("CND_APPROACH_ALPHA", config.solver.approach_alpha);
  apply_int("CND_MAX_STANDARD_ITERATIONS", config.solver.max_standard_iterations);
  apply_ld("CND_BUDGET_MULTIPLIER", config.solver.budget_function_multiplier);
  apply_dbl("CND_BUDGET_UPPER_BOUND", config.solver.budget_upper_bound);
  apply_ld("CND_LINK_THRESHOLD", config.solver.link_capacity_selection_threshold);
  apply_ld("CND_BUDGET_THRESHOLD", config.solver.budget_threshold);

  // Solver.tapas
  apply_ld("CND_TAPAS_MU", config.solver.tapas.mu);
  apply_ld("CND_TAPAS_V", config.solver.tapas.v);
  apply_int("CND_PAS_PER_ORIGIN", config.solver.tapas.pas_per_origin);
  apply_int("CND_PAS_MULTIPLIER", config.solver.tapas.pas_multiplier);
  apply_int("CND_RGAP_CHECK_INTERVAL", config.solver.tapas.rgap_check_interval);

  // Solver.route_based
  apply_str("CND_SHIFT_METHOD", config.solver.route_based.shift_method);
  apply_sizet("CND_ROUTE_THREADS", config.solver.route_based.route_search_threads);
  apply_int("CND_FULL_ITERATION_COUNT", config.solver.route_based.full_iteration_count);
  apply_int("CND_ORIGIN_ITERATION_COUNT", config.solver.route_based.origin_iteration_count);
  apply_ld("CND_EMA_ALPHA", config.solver.route_based.ema_alpha);

  // Output
  apply_bool("CND_LOADER_VERBOSE", config.output.verbose);
  apply_str("CND_PROGRESS_FORMAT", config.output.progress_format);
  apply_bool("CND_PRINT_CONFIG", config.output.print_effective_config);

  // Metrics
  apply_bool("CND_METRICS_ENABLE_TRACE", config.metrics.enable_trace);
  apply_bool("CND_METRICS_ENABLE_RELATIVE_GAP", config.metrics.enable_relative_gap);
  apply_int("CND_METRICS_RELATIVE_GAP_PERIOD", config.metrics.relative_gap_sample_period);
  apply_int("CND_METRICS_FLUSH_EVERY", config.metrics.flush_every_n_points);
  apply_bool("CND_METRICS_WRITE_METADATA", config.metrics.write_metadata_json);
  apply_bool("CND_METRICS_WRITE_SUMMARY", config.metrics.write_summary_csv);
  apply_str("CND_METRICS_OUTPUT_ROOT", config.metrics.output_root);
  apply_str("CND_METRICS_RUN_ID", config.metrics.run_id);
  apply_str("CND_METRICS_SCENARIO_NAME", config.metrics.scenario_name);

  if (const auto v = GetEnvValue("CND_METRICS_NO_DATASET_SUBDIR")) {
    config.metrics.append_dataset_subdir = !ParseBool(*v, "CND_METRICS_NO_DATASET_SUBDIR");
  }
}

inline void ApplyCliOverrides(const ConfigUtils::CliOptions& cli, CndpConfig& config) {
  using namespace ConfigUtils;

  for (const auto& [raw_key, value] : cli.values) {
    if (raw_key == "config" || raw_key == "help" || raw_key == "h") continue;

    const std::string key = NormalizeKey(raw_key);

    // Metrics overrides (prefixed with metrics_)
    if (key.rfind("metrics_", 0) == 0) {
      const std::string mkey = key.substr(8);  // strip "metrics_"
      if (mkey == "enable_trace") { config.metrics.enable_trace = ParseBool(value, key); }
      else if (mkey == "enable_relative_gap") { config.metrics.enable_relative_gap = ParseBool(value, key); }
      else if (mkey == "relative_gap_sample_period") { config.metrics.relative_gap_sample_period = ParseNumber<int>(value, key); }
      else if (mkey == "flush_every_n_points") { config.metrics.flush_every_n_points = ParseNumber<int>(value, key); }
      else if (mkey == "write_metadata_json") { config.metrics.write_metadata_json = ParseBool(value, key); }
      else if (mkey == "write_summary_csv") { config.metrics.write_summary_csv = ParseBool(value, key); }
      else if (mkey == "output_root") { config.metrics.output_root = value; }
      else if (mkey == "run_id") { config.metrics.run_id = value; }
      else if (mkey == "scenario_name" || mkey == "scenario") { config.metrics.scenario_name = value; }
      else if (mkey == "no_dataset_subdir") { config.metrics.append_dataset_subdir = !ParseBool(value, key); }
      else { throw std::runtime_error("Unknown metrics option: '" + raw_key + "'"); }
      continue;
    }

    // Network
    if (key == "dataset") { config.network.dataset = value; }
    else if (key == "constraints_file" || key == "constraints_path" || key == "constraints") { config.network.constraints_file = value; }
    // Solver
    else if (key == "approach" || key == "approach_type") { config.solver.approach = value; }
    else if (key == "approach_alpha" || key == "alpha") { config.solver.approach_alpha = ParseNumber<long double>(value, key); }
    else if (key == "max_standard_iterations" || key == "max_iters") { config.solver.max_standard_iterations = ParseNumber<int>(value, key); }
    else if (key == "budget_function_multiplier" || key == "budget_multiplier") { config.solver.budget_function_multiplier = ParseNumber<long double>(value, key); }
    else if (key == "budget_upper_bound") { config.solver.budget_upper_bound = ParseNumber<double>(value, key); }
    else if (key == "link_capacity_selection_threshold" || key == "link_threshold") { config.solver.link_capacity_selection_threshold = ParseNumber<long double>(value, key); }
    else if (key == "budget_threshold") { config.solver.budget_threshold = ParseNumber<long double>(value, key); }
    // Solver.tapas
    else if (key == "tapas_mu" || key == "mu") { config.solver.tapas.mu = ParseNumber<long double>(value, key); }
    else if (key == "tapas_v" || key == "v") { config.solver.tapas.v = ParseNumber<long double>(value, key); }
    else if (key == "pas_per_origin") { config.solver.tapas.pas_per_origin = ParseNumber<int>(value, key); }
    else if (key == "pas_multiplier") { config.solver.tapas.pas_multiplier = ParseNumber<int>(value, key); }
    else if (key == "rgap_check_interval") { config.solver.tapas.rgap_check_interval = ParseNumber<int>(value, key); }
    // Solver.route_based
    else if (key == "shift_method") { config.solver.route_based.shift_method = value; }
    else if (key == "route_search_threads" || key == "route_threads") { config.solver.route_based.route_search_threads = ParseNumber<std::size_t>(value, key); }
    else if (key == "full_iteration_count") { config.solver.route_based.full_iteration_count = ParseNumber<int>(value, key); }
    else if (key == "origin_iteration_count") { config.solver.route_based.origin_iteration_count = ParseNumber<int>(value, key); }
    else if (key == "ema_alpha") { config.solver.route_based.ema_alpha = ParseNumber<long double>(value, key); }
    // Output
    else if (key == "loader_verbose") { config.output.verbose = ParseBool(value, key); }
    else if (key == "progress_format") { config.output.progress_format = value; }
    else if (key == "print_effective_config" || key == "print_config") { config.output.print_effective_config = ParseBool(value, key); }
    else if (key == "quiet") { config.output.quiet = value; }
    else { throw std::runtime_error("Unknown option: '" + raw_key + "'"); }
  }
}

// ---------------------------------------------------------------------------
// Env var + CLI override functions (TapConfig)
// ---------------------------------------------------------------------------

inline void ApplyEnvironmentOverrides(TapConfig& config) {
  using namespace ConfigUtils;

  auto apply_str = [](const char* env_name, std::string& field) {
    if (const auto v = GetEnvValue(env_name)) field = *v;
  };
  auto apply_ld = [](const char* env_name, long double& field) {
    if (const auto v = GetEnvValue(env_name)) field = ParseNumber<long double>(*v, env_name);
  };
  auto apply_int = [](const char* env_name, int& field) {
    if (const auto v = GetEnvValue(env_name)) field = ParseNumber<int>(*v, env_name);
  };
  auto apply_sizet = [](const char* env_name, std::size_t& field) {
    if (const auto v = GetEnvValue(env_name)) field = ParseNumber<std::size_t>(*v, env_name);
  };
  auto apply_bool = [](const char* env_name, bool& field) {
    if (const auto v = GetEnvValue(env_name)) field = ParseBool(*v, env_name);
  };

  apply_str("CND_DATASET", config.network.dataset);
  apply_str("CND_APPROACH", config.solver.approach);
  apply_ld("CND_APPROACH_ALPHA", config.solver.approach_alpha);
  apply_int("CND_MAX_STANDARD_ITERATIONS", config.solver.max_standard_iterations);
  apply_str("CND_SHIFT_METHOD", config.solver.route_based.shift_method);
  apply_sizet("CND_ROUTE_THREADS", config.solver.route_based.route_search_threads);
  apply_int("CND_FULL_ITERATION_COUNT", config.solver.route_based.full_iteration_count);
  apply_int("CND_ORIGIN_ITERATION_COUNT", config.solver.route_based.origin_iteration_count);
  apply_ld("CND_EMA_ALPHA", config.solver.route_based.ema_alpha);
  apply_ld("CND_TAPAS_MU", config.solver.tapas.mu);
  apply_ld("CND_TAPAS_V", config.solver.tapas.v);
  apply_int("CND_PAS_PER_ORIGIN", config.solver.tapas.pas_per_origin);
  apply_int("CND_PAS_MULTIPLIER", config.solver.tapas.pas_multiplier);
  apply_int("CND_RGAP_CHECK_INTERVAL", config.solver.tapas.rgap_check_interval);
  apply_str("CND_OUTPUT_ROOT", config.output_root);
  apply_bool("CND_PRINT_CONFIG", config.output.print_effective_config);
}

inline void ApplyCliOverrides(const ConfigUtils::CliOptions& cli, TapConfig& config) {
  using namespace ConfigUtils;

  for (const auto& [raw_key, value] : cli.values) {
    if (raw_key == "config" || raw_key == "help" || raw_key == "h") continue;

    const std::string key = NormalizeKey(raw_key);

    if (key == "dataset") { config.network.dataset = value; }
    else if (key == "approach" || key == "approach_type") { config.solver.approach = value; }
    else if (key == "approach_alpha" || key == "alpha") { config.solver.approach_alpha = ParseNumber<long double>(value, key); }
    else if (key == "max_standard_iterations" || key == "max_iters") { config.solver.max_standard_iterations = ParseNumber<int>(value, key); }
    else if (key == "shift_method") { config.solver.route_based.shift_method = value; }
    else if (key == "route_search_threads" || key == "route_threads") { config.solver.route_based.route_search_threads = ParseNumber<std::size_t>(value, key); }
    else if (key == "full_iteration_count") { config.solver.route_based.full_iteration_count = ParseNumber<int>(value, key); }
    else if (key == "origin_iteration_count") { config.solver.route_based.origin_iteration_count = ParseNumber<int>(value, key); }
    else if (key == "ema_alpha") { config.solver.route_based.ema_alpha = ParseNumber<long double>(value, key); }
    else if (key == "tapas_mu" || key == "mu") { config.solver.tapas.mu = ParseNumber<long double>(value, key); }
    else if (key == "tapas_v" || key == "v") { config.solver.tapas.v = ParseNumber<long double>(value, key); }
    else if (key == "pas_per_origin") { config.solver.tapas.pas_per_origin = ParseNumber<int>(value, key); }
    else if (key == "pas_multiplier") { config.solver.tapas.pas_multiplier = ParseNumber<int>(value, key); }
    else if (key == "rgap_check_interval") { config.solver.tapas.rgap_check_interval = ParseNumber<int>(value, key); }
    else if (key == "output_root") { config.output_root = value; }
    else if (key == "print_effective_config" || key == "print_config") { config.output.print_effective_config = ParseBool(value, key); }
    else if (key == "quiet") { config.output.quiet = value; }
    else { throw std::runtime_error("Unknown option: '" + raw_key + "'"); }
  }
}

}  // namespace TrafficAssignment::Config

#endif  // TOML_CONFIG_LOADER_H

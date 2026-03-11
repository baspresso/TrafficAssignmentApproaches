#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <regex>
#include <stdexcept>
#include <string>

#include "../include/common/ConfigUtils.h"
#include "../include/cnd/BilevelCND.h"
#include "../include/cnd/DirectedConstraintLoader.h"
#include "../include/cnd/OptimizationStep.h"
#include "../include/tap/algorithms/route_based/RouteBasedApproach.h"
#include "../include/tap/algorithms/tapas_based/TapasApproach.h"
#include "../include/tap/core/NetworkBuilder.h"

namespace fs = std::filesystem;

using namespace TrafficAssignment::ConfigUtils;

namespace {

struct RunConfig {
  std::string dataset = "SiouxFalls";
  std::string constraints_file;
  std::string approach = "RouteBased";
  long double approach_alpha = 1e-14L;
  std::string shift_method = "NewtonStep";
  std::size_t route_search_threads = 1;

  // TAP approach tuning (0 = use approach default)
  int max_standard_iterations = 0;
  int full_iteration_count = 0;        // RouteBased-specific
  int origin_iteration_count = 0;      // RouteBased-specific
  long double ema_alpha = 0.0L;        // RouteBased-specific
  int pas_per_origin = 0;              // Tapas-specific
  int pas_multiplier = 0;              // Tapas-specific
  int rgap_check_interval = 0;         // Tapas-specific (0 = use default)

  long double link_capacity_selection_threshold = 1e-3L;
  long double budget_threshold = 1e-1L;
  long double budget_function_multiplier = 5.0L;
  double budget_upper_bound = 100000.0;

  bool loader_verbose = true;
  bool print_effective_config = true;
  TrafficAssignment::CndMetricsConfig metrics;

  // Pipeline step configs parsed from [step.N] sections
  std::map<int, TrafficAssignment::OptimizationStepConfig> step_sections;
};


void ApplyRunOption(RunConfig& config,
                    const std::string& raw_key,
                    const std::string& raw_value) {
  const std::string key = NormalizeKey(raw_key);
  if (key == "dataset") {
    config.dataset = raw_value;
    return;
  }
  if (key == "constraints_file" || key == "constraints_path" || key == "constraints") {
    config.constraints_file = raw_value;
    return;
  }
  if (key == "approach" || key == "approach_type") {
    config.approach = raw_value;
    return;
  }
  if (key == "approach_alpha" || key == "alpha") {
    config.approach_alpha = ParseNumber<long double>(raw_value, key);
    return;
  }
  if (key == "shift_method") {
    config.shift_method = raw_value;
    return;
  }
  if (key == "route_search_threads" || key == "route_threads") {
    config.route_search_threads = ParseNumber<std::size_t>(raw_value, key);
    return;
  }
  if (key == "max_standard_iterations" || key == "max_iters") {
    config.max_standard_iterations = ParseNumber<int>(raw_value, key);
    return;
  }
  if (key == "full_iteration_count") {
    config.full_iteration_count = ParseNumber<int>(raw_value, key);
    return;
  }
  if (key == "origin_iteration_count") {
    config.origin_iteration_count = ParseNumber<int>(raw_value, key);
    return;
  }
  if (key == "ema_alpha") {
    config.ema_alpha = ParseNumber<long double>(raw_value, key);
    return;
  }
  if (key == "pas_per_origin") {
    config.pas_per_origin = ParseNumber<int>(raw_value, key);
    return;
  }
  if (key == "pas_multiplier") {
    config.pas_multiplier = ParseNumber<int>(raw_value, key);
    return;
  }
  if (key == "rgap_check_interval") {
    config.rgap_check_interval = ParseNumber<int>(raw_value, key);
    return;
  }
  if (key == "link_capacity_selection_threshold" || key == "link_threshold") {
    config.link_capacity_selection_threshold = ParseNumber<long double>(raw_value, key);
    return;
  }
  if (key == "budget_threshold") {
    config.budget_threshold = ParseNumber<long double>(raw_value, key);
    return;
  }
  if (key == "budget_function_multiplier" || key == "budget_multiplier") {
    config.budget_function_multiplier = ParseNumber<long double>(raw_value, key);
    return;
  }
  if (key == "budget_upper_bound") {
    config.budget_upper_bound = ParseNumber<double>(raw_value, key);
    return;
  }
  if (key == "loader_verbose") {
    config.loader_verbose = ParseBool(raw_value, key);
    return;
  }
  if (key == "print_effective_config" || key == "print_config") {
    config.print_effective_config = ParseBool(raw_value, key);
    return;
  }

  throw std::runtime_error("Unknown run option: '" + raw_key + "'");
}

void ApplyMetricsOption(RunConfig& config,
                        const std::string& raw_key,
                        const std::string& raw_value) {
  const std::string key = NormalizeKey(raw_key);
  if (key == "enable_trace") {
    config.metrics.enable_trace = ParseBool(raw_value, key);
    return;
  }
  if (key == "enable_relative_gap") {
    config.metrics.enable_relative_gap = ParseBool(raw_value, key);
    return;
  }
  if (key == "relative_gap_sample_period") {
    config.metrics.relative_gap_sample_period = ParseNumber<int>(raw_value, key);
    return;
  }
  if (key == "flush_every_n_points") {
    config.metrics.flush_every_n_points = ParseNumber<int>(raw_value, key);
    return;
  }
  if (key == "write_metadata_json") {
    config.metrics.write_metadata_json = ParseBool(raw_value, key);
    return;
  }
  if (key == "write_summary_csv") {
    config.metrics.write_summary_csv = ParseBool(raw_value, key);
    return;
  }
  if (key == "output_root") {
    config.metrics.output_root = raw_value;
    return;
  }
  if (key == "run_id") {
    config.metrics.run_id = raw_value;
    return;
  }
  if (key == "scenario_name" || key == "scenario") {
    config.metrics.scenario_name = raw_value;
    return;
  }
  if (key == "no_dataset_subdir") {
    config.metrics.append_dataset_subdir = !ParseBool(raw_value, key);
    return;
  }

  throw std::runtime_error("Unknown metrics option: '" + raw_key + "'");
}

void ApplyStepOption(TrafficAssignment::OptimizationStepConfig& step,
                     const std::string& raw_key,
                     const std::string& raw_value) {
  const std::string key = NormalizeKey(raw_key);
  if (key == "type") {
    step.type = raw_value;
  } else if (key == "name") {
    step.name = raw_value;
  } else if (key == "max_iterations" || key == "max_iters") {
    step.max_iterations = ParseNumber<int>(raw_value, key);
  } else if (key == "tolerance") {
    step.tolerance = ParseNumber<double>(raw_value, key);
  } else if (key == "algorithm") {
    step.algorithm = raw_value;
  } else if (key == "local_algorithm") {
    step.local_algorithm = raw_value;
  } else if (key == "local_max_iterations" || key == "local_max_iters") {
    step.local_max_iterations = ParseNumber<int>(raw_value, key);
  } else if (key == "local_tolerance") {
    step.local_tolerance = ParseNumber<double>(raw_value, key);
  } else if (key == "population_size") {
    step.population_size = ParseNumber<int>(raw_value, key);
  } else {
    throw std::runtime_error("Unknown step option: '" + raw_key + "'");
  }
}

fs::path ResolveConstraintsPath(const RunConfig& config, const fs::path& project_root) {
  if (config.constraints_file.empty()) {
    return project_root / "data" / "TransportationNetworks" / config.dataset /
           (config.dataset + "_constraints.csv");
  }
  return ResolvePath(fs::path(config.constraints_file), project_root);
}

void ApplyConfigFile(const fs::path& config_path, RunConfig& config) {
  if (!fs::exists(config_path)) {
    throw std::runtime_error("Config file does not exist: " + config_path.string());
  }

  std::ifstream file(config_path);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open config file: " + config_path.string());
  }

  std::string line;
  std::string current_section;
  int current_step_index = -1;
  int line_number = 0;
  const std::regex step_pattern(R"(step\.(\d+))");
  while (std::getline(file, line)) {
    ++line_number;
    const std::string trimmed = TrimCopy(line);
    if (trimmed.empty() || trimmed[0] == '#' || trimmed[0] == ';') {
      continue;
    }

    if (trimmed.front() == '[' && trimmed.back() == ']') {
      current_section = ToLowerCopy(TrimCopy(trimmed.substr(1, trimmed.size() - 2)));
      current_step_index = -1;
      std::smatch match;
      if (std::regex_match(current_section, match, step_pattern)) {
        current_step_index = std::stoi(match[1].str());
        if (config.step_sections.find(current_step_index) == config.step_sections.end()) {
          config.step_sections[current_step_index] = TrafficAssignment::OptimizationStepConfig();
        }
      }
      continue;
    }

    const auto eq_pos = trimmed.find('=');
    if (eq_pos == std::string::npos) {
      throw std::runtime_error(
        "Invalid config line " + std::to_string(line_number) + " in " +
        config_path.string() + ": expected key=value format."
      );
    }

    const std::string key = TrimCopy(trimmed.substr(0, eq_pos));
    const std::string value = TrimCopy(trimmed.substr(eq_pos + 1));
    if (current_section == "run") {
      ApplyRunOption(config, key, value);
    } else if (current_section == "metrics") {
      ApplyMetricsOption(config, key, value);
    } else if (current_step_index >= 0) {
      ApplyStepOption(config.step_sections[current_step_index], key, value);
    } else {
      throw std::runtime_error(
        "Unsupported config section '" + current_section + "' at line " +
        std::to_string(line_number) + " in " + config_path.string()
      );
    }
  }
}

void ApplyEnvironmentOverrides(RunConfig& config) {
  auto apply_run = [&config](const char* env_name, const std::string& key) {
    if (const auto env_value = GetEnvValue(env_name)) {
      ApplyRunOption(config, key, *env_value);
    }
  };
  auto apply_metrics = [&config](const char* env_name, const std::string& key) {
    if (const auto env_value = GetEnvValue(env_name)) {
      ApplyMetricsOption(config, key, *env_value);
    }
  };

  apply_run("CND_DATASET", "dataset");
  apply_run("CND_CONSTRAINTS_FILE", "constraints_file");
  apply_run("CND_APPROACH", "approach");
  apply_run("CND_APPROACH_ALPHA", "approach_alpha");
  apply_run("CND_SHIFT_METHOD", "shift_method");
  apply_run("CND_ROUTE_THREADS", "route_search_threads");
  apply_run("CND_MAX_STANDARD_ITERATIONS", "max_standard_iterations");
  apply_run("CND_FULL_ITERATION_COUNT", "full_iteration_count");
  apply_run("CND_ORIGIN_ITERATION_COUNT", "origin_iteration_count");
  apply_run("CND_EMA_ALPHA", "ema_alpha");
  apply_run("CND_PAS_PER_ORIGIN", "pas_per_origin");
  apply_run("CND_PAS_MULTIPLIER", "pas_multiplier");
  apply_run("CND_RGAP_CHECK_INTERVAL", "rgap_check_interval");
  apply_run("CND_LINK_THRESHOLD", "link_capacity_selection_threshold");
  apply_run("CND_BUDGET_THRESHOLD", "budget_threshold");
  apply_run("CND_BUDGET_MULTIPLIER", "budget_function_multiplier");
  apply_run("CND_BUDGET_UPPER_BOUND", "budget_upper_bound");
  apply_run("CND_LOADER_VERBOSE", "loader_verbose");
  apply_run("CND_PRINT_CONFIG", "print_effective_config");

  apply_metrics("CND_METRICS_ENABLE_TRACE", "enable_trace");
  apply_metrics("CND_METRICS_ENABLE_RELATIVE_GAP", "enable_relative_gap");
  apply_metrics("CND_METRICS_RELATIVE_GAP_PERIOD", "relative_gap_sample_period");
  apply_metrics("CND_METRICS_FLUSH_EVERY", "flush_every_n_points");
  apply_metrics("CND_METRICS_WRITE_METADATA", "write_metadata_json");
  apply_metrics("CND_METRICS_WRITE_SUMMARY", "write_summary_csv");
  apply_metrics("CND_METRICS_OUTPUT_ROOT", "output_root");
  apply_metrics("CND_METRICS_RUN_ID", "run_id");
  apply_metrics("CND_METRICS_SCENARIO_NAME", "scenario_name");
  apply_metrics("CND_METRICS_NO_DATASET_SUBDIR", "no_dataset_subdir");
}

void ApplyCliOverrides(const CliOptions& cli, RunConfig& config) {
  for (const auto& [raw_key, value] : cli.values) {
    if (raw_key == "config" || raw_key == "help" || raw_key == "h") {
      continue;
    }
    if (raw_key.rfind("metrics_", 0) == 0) {
      ApplyMetricsOption(config, raw_key.substr(std::string("metrics_").size()), value);
      continue;
    }
    ApplyRunOption(config, raw_key, value);
  }
}

std::shared_ptr<TrafficAssignment::TrafficAssignmentApproach<long double>>
CreateApproach(const RunConfig& config, TrafficAssignment::Network<long double>& network) {
  const std::string approach = ToLowerCopy(config.approach);
  if (approach == "routebased" || approach == "route_based") {
    return std::make_shared<TrafficAssignment::RouteBasedApproach<long double>>(
      network,
      config.approach_alpha,
      config.shift_method,
      config.route_search_threads,
      config.max_standard_iterations > 0 ? config.max_standard_iterations : 200,
      config.full_iteration_count > 0 ? config.full_iteration_count : 3,
      config.origin_iteration_count > 0 ? config.origin_iteration_count : 1,
      config.ema_alpha > 0.0L ? config.ema_alpha : 0.7L
    );
  }
  if (approach == "tapas") {
    return std::make_shared<TrafficAssignment::TapasApproach<long double>>(
      network,
      config.approach_alpha,
      config.shift_method,
      config.max_standard_iterations > 0 ? config.max_standard_iterations : 20,
      config.pas_per_origin > 0 ? config.pas_per_origin : 1,
      config.pas_multiplier > 0 ? config.pas_multiplier : 5,
      config.rgap_check_interval > 0 ? config.rgap_check_interval : 5
    );
  }
  throw std::runtime_error(
    "Unsupported approach '" + config.approach +
    "'. Supported: RouteBased, Tapas."
  );
}

void PrintHelp() {
  std::cout
    << "Usage: main.exe --config <path> [options]\n\n"
    << "Layering: defaults -> config file -> environment -> CLI\n\n"
    << "Optimization steps are defined via [step.N] sections in the INI config file.\n"
    << "Example:\n"
    << "  [step.0]\n"
    << "  type = nlopt\n"
    << "  algorithm = LN_COBYLA\n"
    << "  max_iterations = 200\n"
    << "  tolerance = 1e-4\n\n"
    << "  [step.1]\n"
    << "  type = optimality_condition\n"
    << "  max_iterations = 30\n\n"
    << "Step options:\n"
    << "  type                             nlopt | optimality_condition | optimlib\n"
    << "  algorithm                        NLopt: LN_COBYLA, GN_ISRES, ...\n"
    << "                                   OptimLib: DE, PSO, NM, DE_PRMM, PSO_DV, GD\n"
    << "  local_algorithm                  Local sub-algorithm for AUGLAG (NLopt)\n"
    << "  max_iterations                   Max iterations/generations for this step\n"
    << "  tolerance                        Convergence tolerance\n"
    << "  population_size                  Population size for OptimLib (0=auto)\n"
    << "  name                             Display name for this step\n\n"
    << "General options (CLI or [run] section):\n"
    << "  --config <path>                  INI config file path (required)\n"
    << "  --dataset <name>                 Dataset name (e.g. SiouxFalls)\n"
    << "  --constraints-file <path>        Constraints CSV path\n"
    << "  --approach <RouteBased|Tapas>\n"
    << "  --shift-method <name>            Shift method (RouteBased/Tapas)\n"
    << "  --route-threads <n>              Route search thread count\n"
    << "  --max-standard-iterations <n>    TAP iteration limit (default: 200/20)\n"
    << "  --full-iteration-count <n>       RouteBased OD queue passes (default: 3)\n"
    << "  --origin-iteration-count <n>     RouteBased per-origin passes (default: 1)\n"
    << "  --ema-alpha <value>              RouteBased EMA smoothing (default: 0.7)\n"
    << "  --pas-per-origin <n>             Tapas PAS per origin (default: 1)\n"
    << "  --pas-multiplier <n>             Tapas PAS multiplier (default: 5)\n"
    << "  --link-threshold <value>\n"
    << "  --budget-threshold <value>\n"
    << "  --budget-multiplier <value>\n"
    << "  --budget-upper-bound <value>\n"
    << "  --loader-verbose <true|false>\n"
    << "  --print-config <true|false>\n\n"
    << "Metrics CLI overrides (prefix metrics_):\n"
    << "  --metrics_enable_trace <true|false>\n"
    << "  --metrics_enable_relative_gap <true|false>\n"
    << "  --metrics_relative_gap_sample_period <n>\n"
    << "  --metrics_flush_every_n_points <n>\n"
    << "  --metrics_write_metadata_json <true|false>\n"
    << "  --metrics_write_summary_csv <true|false>\n"
    << "  --metrics_output_root <path>\n"
    << "  --metrics_run_id <id>\n"
    << "  --metrics_scenario_name <name>\n"
    << "  --metrics_no_dataset_subdir <true|false>  Skip dataset subdirectory in output\n";
}

void PrintEffectiveConfig(const RunConfig& config,
                          const fs::path& project_root,
                          const fs::path& constraints_path,
                          const std::optional<fs::path>& config_path) {
  std::cout << "\nEffective run configuration:\n";
  std::cout << "  project_root: " << project_root.string() << '\n';
  std::cout << "  config_file: " << (config_path ? config_path->string() : "<none>") << '\n';
  std::cout << "  dataset: " << config.dataset << '\n';
  std::cout << "  constraints_file: " << constraints_path.string() << '\n';
  std::cout << "  approach: " << config.approach << '\n';
  std::cout << "  shift_method: " << config.shift_method << '\n';
  std::cout << "  route_search_threads: " << config.route_search_threads << '\n';
  if (config.max_standard_iterations > 0)
    std::cout << "  max_standard_iterations: " << config.max_standard_iterations << '\n';
  if (config.full_iteration_count > 0)
    std::cout << "  full_iteration_count: " << config.full_iteration_count << '\n';
  if (config.origin_iteration_count > 0)
    std::cout << "  origin_iteration_count: " << config.origin_iteration_count << '\n';
  if (config.ema_alpha > 0.0L)
    std::cout << "  ema_alpha: " << config.ema_alpha << '\n';
  if (config.pas_per_origin > 0)
    std::cout << "  pas_per_origin: " << config.pas_per_origin << '\n';
  if (config.pas_multiplier > 0)
    std::cout << "  pas_multiplier: " << config.pas_multiplier << '\n';
  if (config.rgap_check_interval > 0)
    std::cout << "  rgap_check_interval: " << config.rgap_check_interval << '\n';
  std::cout << "  budget_upper_bound: " << config.budget_upper_bound << '\n';
  std::cout << "  metrics.output_root: " << config.metrics.output_root << '\n';
  if (!config.metrics.append_dataset_subdir) {
    std::cout << "  metrics.no_dataset_subdir: true\n";
  }
  if (!config.metrics.scenario_name.empty()) {
    std::cout << "  metrics.scenario_name: " << config.metrics.scenario_name << '\n';
  }
  std::cout << "  pipeline steps: " << config.step_sections.size() << '\n';
  for (const auto& [index, step] : config.step_sections) {
    std::cout << "    [step." << index << "] type=" << step.type;
    if (!step.algorithm.empty()) std::cout << " algorithm=" << step.algorithm;
    std::cout << " max_iter=" << step.max_iterations;
    if (step.tolerance > 0) std::cout << " tol=" << step.tolerance;
    if (step.population_size > 0) std::cout << " pop_size=" << step.population_size;
    std::cout << '\n';
  }
  std::cout << std::endl;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const CliOptions cli = ParseCliOptions(argc, argv);
    if (cli.help_requested) {
      PrintHelp();
      return 0;
    }

    const fs::path project_root = FindProjectRoot();
    RunConfig config;

    std::optional<fs::path> config_path;
    const auto cli_config_it = cli.values.find("config");
    if (cli_config_it != cli.values.end()) {
      config_path = ResolvePath(fs::path(cli_config_it->second), project_root);
    } else if (const auto env_config = GetEnvValue("CND_CONFIG")) {
      config_path = ResolvePath(fs::path(*env_config), project_root);
    }

    if (config_path.has_value()) {
      ApplyConfigFile(*config_path, config);
    }

    ApplyEnvironmentOverrides(config);
    ApplyCliOverrides(cli, config);

    const fs::path constraints_path = ResolveConstraintsPath(config, project_root);
    if (!fs::exists(constraints_path)) {
      throw std::runtime_error(
        "Constraint file does not exist: " + constraints_path.string()
      );
    }

    if (config.print_effective_config) {
      PrintEffectiveConfig(config, project_root, constraints_path, config_path);
    }

    TrafficAssignment::NetworkBuilder builder;
    auto network = builder.BuildFromDataset<long double>(config.dataset);
    auto approach = CreateApproach(config, network);

    TrafficAssignment::DirectedConstraintLoader loader;
    loader.SetVerbose(config.loader_verbose);
    auto constraints = loader.LoadFromFile(constraints_path.string());

    if (config.step_sections.empty()) {
      throw std::runtime_error(
        "No [step.N] sections found in config. At least one optimization step "
        "must be defined. Example:\n"
        "  [step.0]\n"
        "  type = nlopt\n"
        "  algorithm = LN_COBYLA\n"
        "  max_iterations = 100\n"
      );
    }

    std::vector<TrafficAssignment::OptimizationStepConfig> steps;
    for (const auto& [index, step_config] : config.step_sections) {
      steps.push_back(step_config);
    }

    std::cout << "Creating BilevelCND solver..." << std::endl;
    TrafficAssignment::BilevelCND<long double> cnd(
      network,
      approach,
      constraints,
      steps,
      config.link_capacity_selection_threshold,
      config.budget_threshold,
      config.budget_function_multiplier,
      config.budget_upper_bound,
      config.metrics,
      config.route_search_threads
    );

    std::cout << "      BilevelCND solver created" << std::endl;
    std::cout << "      Design variables: " << constraints.size() << std::endl;

    std::cout << "\nRunning bilevel optimization..." << std::endl;

    cnd.ComputeNetworkDesign();
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "Fatal error: " << ex.what() << std::endl;
  } catch (...) {
    std::cerr << "Fatal error: unknown exception" << std::endl;
  }

  return 1;
}

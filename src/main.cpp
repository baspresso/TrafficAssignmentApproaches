#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include <Eigen/Dense>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <nlopt.hpp>

#include "../include/cnd/BilevelCND.h"
#include "../include/cnd/DirectedConstraintLoader.h"
#include "../include/tap/algorithms/demand_based/PumpOutDemandBasedApproach.h"
#include "../include/tap/algorithms/route_based/RouteBasedApproach.h"
#include "../include/tap/algorithms/tapas_based/TapasApproach.h"
#include "../include/tap/core/NetworkBuilder.h"

namespace fs = std::filesystem;

namespace {

struct CliOptions {
  std::unordered_map<std::string, std::string> values;
  bool help_requested = false;
};

struct RunConfig {
  std::string dataset = "SiouxFalls";
  std::string constraints_file;
  std::string approach = "RouteBased";
  long double approach_alpha = 1e-14L;
  std::string shift_method = "NewtonStep";
  std::size_t route_search_threads = 1;

  int max_standard_iterations = 100;
  int max_optimality_condition_iterations = 10;
  long double optimization_tolerance = 1e-4L;
  long double link_capacity_selection_threshold = 1e-3L;
  long double budget_threshold = 1e-1L;
  long double budget_function_multiplier = 5.0L;
  double budget_upper_bound = 100000.0;
  std::string nlopt_algorithm = "LN_COBYLA";

  bool loader_verbose = true;
  bool print_effective_config = true;
  TrafficAssignment::CndMetricsConfig metrics;
};

std::string ToLowerCopy(std::string text) {
  std::transform(
    text.begin(),
    text.end(),
    text.begin(),
    [](unsigned char c) { return static_cast<char>(std::tolower(c)); }
  );
  return text;
}

std::string ToUpperCopy(std::string text) {
  std::transform(
    text.begin(),
    text.end(),
    text.begin(),
    [](unsigned char c) { return static_cast<char>(std::toupper(c)); }
  );
  return text;
}

std::string NormalizeKey(std::string key) {
  key = ToLowerCopy(std::move(key));
  std::replace(key.begin(), key.end(), '-', '_');
  return key;
}

std::string TrimCopy(const std::string& text) {
  const auto first = text.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) {
    return "";
  }
  const auto last = text.find_last_not_of(" \t\r\n");
  return text.substr(first, last - first + 1);
}

bool ParseBool(const std::string& raw_value, const std::string& field_name) {
  const std::string value = ToLowerCopy(raw_value);
  if (value == "1" || value == "true" || value == "yes" || value == "on") {
    return true;
  }
  if (value == "0" || value == "false" || value == "no" || value == "off") {
    return false;
  }
  throw std::runtime_error(
    "Invalid boolean value for '" + field_name + "': '" + raw_value + "'"
  );
}

template <typename T>
T ParseNumber(const std::string& raw_value, const std::string& field_name) {
  std::stringstream stream(raw_value);
  T value {};
  stream >> value;
  if (stream.fail() || !stream.eof()) {
    throw std::runtime_error(
      "Invalid numeric value for '" + field_name + "': '" + raw_value + "'"
    );
  }
  return value;
}

std::optional<std::string> GetEnvValue(const char* name) {
  const char* value = std::getenv(name);
  if (value == nullptr || value[0] == '\0') {
    return std::nullopt;
  }
  return std::string(value);
}

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
  if (key == "max_standard_iterations" || key == "max_standard_iters") {
    config.max_standard_iterations = ParseNumber<int>(raw_value, key);
    return;
  }
  if (key == "max_optimality_condition_iterations" ||
      key == "max_optimality_iters" ||
      key == "max_optcond_iters") {
    config.max_optimality_condition_iterations = ParseNumber<int>(raw_value, key);
    return;
  }
  if (key == "optimization_tolerance" || key == "tolerance") {
    config.optimization_tolerance = ParseNumber<long double>(raw_value, key);
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
  if (key == "nlopt_algorithm" || key == "algorithm") {
    config.nlopt_algorithm = raw_value;
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

  throw std::runtime_error("Unknown metrics option: '" + raw_key + "'");
}

CliOptions ParseCliOptions(int argc, char** argv) {
  CliOptions cli;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      cli.help_requested = true;
      continue;
    }
    if (arg.rfind("--", 0) != 0) {
      throw std::runtime_error("Unsupported argument format: '" + arg + "'");
    }

    arg = arg.substr(2);
    std::string key;
    std::string value;
    const auto eq_pos = arg.find('=');
    if (eq_pos != std::string::npos) {
      key = arg.substr(0, eq_pos);
      value = arg.substr(eq_pos + 1);
    } else {
      key = arg;
      if ((i + 1) < argc && std::string(argv[i + 1]).rfind("--", 0) != 0) {
        value = argv[++i];
      } else {
        value = "true";
      }
    }

    cli.values[NormalizeKey(key)] = value;
  }
  return cli;
}

fs::path FindProjectRoot() {
  fs::path current = fs::current_path();
  while (true) {
    if (fs::exists(current / "CMakeLists.txt") &&
        fs::exists(current / "data" / "TransportationNetworks")) {
      return current;
    }
    if (!current.has_parent_path() || current == current.parent_path()) {
      return fs::current_path();
    }
    current = current.parent_path();
  }
}

fs::path ResolvePath(const fs::path& raw_path, const fs::path& project_root) {
  if (raw_path.empty()) {
    return raw_path;
  }
  if (raw_path.is_absolute()) {
    return raw_path;
  }
  if (fs::exists(project_root / raw_path)) {
    return project_root / raw_path;
  }
  return fs::absolute(raw_path);
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
  int line_number = 0;
  while (std::getline(file, line)) {
    ++line_number;
    const std::string trimmed = TrimCopy(line);
    if (trimmed.empty() || trimmed[0] == '#' || trimmed[0] == ';') {
      continue;
    }

    if (trimmed.front() == '[' && trimmed.back() == ']') {
      current_section = ToLowerCopy(TrimCopy(trimmed.substr(1, trimmed.size() - 2)));
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
  apply_run("CND_MAX_STANDARD_ITERS", "max_standard_iterations");
  apply_run("CND_MAX_OPTCOND_ITERS", "max_optimality_condition_iterations");
  apply_run("CND_TOLERANCE", "optimization_tolerance");
  apply_run("CND_LINK_THRESHOLD", "link_capacity_selection_threshold");
  apply_run("CND_BUDGET_THRESHOLD", "budget_threshold");
  apply_run("CND_BUDGET_MULTIPLIER", "budget_function_multiplier");
  apply_run("CND_BUDGET_UPPER_BOUND", "budget_upper_bound");
  apply_run("CND_NLOPT_ALGORITHM", "nlopt_algorithm");
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

nlopt::algorithm ParseNloptAlgorithm(const std::string& raw_value) {
  const std::string value = ToUpperCopy(raw_value);
  if (value == "LN_COBYLA") return nlopt::LN_COBYLA;
  if (value == "LN_BOBYQA") return nlopt::LN_BOBYQA;
  if (value == "LN_NELDERMEAD") return nlopt::LN_NELDERMEAD;
  if (value == "LN_SBPLX") return nlopt::LN_SBPLX;
  if (value == "LN_PRAXIS") return nlopt::LN_PRAXIS;
  if (value == "LN_NEWUOA") return nlopt::LN_NEWUOA;
  if (value == "LN_NEWUOA_BOUND") return nlopt::LN_NEWUOA_BOUND;
  if (value == "GN_ISRES") return nlopt::GN_ISRES;
  if (value == "GN_AGS") return nlopt::GN_AGS;
  if (value == "GN_ESCH") return nlopt::GN_ESCH;
  if (value == "GN_CRS2_LM") return nlopt::GN_CRS2_LM;
  if (value == "AUGLAG") return nlopt::AUGLAG;
  if (value == "AUGLAG_EQ") return nlopt::AUGLAG_EQ;
  throw std::runtime_error(
    "Unsupported nlopt algorithm '" + raw_value +
    "'. Supported: LN_COBYLA, LN_BOBYQA, LN_NELDERMEAD, LN_SBPLX, LN_PRAXIS, "
    "LN_NEWUOA, LN_NEWUOA_BOUND, GN_ISRES, GN_AGS, GN_ESCH, GN_CRS2_LM, AUGLAG, AUGLAG_EQ."
  );
}

std::shared_ptr<TrafficAssignment::TrafficAssignmentApproach<long double>>
CreateApproach(const RunConfig& config, TrafficAssignment::Network<long double>& network) {
  const std::string approach = ToLowerCopy(config.approach);
  if (approach == "routebased" || approach == "route_based") {
    return std::make_shared<TrafficAssignment::RouteBasedApproach<long double>>(
      network,
      config.approach_alpha,
      config.shift_method,
      config.route_search_threads
    );
  }
  if (approach == "tapas") {
    return std::make_shared<TrafficAssignment::TapasApproach<long double>>(
      network,
      config.approach_alpha,
      config.shift_method
    );
  }
  if (approach == "demandbased" || approach == "demand_based" || approach == "pumpout") {
    return std::make_shared<TrafficAssignment::PumpOutDemandBasedApproach<long double>>(
      network,
      config.approach_alpha
    );
  }

  throw std::runtime_error(
    "Unsupported approach '" + config.approach +
    "'. Supported: RouteBased, Tapas, DemandBased."
  );
}

void PrintHelp() {
  std::cout
    << "Usage: main.exe [options]\n\n"
    << "Layering: defaults -> config file -> environment -> CLI\n\n"
    << "General options:\n"
    << "  --config <path>                  INI config file path\n"
    << "  --dataset <name>                 Dataset name (e.g. SiouxFalls)\n"
    << "  --constraints-file <path>        Constraints CSV path\n"
    << "  --approach <RouteBased|Tapas|DemandBased>\n"
    << "  --shift-method <name>            Shift method (RouteBased/Tapas)\n"
    << "  --route-threads <n>              Route search thread count\n"
    << "  --max-standard-iters <n>\n"
    << "  --max-optcond-iters <n>\n"
    << "  --tolerance <value>\n"
    << "  --link-threshold <value>\n"
    << "  --budget-threshold <value>\n"
    << "  --budget-multiplier <value>\n"
    << "  --budget-upper-bound <value>\n"
    << "  --nlopt-algorithm <LN_COBYLA|LN_BOBYQA|LN_NELDERMEAD|LN_SBPLX|LN_PRAXIS|LN_NEWUOA|LN_NEWUOA_BOUND|GN_ISRES|GN_AGS|GN_ESCH|GN_CRS2_LM|AUGLAG|AUGLAG_EQ>\n"
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
    << "  --metrics_run_id <id>\n\n"
    << "Environment variables (same precedence layer before CLI):\n"
    << "  CND_CONFIG, CND_DATASET, CND_CONSTRAINTS_FILE, CND_APPROACH,\n"
    << "  CND_SHIFT_METHOD, CND_ROUTE_THREADS, CND_MAX_STANDARD_ITERS,\n"
    << "  CND_MAX_OPTCOND_ITERS, CND_TOLERANCE, CND_LINK_THRESHOLD,\n"
    << "  CND_BUDGET_THRESHOLD, CND_BUDGET_MULTIPLIER, CND_BUDGET_UPPER_BOUND,\n"
    << "  CND_NLOPT_ALGORITHM, CND_LOADER_VERBOSE, CND_PRINT_CONFIG,\n"
    << "  CND_METRICS_ENABLE_TRACE, CND_METRICS_ENABLE_RELATIVE_GAP,\n"
    << "  CND_METRICS_RELATIVE_GAP_PERIOD, CND_METRICS_FLUSH_EVERY,\n"
    << "  CND_METRICS_WRITE_METADATA, CND_METRICS_WRITE_SUMMARY,\n"
    << "  CND_METRICS_OUTPUT_ROOT, CND_METRICS_RUN_ID\n";
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
  std::cout << "  max_standard_iterations: " << config.max_standard_iterations << '\n';
  std::cout << "  max_optimality_condition_iterations: "
            << config.max_optimality_condition_iterations << '\n';
  std::cout << "  optimization_tolerance: " << config.optimization_tolerance << '\n';
  std::cout << "  budget_upper_bound: " << config.budget_upper_bound << '\n';
  std::cout << "  nlopt_algorithm: " << config.nlopt_algorithm << '\n';
  std::cout << "  metrics.output_root: " << config.metrics.output_root << '\n';
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

    const nlopt::algorithm algorithm = ParseNloptAlgorithm(config.nlopt_algorithm);
    TrafficAssignment::BilevelCND<long double> cnd(
      network,
      approach,
      constraints,
      config.max_standard_iterations,
      config.max_optimality_condition_iterations,
      config.optimization_tolerance,
      config.link_capacity_selection_threshold,
      config.budget_threshold,
      config.budget_function_multiplier,
      config.budget_upper_bound,
      algorithm,
      config.metrics,
      config.route_search_threads
    );

    std::cout << "Creating BilevelCND solver..." << std::endl;
    std::cout << "      BilevelCND solver created" << std::endl;
    std::cout << "      Design variables: " << constraints.size() << std::endl;

    std::cout << "\nRunning bilevel optimization..." << std::endl;
    std::cout << "      (Requires complete Network and TA implementation)" << std::endl;

    cnd.ComputeNetworkDesign();
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "Fatal error: " << ex.what() << std::endl;
  } catch (...) {
    std::cerr << "Fatal error: unknown exception" << std::endl;
  }

  return 1;
}

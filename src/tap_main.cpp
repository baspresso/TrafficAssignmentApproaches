#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

#include "../include/common/ConfigUtils.h"
#include "../include/tap/algorithms/route_based/RouteBasedApproach.h"
#include "../include/tap/algorithms/tapas/TapasApproach.h"
#include "../include/tap/core/NetworkBuilder.h"

namespace fs = std::filesystem;

using namespace TrafficAssignment::ConfigUtils;

namespace {

struct TapRunConfig {
  std::string dataset = "SiouxFalls";
  std::string approach = "Tapas";
  long double approach_alpha = 1e-14L;
  std::string shift_method = "NewtonStep";
  std::size_t route_search_threads = 1;
  std::string output_root = "performance_results";
  bool print_effective_config = true;

  // TAP approach tuning (0 = use approach default)
  int max_standard_iterations = 0;
  int full_iteration_count = 0;        // RouteBased-specific
  int origin_iteration_count = 0;      // RouteBased-specific
  long double ema_alpha = 0.0L;        // RouteBased-specific
  int pas_per_origin = 0;              // Tapas-specific
  int pas_multiplier = 0;              // Tapas-specific
  int rgap_check_interval = 0;         // Tapas-specific (0 = use default)
  long double tapas_mu = 0.0L;         // Tapas PAS cost-effectiveness (0 = use default 0.5)
  long double tapas_v = 0.0L;          // Tapas PAS flow-effectiveness (0 = use default 0.25)
};

void ApplyTapRunOption(TapRunConfig& config,
                       const std::string& raw_key,
                       const std::string& raw_value) {
  const std::string key = NormalizeKey(raw_key);
  if (key == "dataset") {
    config.dataset = raw_value;
  } else if (key == "approach" || key == "approach_type") {
    config.approach = raw_value;
  } else if (key == "approach_alpha" || key == "alpha") {
    config.approach_alpha = ParseNumber<long double>(raw_value, key);
  } else if (key == "shift_method") {
    config.shift_method = raw_value;
  } else if (key == "route_search_threads" || key == "route_threads") {
    config.route_search_threads = ParseNumber<std::size_t>(raw_value, key);
  } else if (key == "max_standard_iterations" || key == "max_iters") {
    config.max_standard_iterations = ParseNumber<int>(raw_value, key);
  } else if (key == "full_iteration_count") {
    config.full_iteration_count = ParseNumber<int>(raw_value, key);
  } else if (key == "origin_iteration_count") {
    config.origin_iteration_count = ParseNumber<int>(raw_value, key);
  } else if (key == "ema_alpha") {
    config.ema_alpha = ParseNumber<long double>(raw_value, key);
  } else if (key == "pas_per_origin") {
    config.pas_per_origin = ParseNumber<int>(raw_value, key);
  } else if (key == "pas_multiplier") {
    config.pas_multiplier = ParseNumber<int>(raw_value, key);
  } else if (key == "rgap_check_interval") {
    config.rgap_check_interval = ParseNumber<int>(raw_value, key);
  } else if (key == "tapas_mu" || key == "mu") {
    config.tapas_mu = ParseNumber<long double>(raw_value, key);
  } else if (key == "tapas_v" || key == "v") {
    config.tapas_v = ParseNumber<long double>(raw_value, key);
  } else if (key == "output_root") {
    config.output_root = raw_value;
  } else if (key == "print_effective_config" || key == "print_config") {
    config.print_effective_config = ParseBool(raw_value, key);
  } else {
    throw std::runtime_error("Unknown option: '" + raw_key + "'");
  }
}

void ApplyConfigFile(const fs::path& config_path, TapRunConfig& config) {
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
      ApplyTapRunOption(config, key, value);
    } else {
      throw std::runtime_error(
        "Unsupported config section '" + current_section + "' at line " +
        std::to_string(line_number) + " in " + config_path.string() +
        ". Only [run] is supported for TAP."
      );
    }
  }
}

void ApplyEnvironmentOverrides(TapRunConfig& config) {
  auto apply = [&config](const char* env_name, const std::string& key) {
    if (const auto env_value = GetEnvValue(env_name)) {
      ApplyTapRunOption(config, key, *env_value);
    }
  };

  apply("CND_DATASET", "dataset");
  apply("CND_APPROACH", "approach");
  apply("CND_APPROACH_ALPHA", "approach_alpha");
  apply("CND_SHIFT_METHOD", "shift_method");
  apply("CND_ROUTE_THREADS", "route_search_threads");
  apply("CND_MAX_STANDARD_ITERATIONS", "max_standard_iterations");
  apply("CND_FULL_ITERATION_COUNT", "full_iteration_count");
  apply("CND_ORIGIN_ITERATION_COUNT", "origin_iteration_count");
  apply("CND_EMA_ALPHA", "ema_alpha");
  apply("CND_PAS_PER_ORIGIN", "pas_per_origin");
  apply("CND_PAS_MULTIPLIER", "pas_multiplier");
  apply("CND_RGAP_CHECK_INTERVAL", "rgap_check_interval");
  apply("CND_TAPAS_MU", "tapas_mu");
  apply("CND_TAPAS_V", "tapas_v");
  apply("CND_OUTPUT_ROOT", "output_root");
  apply("CND_PRINT_CONFIG", "print_effective_config");
}

void ApplyCliOverrides(const CliOptions& cli, TapRunConfig& config) {
  for (const auto& [raw_key, value] : cli.values) {
    if (raw_key == "config" || raw_key == "help" || raw_key == "h") {
      continue;
    }
    ApplyTapRunOption(config, raw_key, value);
  }
}

std::shared_ptr<TrafficAssignment::TrafficAssignmentApproach<long double>>
CreateApproach(const TapRunConfig& config, TrafficAssignment::Network<long double>& network) {
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
  if (approach == "tapas" || approach == "tasktapas" || approach == "task_tapas" || approach == "task") {
    return std::make_shared<TrafficAssignment::TapasApproach<long double>>(
      network,
      config.approach_alpha,
      config.max_standard_iterations > 0 ? config.max_standard_iterations : 200,
      config.tapas_mu > 0.0L ? config.tapas_mu : 0.5L,
      config.tapas_v > 0.0L ? config.tapas_v : 0.25L
    );
  }
  throw std::runtime_error(
    "Unsupported approach '" + config.approach +
    "'. Supported: RouteBased, Tapas."
  );
}

void PrintHelp() {
  std::cout
    << "Usage: tap_runner.exe [options]\n\n"
    << "Standalone Traffic Assignment solver.\n"
    << "All options have defaults, so no arguments are required.\n\n"
    << "Layering: defaults -> config file -> environment -> CLI\n\n"
    << "Options:\n"
    << "  --config <path>                  INI config file path\n"
    << "  --dataset <name>                 Dataset name (default: SiouxFalls)\n"
    << "  --approach <name>                RouteBased | Tapas (default: Tapas)\n"
    << "  --shift-method <name>            Shift method (default: NewtonStep)\n"
    << "  --approach-alpha <value>         Convergence threshold (default: 1e-14)\n"
    << "  --route-threads <n>              Route search thread count (default: 1)\n"
    << "  --max-standard-iterations <n>    TAP iteration limit (default: 200/20)\n"
    << "  --full-iteration-count <n>       RouteBased OD queue passes (default: 3)\n"
    << "  --origin-iteration-count <n>     RouteBased per-origin passes (default: 1)\n"
    << "  --ema-alpha <value>              RouteBased EMA smoothing (default: 0.7)\n"
    << "  --tapas-mu <value>               Tapas PAS cost-effectiveness (default: 0.5)\n"
    << "  --tapas-v <value>                Tapas PAS flow-effectiveness (default: 0.25)\n"
    << "  --output-root <path>             Output directory (default: performance_results)\n"
    << "  --print-config <true|false>      Print effective config (default: true)\n"
    << "  --help                           Show this help\n\n"
    << "Environment variables (CND_ prefix):\n"
    << "  CND_DATASET, CND_APPROACH, CND_APPROACH_ALPHA, CND_SHIFT_METHOD,\n"
    << "  CND_ROUTE_THREADS, CND_MAX_STANDARD_ITERATIONS,\n"
    << "  CND_FULL_ITERATION_COUNT, CND_ORIGIN_ITERATION_COUNT, CND_EMA_ALPHA,\n"
    << "  CND_PAS_PER_ORIGIN, CND_PAS_MULTIPLIER, CND_RGAP_CHECK_INTERVAL,\n"
    << "  CND_TAPAS_MU, CND_TAPAS_V,\n"
    << "  CND_OUTPUT_ROOT, CND_PRINT_CONFIG\n\n"
    << "Examples:\n"
    << "  tap_runner.exe --dataset SiouxFalls --approach Tapas\n"
    << "  tap_runner.exe --config configs/tap.siouxfalls.ini\n"
    << "  tap_runner.exe --dataset Anaheim --approach RouteBased --shift-method NewtonStep\n"
    << "\nSupported approaches: RouteBased, Tapas\n";
}

void PrintEffectiveConfig(const TapRunConfig& config,
                          const fs::path& project_root,
                          const std::optional<fs::path>& config_path) {
  std::cout << "\nEffective TAP configuration:\n";
  std::cout << "  project_root: " << project_root.string() << '\n';
  std::cout << "  config_file: " << (config_path ? config_path->string() : "<none>") << '\n';
  std::cout << "  dataset: " << config.dataset << '\n';
  std::cout << "  approach: " << config.approach << '\n';
  std::cout << "  shift_method: " << config.shift_method << '\n';
  std::cout << "  approach_alpha: " << config.approach_alpha << '\n';
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
  if (config.tapas_mu > 0.0L)
    std::cout << "  tapas_mu: " << config.tapas_mu << '\n';
  if (config.tapas_v > 0.0L)
    std::cout << "  tapas_v: " << config.tapas_v << '\n';
  std::cout << "  output_root: " << config.output_root << '\n';
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
    TapRunConfig config;

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

    if (config.print_effective_config) {
      PrintEffectiveConfig(config, project_root, config_path);
    }

    std::cout << "Building network from dataset '" << config.dataset << "'..." << std::endl;
    TrafficAssignment::NetworkBuilder builder;
    auto network = builder.BuildFromDataset<long double>(config.dataset);

    std::cout << "Creating approach: " << config.approach << " (" << config.shift_method << ")" << std::endl;
    auto approach = CreateApproach(config, network);

    // Configure output path for statistics
    const fs::path output_root = ResolvePath(fs::path(config.output_root), project_root);
    approach->SetOutputRoot(output_root.string());

    std::cout << "\nRunning Traffic Assignment...\n" << std::endl;
    approach->ComputeTrafficFlows(true);

    // Print final summary
    std::cout << std::setprecision(15);
    std::cout << "\n=== Traffic Assignment Results ===" << std::endl;
    std::cout << "  Approach:          " << approach->GetApproachName() << std::endl;
    std::cout << "  Dataset:           " << config.dataset << std::endl;
    std::cout << "  RelativeGap:       " << network.RelativeGap() << std::endl;
    std::cout << "  ObjectiveFunction: " << network.ObjectiveFunction() << std::endl;
    std::cout << "  TotalTravelTime:   " << network.TotalTravelTime() << std::endl;

    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "Fatal error: " << ex.what() << std::endl;
  } catch (...) {
    std::cerr << "Fatal error: unknown exception" << std::endl;
  }

  return 1;
}

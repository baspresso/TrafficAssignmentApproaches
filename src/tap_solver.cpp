#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

#ifdef _WIN32
  #include <io.h>
  #include <cstdio>
#else
  #include <unistd.h>
#endif

#include "../include/common/TomlConfigLoader.h"
#include "../include/tap/algorithms/route_based/RouteBasedApproach.h"
#include "../include/tap/algorithms/tapas/TapasApproach.h"
#include "../include/tap/core/NetworkBuilder.h"

namespace fs = std::filesystem;

using namespace TrafficAssignment::ConfigUtils;
using namespace TrafficAssignment::Config;

namespace {

std::shared_ptr<TrafficAssignment::TrafficAssignmentApproach<long double>>
CreateApproach(const TapConfig& config, TrafficAssignment::Network<long double>& network) {
  const std::string approach = ToLowerCopy(config.solver.approach);
  if (approach == "routebased" || approach == "route_based") {
    return std::make_shared<TrafficAssignment::RouteBasedApproach<long double>>(
      network,
      config.solver.approach_alpha,
      config.solver.route_based.shift_method,
      config.solver.route_based.route_search_threads,
      config.solver.max_standard_iterations > 0 ? config.solver.max_standard_iterations : 200,
      config.solver.route_based.full_iteration_count > 0 ? config.solver.route_based.full_iteration_count : 3,
      config.solver.route_based.origin_iteration_count > 0 ? config.solver.route_based.origin_iteration_count : 1,
      config.solver.route_based.ema_alpha > 0.0L ? config.solver.route_based.ema_alpha : 0.7L
    );
  }
  if (approach == "tapas" || approach == "tasktapas" || approach == "task_tapas" || approach == "task") {
    return std::make_shared<TrafficAssignment::TapasApproach<long double>>(
      network,
      config.solver.approach_alpha,
      config.solver.max_standard_iterations > 0 ? config.solver.max_standard_iterations : 200,
      config.solver.tapas.mu > 0.0L ? config.solver.tapas.mu : 0.5L,
      config.solver.tapas.v > 0.0L ? config.solver.tapas.v : 0.25L
    );
  }
  throw std::runtime_error(
    "Unsupported approach '" + config.solver.approach +
    "'. Supported: RouteBased, Tapas."
  );
}

void PrintHelp() {
  std::cout
    << "Usage: tap_solver.exe [options]\n\n"
    << "Standalone Traffic Assignment solver.\n"
    << "All options have defaults, so no arguments are required.\n\n"
    << "Layering: defaults -> TOML config file -> environment -> CLI\n\n"
    << "Options:\n"
    << "  --config <path>                  TOML config file path\n"
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
    << "  --quiet <auto|true|false>        Quiet mode (auto=quiet when piped, default: auto)\n"
    << "  --help                           Show this help\n\n"
    << "Environment variables (CND_ prefix):\n"
    << "  CND_DATASET, CND_APPROACH, CND_APPROACH_ALPHA, CND_SHIFT_METHOD,\n"
    << "  CND_ROUTE_THREADS, CND_MAX_STANDARD_ITERATIONS,\n"
    << "  CND_FULL_ITERATION_COUNT, CND_ORIGIN_ITERATION_COUNT, CND_EMA_ALPHA,\n"
    << "  CND_PAS_PER_ORIGIN, CND_PAS_MULTIPLIER, CND_RGAP_CHECK_INTERVAL,\n"
    << "  CND_TAPAS_MU, CND_TAPAS_V,\n"
    << "  CND_OUTPUT_ROOT, CND_PRINT_CONFIG\n\n"
    << "Examples:\n"
    << "  tap_solver.exe --dataset SiouxFalls --approach Tapas\n"
    << "  tap_solver.exe --config configs/tap.siouxfalls.toml\n"
    << "  tap_solver.exe --dataset Anaheim --approach RouteBased --shift-method NewtonStep\n"
    << "\nSupported approaches: RouteBased, Tapas\n";
}

void PrintEffectiveConfig(const TapConfig& config,
                          const fs::path& project_root,
                          const std::optional<fs::path>& config_path) {
  std::cout << "\nEffective TAP configuration:\n";
  std::cout << "  project_root: " << project_root.string() << '\n';
  std::cout << "  config_file: " << (config_path ? config_path->string() : "<none>") << '\n';
  std::cout << "  dataset: " << config.network.dataset << '\n';
  std::cout << "  approach: " << config.solver.approach << '\n';
  std::cout << "  shift_method: " << config.solver.route_based.shift_method << '\n';
  std::cout << "  approach_alpha: " << config.solver.approach_alpha << '\n';
  std::cout << "  route_search_threads: " << config.solver.route_based.route_search_threads << '\n';
  if (config.solver.max_standard_iterations > 0)
    std::cout << "  max_standard_iterations: " << config.solver.max_standard_iterations << '\n';
  if (config.solver.route_based.full_iteration_count > 0)
    std::cout << "  full_iteration_count: " << config.solver.route_based.full_iteration_count << '\n';
  if (config.solver.route_based.origin_iteration_count > 0)
    std::cout << "  origin_iteration_count: " << config.solver.route_based.origin_iteration_count << '\n';
  if (config.solver.route_based.ema_alpha > 0.0L)
    std::cout << "  ema_alpha: " << config.solver.route_based.ema_alpha << '\n';
  if (config.solver.tapas.pas_per_origin > 0)
    std::cout << "  pas_per_origin: " << config.solver.tapas.pas_per_origin << '\n';
  if (config.solver.tapas.pas_multiplier > 0)
    std::cout << "  pas_multiplier: " << config.solver.tapas.pas_multiplier << '\n';
  if (config.solver.tapas.rgap_check_interval > 0)
    std::cout << "  rgap_check_interval: " << config.solver.tapas.rgap_check_interval << '\n';
  if (config.solver.tapas.mu > 0.0L)
    std::cout << "  tapas_mu: " << config.solver.tapas.mu << '\n';
  if (config.solver.tapas.v > 0.0L)
    std::cout << "  tapas_v: " << config.solver.tapas.v << '\n';
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
    TapConfig config;

    std::optional<fs::path> config_path;
    const auto cli_config_it = cli.values.find("config");
    if (cli_config_it != cli.values.end()) {
      config_path = ResolvePath(fs::path(cli_config_it->second), project_root);
    } else if (const auto env_config = GetEnvValue("CND_CONFIG")) {
      config_path = ResolvePath(fs::path(*env_config), project_root);
    }

    if (config_path.has_value()) {
      config = LoadTapConfig(*config_path);
    }

    ApplyEnvironmentOverrides(config);
    ApplyCliOverrides(cli, config);

    // Resolve quiet mode: auto = quiet when stdout is piped
    bool quiet;
    if (config.output.quiet == "true") {
      quiet = true;
    } else if (config.output.quiet == "false") {
      quiet = false;
    } else {
      #ifdef _WIN32
        quiet = !_isatty(_fileno(stdout));
      #else
        quiet = !isatty(STDOUT_FILENO);
      #endif
    }

    if (config.output.print_effective_config && !quiet) {
      PrintEffectiveConfig(config, project_root, config_path);
    }

    if (!quiet) {
      std::cout << "Building network from dataset '" << config.network.dataset << "'..." << std::endl;
    }
    TrafficAssignment::NetworkBuilder builder;
    auto network = builder.BuildFromDataset<long double>(config.network.dataset);

    if (!quiet) {
      std::cout << "Creating approach: " << config.solver.approach << " (" << config.solver.route_based.shift_method << ")" << std::endl;
    }
    auto approach = CreateApproach(config, network);

    // Configure output path for statistics
    const fs::path output_root = ResolvePath(fs::path(config.output_root), project_root);
    approach->SetOutputRoot(output_root.string());

    if (!quiet) {
      std::cout << "\nRunning Traffic Assignment...\n" << std::endl;
    }
    approach->ComputeTrafficFlows(true);

    // Always emit structured [RESULT] line
    std::cout << "[RESULT]"
              << " approach=" << approach->GetApproachName()
              << " dataset=" << config.network.dataset
              << " relative_gap=" << std::setprecision(15) << network.RelativeGap()
              << " objective_function=" << network.ObjectiveFunction()
              << " total_travel_time=" << network.TotalTravelTime()
              << std::endl;

    // Print human-readable results when not quiet
    if (!quiet) {
      std::cout << std::setprecision(15);
      std::cout << "\n=== Traffic Assignment Results ===" << std::endl;
      std::cout << "  Approach:          " << approach->GetApproachName() << std::endl;
      std::cout << "  Dataset:           " << config.network.dataset << std::endl;
      std::cout << "  RelativeGap:       " << network.RelativeGap() << std::endl;
      std::cout << "  ObjectiveFunction: " << network.ObjectiveFunction() << std::endl;
      std::cout << "  TotalTravelTime:   " << network.TotalTravelTime() << std::endl;
    }

    // Write link flow distribution
    {
      const auto flows_dir = output_root / config.network.dataset;
      std::filesystem::create_directories(flows_dir);
      const auto flows_path = flows_dir / "link_flows.csv";
      std::ofstream flows_file(flows_path, std::ios::out);
      flows_file << "link_id,init_node,term_node,capacity,free_flow_time,flow,cost\n";
      flows_file << std::setprecision(10);
      const auto& links = network.links();
      for (int i = 0; i < network.number_of_links(); ++i) {
        const auto& link = links[i];
        flows_file << i << ","
                   << link.init << ","
                   << link.term << ","
                   << link.capacity << ","
                   << link.free_flow_time << ","
                   << link.flow << ","
                   << link.Delay() << "\n";
      }
      if (!quiet) {
        std::cout << "  Link flows:        " << flows_path.string() << std::endl;
      }
    }

    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "Fatal error: " << ex.what() << std::endl;
  } catch (...) {
    std::cerr << "Fatal error: unknown exception" << std::endl;
  }

  return 1;
}

#include <filesystem>
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
#include "../include/cnd/BilevelCND.h"
#include "../include/cnd/DirectedConstraintLoader.h"
#include "../include/tap/algorithms/route_based/RouteBasedApproach.h"
#include "../include/tap/algorithms/tapas/TapasApproach.h"
#include "../include/tap/core/NetworkBuilder.h"

namespace fs = std::filesystem;

using namespace TrafficAssignment::ConfigUtils;
using namespace TrafficAssignment::Config;

namespace {

fs::path ResolveConstraintsPath(const CndpConfig& config, const fs::path& project_root) {
  if (config.network.constraints_file.empty()) {
    return project_root / "data" / "TransportationNetworks" / config.network.dataset /
           (config.network.dataset + "_constraints.csv");
  }
  return ResolvePath(fs::path(config.network.constraints_file), project_root);
}

std::shared_ptr<TrafficAssignment::TrafficAssignmentApproach<long double>>
CreateApproach(const CndpConfig& config, TrafficAssignment::Network<long double>& network) {
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
    << "Usage: cndp_solver.exe --config <path> [options]\n\n"
    << "Layering: defaults -> TOML config file -> environment -> CLI\n\n"
    << "Optimization steps are defined via [[pipeline]] in the TOML config file.\n"
    << "Example:\n"
    << "  [[pipeline]]\n"
    << "  type = \"nlopt\"\n"
    << "  algorithm = \"LN_COBYLA\"\n"
    << "  max_iterations = 200\n"
    << "  tolerance = 1e-4\n\n"
    << "  [[pipeline]]\n"
    << "  type = \"optimality_condition\"\n"
    << "  max_iterations = 30\n\n"
    << "Step options:\n"
    << "  type                             nlopt | optimality_condition | optimlib | gradient_descent\n"
    << "  algorithm                        NLopt: LN_COBYLA, GN_ISRES, ...\n"
    << "                                   OptimLib: DE, PSO, NM, DE_PRMM, PSO_DV, GD\n"
    << "  local_algorithm                  Local sub-algorithm for AUGLAG (NLopt)\n"
    << "  max_iterations                   Max iterations/generations for this step\n"
    << "  tolerance                        Convergence tolerance\n"
    << "  population_size                  Population size for OptimLib (0=auto)\n"
    << "  step_size                        Initial step size for gradient descent (default: 1.0)\n"
    << "  fd_epsilon                       Finite difference epsilon for gradient descent (default: 1e-4)\n"
    << "  gradient_method                  Gradient estimator: finite_difference (default), spsa, sensitivity\n"
    << "  stochastic_optimizer             Stochastic optimizer: sgd (default), momentum, adam\n"
    << "  name                             Display name for this step\n\n"
    << "General options (CLI or TOML sections):\n"
    << "  --config <path>                  TOML config file path (required)\n"
    << "  --dataset <name>                 Dataset name (e.g. SiouxFalls)\n"
    << "  --constraints-file <path>        Constraints CSV path\n"
    << "  --approach <RouteBased|Tapas>\n"
    << "  --shift-method <name>            Shift method (RouteBased)\n"
    << "  --route-threads <n>              Route search thread count\n"
    << "  --max-standard-iterations <n>    TAP iteration limit (default: 200/20)\n"
    << "  --full-iteration-count <n>       RouteBased OD queue passes (default: 3)\n"
    << "  --origin-iteration-count <n>     RouteBased per-origin passes (default: 1)\n"
    << "  --ema-alpha <value>              RouteBased EMA smoothing (default: 0.7)\n"
    << "  --tapas-mu <value>               Tapas PAS cost-effectiveness (default: 0.5)\n"
    << "  --tapas-v <value>                Tapas PAS flow-effectiveness (default: 0.25)\n"
    << "  --link-threshold <value>\n"
    << "  --budget-threshold <value>\n"
    << "  --budget-multiplier <value>\n"
    << "  --budget-upper-bound <value>\n"
    << "  --progress-format <auto|bar|line|none>  Progress output format (default: auto)\n"
    << "  --loader-verbose <true|false>\n"
    << "  --print-config <true|false>\n"
    << "  --quiet <auto|true|false>        Quiet mode (auto=quiet when piped, default: auto)\n\n"
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

void PrintEffectiveConfig(const CndpConfig& config,
                          const fs::path& project_root,
                          const fs::path& constraints_path,
                          const std::optional<fs::path>& config_path) {
  std::cout << "\nEffective run configuration:\n";
  std::cout << "  project_root: " << project_root.string() << '\n';
  std::cout << "  config_file: " << (config_path ? config_path->string() : "<none>") << '\n';
  std::cout << "  dataset: " << config.network.dataset << '\n';
  std::cout << "  constraints_file: " << constraints_path.string() << '\n';
  std::cout << "  approach: " << config.solver.approach << '\n';
  std::cout << "  shift_method: " << config.solver.route_based.shift_method << '\n';
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
  std::cout << "  budget_upper_bound: " << config.solver.budget_upper_bound << '\n';
  std::cout << "  metrics.output_root: " << config.metrics.output_root << '\n';
  if (!config.metrics.append_dataset_subdir) {
    std::cout << "  metrics.no_dataset_subdir: true\n";
  }
  if (!config.metrics.scenario_name.empty()) {
    std::cout << "  metrics.scenario_name: " << config.metrics.scenario_name << '\n';
  }
  std::cout << "  pipeline steps: " << config.pipeline.size() << '\n';
  for (std::size_t i = 0; i < config.pipeline.size(); ++i) {
    const auto& step = config.pipeline[i];
    std::cout << "    [[pipeline]] #" << i << " type=" << step.type;
    if (!step.algorithm.empty()) std::cout << " algorithm=" << step.algorithm;
    std::cout << " max_iter=" << step.max_iterations;
    if (step.tolerance > 0) std::cout << " tol=" << step.tolerance;
    if (step.population_size > 0) std::cout << " pop_size=" << step.population_size;
    if (step.type == "gradient_descent") {
      std::cout << " step_size=" << step.step_size << " fd_eps=" << step.fd_epsilon;
      if (!step.gradient_method.empty()) {
        std::cout << " gradient_method=" << step.gradient_method;
      }
      if (!step.stochastic_optimizer.empty()) {
        std::cout << " optimizer=" << step.stochastic_optimizer;
      }
    }
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
    CndpConfig config;

    std::optional<fs::path> config_path;
    const auto cli_config_it = cli.values.find("config");
    if (cli_config_it != cli.values.end()) {
      config_path = ResolvePath(fs::path(cli_config_it->second), project_root);
    } else if (const auto env_config = GetEnvValue("CND_CONFIG")) {
      config_path = ResolvePath(fs::path(*env_config), project_root);
    }

    if (config_path.has_value()) {
      config = LoadCndpConfig(*config_path);
    }

    ApplyEnvironmentOverrides(config);
    ApplyCliOverrides(cli, config);

    // Resolve "auto" progress format: use bar for TTY, line for piped stdout
    if (config.output.progress_format == "auto") {
      #ifdef _WIN32
        config.output.progress_format = _isatty(_fileno(stdout)) ? "bar" : "line";
      #else
        config.output.progress_format = isatty(STDOUT_FILENO) ? "bar" : "line";
      #endif
    }

    // Resolve quiet mode: auto = quiet when stdout is piped
    bool quiet;
    if (config.output.quiet == "true") {
      quiet = true;
    } else if (config.output.quiet == "false") {
      quiet = false;
    } else {
      // auto: quiet when stdout is not a TTY (piped)
      #ifdef _WIN32
        quiet = !_isatty(_fileno(stdout));
      #else
        quiet = !isatty(STDOUT_FILENO);
      #endif
    }

    // When quiet, force line progress format (structured output only)
    if (quiet) {
      config.output.progress_format = "line";
    }

    const fs::path constraints_path = ResolveConstraintsPath(config, project_root);
    if (!fs::exists(constraints_path)) {
      throw std::runtime_error(
        "Constraint file does not exist: " + constraints_path.string()
      );
    }

    if (config.output.print_effective_config && !quiet) {
      PrintEffectiveConfig(config, project_root, constraints_path, config_path);
    }

    TrafficAssignment::NetworkBuilder builder;
    auto network = builder.BuildFromDataset<long double>(config.network.dataset);
    auto approach = CreateApproach(config, network);

    TrafficAssignment::DirectedConstraintLoader loader;
    loader.SetVerbose(config.output.verbose && !quiet);
    auto constraints = loader.LoadFromFile(constraints_path.string());

    // Exclude flow-insensitive links from optimization:
    // - b ~ 0 or power ~ 0: delay independent of flow/capacity
    // - free_flow_time ~ 0: delay always near zero regardless of capacity
    int excluded_count = 0;
    for (int i = 0; i < static_cast<int>(constraints.size()); ++i) {
      const auto& link = network.links()[i];
      if (std::abs(static_cast<double>(link.b)) < 1e-10 ||
          std::abs(static_cast<double>(link.power)) < 1e-10 ||
          std::abs(static_cast<double>(link.free_flow_time)) < 1e-10) {
        constraints[i].upper_bound = constraints[i].lower_bound;
        ++excluded_count;
      }
    }
    int active_count = static_cast<int>(constraints.size()) - excluded_count;
    if (!quiet) {
      std::cout << "      Link filtering: " << excluded_count << " of " << constraints.size()
                << " links excluded (b ~ 0, power ~ 0, or free_flow_time ~ 0), "
                << active_count << " active design variables" << std::endl;
    }

    if (config.pipeline.empty()) {
      throw std::runtime_error(
        "No [[pipeline]] entries found in config. At least one optimization step "
        "must be defined. Example:\n"
        "  [[pipeline]]\n"
        "  type = \"nlopt\"\n"
        "  algorithm = \"LN_COBYLA\"\n"
        "  max_iterations = 100\n"
      );
    }

    if (!quiet) {
      std::cout << "Creating BilevelCND solver..." << std::endl;
    }
    TrafficAssignment::BilevelCND<long double> cnd(
      network,
      approach,
      constraints,
      config.pipeline,
      config.solver.link_capacity_selection_threshold,
      config.solver.budget_threshold,
      config.solver.budget_function_multiplier,
      config.solver.budget_upper_bound,
      config.metrics,
      config.solver.route_based.route_search_threads,
      config.output.progress_format
    );

    cnd.SetVerbose(!quiet);

    if (!quiet) {
      std::cout << "      BilevelCND solver created" << std::endl;
      std::cout << "      Design variables: " << active_count
                << " active (" << excluded_count << " excluded of "
                << constraints.size() << " total)" << std::endl;
    }

    if (!quiet) {
      std::cout << "\nRunning bilevel optimization..." << std::endl;
    }

    cnd.ComputeNetworkDesign();
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "Fatal error: " << ex.what() << std::endl;
  } catch (...) {
    std::cerr << "Fatal error: unknown exception" << std::endl;
  }

  return 1;
}

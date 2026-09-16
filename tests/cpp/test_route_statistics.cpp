#include <traffic_assignment/solve.hpp>

#include "detail/access.hpp"
#include "detail/cnd/CndStatisticsRecorder.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>

namespace ta = traffic_assignment;
namespace fs = std::filesystem;

void Require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

struct TemporaryDirectory {
  fs::path path = fs::temp_directory_path() /
      ("ta_route_statistics_" + std::to_string(
          std::chrono::steady_clock::now().time_since_epoch().count()));
  TemporaryDirectory() { fs::create_directory(path); }
  ~TemporaryDirectory() {
    std::error_code error;
    fs::remove_all(path, error);
  }
};

int main() {
  try {
    TemporaryDirectory temporary;
    ta::NetworkData data;
    data.name = "TwoRoutes";
    data.number_of_nodes = 3;
    data.number_of_zones = 2;
    data.links = {{0, 1, 10, 0, 2, 1, 1}, {0, 2, 10, 0, 1, 1, 1},
                  {2, 1, 10, 0, 1, 1, 1}};
    data.demand = {{0, 10}, {0, 0}};
    auto network = std::make_shared<ta::Network>(data);
    ta::TapOptions options;
    options.approach = "routebased";
    auto route_approach = ta::MakeApproach(network, options);
    auto tapas_approach = ta::MakeApproach(network);
    auto& native_network = ta::detail::Access::Get(*network);
    auto& od_pair = native_network.mutable_od_pairs()[0];
    od_pair.AddNewRoute({0});
    od_pair.AddNewRoute({1, 2});  // Stored with zero flow.
    Require(od_pair.GetRoutesCount() == 2, "fixture has two stored routes");

    const auto flows = network->flows();
    const auto objective = network->beckmann_objective();
    const auto counts = ta::detail::Access::Get(*route_approach)->GetRouteCountsForStatistics();
    Require(counts == std::vector<std::size_t>{1}, "RouteBased excludes zero-flow routes");
    Require(ta::detail::Access::Get(*tapas_approach)->GetRouteCountsForStatistics() ==
                std::vector<std::size_t>{2}, "TAPAS reports raw stored route counts");

    TrafficAssignment::CndStatisticsRecorder<ta::Real> recorder(native_network);
    ta::MetricsConfig metrics;
    metrics.output_root = temporary.path.string();
    metrics.append_dataset_subdir = false;
    metrics.run_id = "final";
    metrics.enable_trace = false;
    metrics.write_summary_csv = false;
    metrics.write_metadata_json = false;
    TrafficAssignment::CndRunMetadata metadata;
    metadata.approach_name = route_approach->GetApproachName();
    TrafficAssignment::CndRunSummary summary;

    recorder.StopRun(summary, counts);
    Require(fs::is_empty(temporary.path), "inactive recorder writes no files");
    recorder.StartRun(metrics, metadata);
    recorder.StopRun(summary, counts);
    const auto path = temporary.path / "BilevelCND_RouteBasedNewtonStep_final_route_counts.csv";
    std::ifstream file(path);
    Require(static_cast<bool>(file), "successful recording writes route counts with other outputs off");
    const std::string contents((std::istreambuf_iterator<char>(file)), {});
    Require(contents == "od_pair_index;init_node;dest_node;routes_count\n0;0;1;1\n",
            "snapshot contains zero-based OD identity and the selected count");
    Require(network->flows() == flows && network->beckmann_objective() == objective,
            "counting and recording leave the assignment unchanged");

    od_pair.SetRoutesFlow({10, 1e-30L});
    Require(ta::detail::Access::Get(*route_approach)->GetRouteCountsForStatistics() ==
                std::vector<std::size_t>{2}, "even tiny positive route flows count");
    // Finalization is idempotent and cannot replace the captured final state.
    recorder.StopRun(summary, {2});
    file.clear();
    file.seekg(0);
    Require(std::string((std::istreambuf_iterator<char>(file)), {}) == contents,
            "stopped recorder keeps its original snapshot");

    metrics.run_id = "failed";
    recorder.StartRun(metrics, metadata);
    summary.status = "failure";
    recorder.StopRun(summary, {});
    Require(!fs::exists(temporary.path / "BilevelCND_RouteBasedNewtonStep_failed_route_counts.csv"),
            "failed runs do not write a final snapshot");
    std::cout << "Route statistics checks passed\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}

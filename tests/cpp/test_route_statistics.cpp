#include <traffic_assignment/solve.hpp>

#include "detail/access.hpp"
#include "detail/cnd/CndStatisticsRecorder.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
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

std::vector<std::vector<std::string>> ReadRows(const fs::path& path, char delimiter) {
  std::ifstream file(path);
  Require(static_cast<bool>(file), "expected statistics file exists");
  std::vector<std::vector<std::string>> rows;
  for (std::string line; std::getline(file, line);) {
    std::istringstream input(line);
    std::vector<std::string> row;
    for (std::string field; std::getline(input, field, delimiter);) row.push_back(field);
    if (!line.empty() && line.back() == delimiter) row.emplace_back();
    rows.push_back(std::move(row));
  }
  return rows;
}

void TestRouteHistory() {
  TemporaryDirectory temporary;
  ta::NetworkData data;
  data.name = "TwoODs";
  data.number_of_nodes = data.number_of_zones = 2;
  data.links = {{0, 1, 10, 0, 2, 1, 1}, {1, 0, 10, 0, 2, 1, 1}};
  data.demand = {{0, 10}, {5, 0}};
  ta::Network network(data);
  TrafficAssignment::CndStatisticsRecorder<ta::Real> recorder(ta::detail::Access::Get(network));
  ta::MetricsConfig metrics;
  metrics.output_root = temporary.path.string();
  metrics.append_dataset_subdir = false;
  metrics.write_summary_csv = metrics.write_metadata_json = false;
  TrafficAssignment::CndRunMetadata metadata;
  metadata.approach_name = "RouteBasedNewtonStep";
  metadata.scenario_name = "History";

  for (int flush_period : {0, 2}) {
    metrics.run_id = "history_" + std::to_string(flush_period);
    metrics.flush_every_n_points = flush_period;
    const auto prefix = "BilevelCND_RouteBasedNewtonStep_" + metrics.run_id;
    const auto history_path = temporary.path / (prefix + "_route_counts_history.csv");
    const auto trace_path = temporary.path / (prefix + "_quality_time.csv");
    recorder.StartRun(metrics, metadata);
    Require(recorder.IsTraceEnabled(), "history enabled with convergence tracing");

    std::vector<std::size_t> counts{1, 2};
    recorder.LogPoint("standard_eval", 1, 20, 15, 5, 100, 0, 0.1, 0, counts);
    Require(ReadRows(history_path, ';').size() == 1, "first snapshot stays buffered");
    counts = {2, 0};
    // Repeated phase/step values still refer to separate trace points.
    recorder.LogPoint("standard_eval", 1, 19, 14, 5, 100, 0, 0.1, 0, counts);
    const std::size_t flushed_rows = flush_period == 0 ? 1 : 5;
    Require(ReadRows(history_path, ';').size() == flushed_rows, "history follows flush interval");
    const auto nan = std::numeric_limits<double>::quiet_NaN();
    recorder.LogPoint("gd_eval", 2, nan, nan, nan, 100, 0, 0.1, nan);
    bool rejected = false;
    try {
      recorder.LogPoint("bad_snapshot", 3, 19, 14, 5, 100, 0, 0.1, 0,
                        std::vector<std::size_t>{1});
    } catch (const std::invalid_argument&) {
      rejected = true;
    }
    Require(rejected, "available snapshot must have one count per OD pair");

    TrafficAssignment::CndRunSummary summary;
    summary.status = "failure";
    recorder.StopRun(summary, {});
    const auto history = ReadRows(history_path, ';');
    const auto trace = ReadRows(trace_path, ',');
    Require(history.size() == 7 && trace.size() == 4, "failure flushes all pending points");
    Require(history[0] == std::vector<std::string>{
      "RunId", "Scenario", "TracePointIndex", "ElapsedTime(s)", "Phase", "Step",
      "od_pair_index", "init_node", "dest_node", "routes_count"}, "history CSV header");
    const std::vector<std::string> expected_counts{"1", "2", "2", "0", "", ""};
    for (std::size_t row = 1; row < history.size(); ++row) {
      const auto point = (row - 1) / 2;
      const auto od = (row - 1) % 2;
      const auto& fields = history[row];
      Require(fields.size() == 10, "history row has all fields, including unavailable counts");
      Require(fields[0] == metrics.run_id && fields[1] == "History", "history identifies run");
      Require(fields[2] == std::to_string(point), "point index follows trace row order");
      Require(fields[3] == trace[point + 1][5] && fields[4] == trace[point + 1][6] &&
              fields[5] == trace[point + 1][7], "history shares trace timestamp, phase and step");
      Require(fields[6] == std::to_string(od) && fields[7] == std::to_string(od) &&
              fields[8] == std::to_string(1 - od), "every point contains stable OD identities");
      Require(fields[9] == expected_counts[row - 1], "counts are immutable snapshots, invalid counts blank");
    }
    Require(!fs::exists(temporary.path / (prefix + "_route_counts.csv")),
            "failed runs retain history without a final snapshot");
    recorder.StopRun(summary, {});
    Require(ReadRows(history_path, ';') == history, "repeated finalization does not duplicate history");

    metrics.run_id = "reused_" + std::to_string(flush_period);
    recorder.StartRun(metrics, metadata);
    recorder.LogPoint("final", 0, 19, 14, 5, 100, 0, 0.1, 0, counts);
    summary.status = "success";
    recorder.StopRun(summary, counts);
    const auto reused_prefix = "BilevelCND_RouteBasedNewtonStep_" + metrics.run_id;
    const auto reused = ReadRows(temporary.path / (reused_prefix + "_route_counts_history.csv"), ';');
    Require(reused.size() == 3 && reused[1][2] == "0" && reused[2][2] == "0",
            "recorder reuse resets point indices and flush state");
    const auto final_counts = ReadRows(temporary.path / (reused_prefix + "_route_counts.csv"), ';');
    Require(reused[1][9] == final_counts[1][3] && reused[2][9] == final_counts[2][3],
            "final history agrees with the preserved final snapshot");
  }
}

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
    recorder.LogPoint("inactive", 0, 20, 15, 5, 100, 0, 0.1, 0, counts);
    Require(fs::is_empty(temporary.path), "inactive recorder writes no files");
    recorder.StartRun(metrics, metadata);
    Require(!recorder.IsTraceEnabled(), "history disabled with convergence tracing");
    recorder.LogPoint("standard_eval", 1, 20, 15, 5, 100, 0, 0.1, 0, counts);
    recorder.StopRun(summary, counts);
    Require(!fs::exists(temporary.path / "BilevelCND_RouteBasedNewtonStep_final_route_counts_history.csv"),
            "disabled tracing writes no route history");
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
    TestRouteHistory();
    std::cout << "Route statistics checks passed\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}

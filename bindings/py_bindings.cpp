// Python bindings for the TAP / CNDP solvers (module traffic_assignment._core).
//
// Everything is instantiated with long double to reproduce the numerics of the
// cndp_solver / tap_solver executables; values cross the Python boundary as
// double. Network objects are constructed directly on the heap and bound as
// non-copyable: OD pairs capture a reference to their owning network, so a
// copied Network would dangle.

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../include/cnd/BilevelCND.h"
#include "../include/cnd/CndStatisticsRecorder.h"
#include "../include/cnd/DirectedConstraintLoader.h"
#include "../include/cnd/OptimizationStep.h"
#include "../include/tap/algorithms/TrafficAssignmentApproach.h"
#include "../include/tap/algorithms/route_based/RouteBasedApproach.h"
#include "../include/tap/algorithms/tapas/TapasApproach.h"
#include "../include/tap/core/Network.h"
#include "../include/tap/core/NetworkBuilder.h"

namespace py = pybind11;

namespace {

using Real = long double;
using TrafficAssignment::BilevelCND;
using TrafficAssignment::CndMetricsConfig;
using TrafficAssignment::DirectedConstraintLoader;
using TrafficAssignment::DirectedLinkCapacityConstraint;
using TrafficAssignment::Link;
using TrafficAssignment::Network;
using TrafficAssignment::NetworkBuilder;
using TrafficAssignment::OptimizationStepConfig;
using TrafficAssignment::RouteBasedApproach;
using TrafficAssignment::TapasApproach;
using TrafficAssignment::TrafficAssignmentApproach;

using DoubleArray = py::array_t<double, py::array::c_style | py::array::forcecast>;
using IndexArray = py::array_t<std::int64_t, py::array::c_style | py::array::forcecast>;

template <typename Getter>
py::array_t<double> LinkArray(const Network<Real>& network, Getter getter) {
  const auto& links = network.links();
  py::array_t<double> result(static_cast<py::ssize_t>(links.size()));
  auto buffer = result.mutable_unchecked<1>();
  for (py::ssize_t i = 0; i < static_cast<py::ssize_t>(links.size()); ++i) {
    buffer(i) = static_cast<double>(getter(links[static_cast<std::size_t>(i)]));
  }
  return result;
}

std::unique_ptr<Network<Real>> BuildNetworkFromDataset(const std::string& dataset,
                                                       const std::string& data_root) {
  NetworkBuilder builder;
  auto [nodes, zones, links] = builder.LoadNetworkData<Real>(dataset, data_root);
  auto trips = builder.LoadTripData<Real>(dataset, zones, data_root);
  auto [adjacency, reverse_adjacency] = builder.BuildAdjacencyLists<Real>(links, nodes);
  return std::make_unique<Network<Real>>(dataset, nodes, zones, std::move(links),
                                         std::move(trips), std::move(adjacency),
                                         std::move(reverse_adjacency));
}

std::unique_ptr<Network<Real>> BuildNetworkFromArrays(
    const std::string& name, int n_nodes, int n_zones,
    const IndexArray& init_node, const IndexArray& term_node,
    const DoubleArray& capacity, const DoubleArray& length,
    const DoubleArray& free_flow_time, const DoubleArray& b,
    const DoubleArray& power, const DoubleArray& speed, const DoubleArray& toll,
    const IndexArray& link_type, const DoubleArray& demand) {
  const auto init = init_node.unchecked<1>();
  const auto term = term_node.unchecked<1>();
  const py::ssize_t n_links = init.shape(0);

  auto check_length = [n_links](py::ssize_t size, const char* field) {
    if (size != n_links) {
      throw std::invalid_argument(std::string(field) +
                                  " must have one entry per link");
    }
  };
  check_length(term.shape(0), "term_node");
  check_length(capacity.shape(0), "capacity");
  check_length(length.shape(0), "length");
  check_length(free_flow_time.shape(0), "free_flow_time");
  check_length(b.shape(0), "b");
  check_length(power.shape(0), "power");
  check_length(speed.shape(0), "speed");
  check_length(toll.shape(0), "toll");
  check_length(link_type.shape(0), "link_type");

  if (n_zones <= 0 || n_nodes < n_zones) {
    throw std::invalid_argument(
        "need 0 < n_zones <= n_nodes (zone nodes are ids 0..n_zones-1)");
  }
  const auto demand_view = demand.unchecked<2>();
  if (demand_view.shape(0) != n_zones || demand_view.shape(1) != n_zones) {
    throw std::invalid_argument("demand must be an n_zones x n_zones matrix");
  }

  const auto cap = capacity.unchecked<1>();
  const auto len = length.unchecked<1>();
  const auto fft = free_flow_time.unchecked<1>();
  const auto b_view = b.unchecked<1>();
  const auto power_view = power.unchecked<1>();
  const auto speed_view = speed.unchecked<1>();
  const auto toll_view = toll.unchecked<1>();
  const auto type_view = link_type.unchecked<1>();

  std::vector<Link<Real>> links;
  links.reserve(static_cast<std::size_t>(n_links));
  for (py::ssize_t i = 0; i < n_links; ++i) {
    if (init(i) < 0 || init(i) >= n_nodes || term(i) < 0 || term(i) >= n_nodes) {
      throw std::invalid_argument("link " + std::to_string(i) +
                                  " has a node id outside [0, n_nodes)");
    }
    links.emplace_back(static_cast<int>(init(i)), static_cast<int>(term(i)),
                       static_cast<Real>(cap(i)), static_cast<Real>(len(i)),
                       static_cast<Real>(fft(i)), static_cast<Real>(b_view(i)),
                       static_cast<Real>(power_view(i)),
                       static_cast<Real>(speed_view(i)),
                       static_cast<Real>(toll_view(i)),
                       static_cast<int>(type_view(i)));
  }

  std::vector<std::vector<Real>> trips(static_cast<std::size_t>(n_zones),
                                       std::vector<Real>(static_cast<std::size_t>(n_zones), Real(0)));
  for (py::ssize_t origin = 0; origin < n_zones; ++origin) {
    for (py::ssize_t dest = 0; dest < n_zones; ++dest) {
      trips[static_cast<std::size_t>(origin)][static_cast<std::size_t>(dest)] =
          static_cast<Real>(demand_view(origin, dest));
    }
  }

  NetworkBuilder builder;
  auto [adjacency, reverse_adjacency] = builder.BuildAdjacencyLists<Real>(links, n_nodes);
  return std::make_unique<Network<Real>>(name, n_nodes, n_zones, std::move(links),
                                         std::move(trips), std::move(adjacency),
                                         std::move(reverse_adjacency));
}

}  // namespace

PYBIND11_MODULE(_core, m) {
  m.doc() = "Native TAP / CNDP solvers (long double instantiation of the C++ core)";

  py::class_<Link<Real>>(m, "Link",
                         "Directed link with BPR delay t(x) = t0*(1 + b*(x/c)^p)")
      .def_readonly("init", &Link<Real>::init, "Source node id (0-based)")
      .def_readonly("term", &Link<Real>::term, "Target node id (0-based)")
      .def_readonly("type", &Link<Real>::type)
      .def_readonly("length", &Link<Real>::length)
      .def_readonly("free_flow_time", &Link<Real>::free_flow_time)
      .def_readonly("b", &Link<Real>::b)
      .def_readonly("power", &Link<Real>::power)
      .def_readonly("speed", &Link<Real>::speed)
      .def_readonly("toll", &Link<Real>::toll)
      .def_readwrite("capacity", &Link<Real>::capacity)
      .def_readwrite("flow", &Link<Real>::flow)
      .def(
          "delay",
          [](const Link<Real>& link, double flow) {
            return static_cast<double>(link.Delay(static_cast<Real>(flow)));
          },
          py::arg("flow") = -1.0,
          "BPR travel time at the given flow (default: current link flow)")
      .def("__repr__", [](const Link<Real>& link) {
        return "<Link " + std::to_string(link.init) + "->" + std::to_string(link.term) +
               " capacity=" + std::to_string(static_cast<double>(link.capacity)) +
               " flow=" + std::to_string(static_cast<double>(link.flow)) + ">";
      });

  py::class_<Network<Real>>(m, "Network",
                            "Transportation network (links, OD pairs, adjacency)")
      .def_property_readonly("name", &Network<Real>::name)
      .def_property_readonly("number_of_links", &Network<Real>::number_of_links)
      .def_property_readonly("number_of_nodes", &Network<Real>::number_of_nodes)
      .def_property_readonly("number_of_zones", &Network<Real>::number_of_zones)
      .def_property_readonly("number_of_od_pairs", &Network<Real>::number_of_od_pairs)
      .def(
          "total_travel_time",
          [](const Network<Real>& network) {
            return static_cast<double>(network.TotalTravelTime());
          },
          "Total system travel time: sum_a f_a * c_a(f_a)")
      .def(
          "beckmann_objective",
          [](const Network<Real>& network) {
            return static_cast<double>(network.ObjectiveFunction());
          },
          "Beckmann user-equilibrium objective: sum_a integral_0^{f_a} t_a(x) dx")
      .def(
          "relative_gap",
          [](const Network<Real>& network) {
            py::gil_scoped_release release;
            return static_cast<double>(network.RelativeGap());
          },
          "RGAP convergence measure (0 at user equilibrium); runs one shortest-path sweep")
      .def("reset", &Network<Real>::Reset,
           "Zero all link flows and clear OD route sets")
      .def(
          "link",
          [](Network<Real>& network, int index) -> Link<Real>& {
            if (index < 0 || index >= network.number_of_links()) {
              throw py::index_error("link index out of range");
            }
            return network.mutable_links()[static_cast<std::size_t>(index)];
          },
          py::arg("index"), py::return_value_policy::reference_internal,
          "Mutable view of one link")
      .def(
          "capacities",
          [](const Network<Real>& network) {
            return LinkArray(network, [](const Link<Real>& link) { return link.capacity; });
          },
          "Per-link capacities as a float array")
      .def(
          "flows",
          [](const Network<Real>& network) {
            return LinkArray(network, [](const Link<Real>& link) { return link.flow; });
          },
          "Per-link flows as a float array")
      .def(
          "free_flow_times",
          [](const Network<Real>& network) {
            return LinkArray(network,
                             [](const Link<Real>& link) { return link.free_flow_time; });
          })
      .def(
          "link_costs",
          [](const Network<Real>& network) {
            return LinkArray(network, [](const Link<Real>& link) { return link.Delay(); });
          },
          "Per-link BPR travel times at current flows")
      .def(
          "link_nodes",
          [](const Network<Real>& network) {
            const auto& links = network.links();
            py::array_t<std::int64_t> result(
                {static_cast<py::ssize_t>(links.size()), static_cast<py::ssize_t>(2)});
            auto buffer = result.mutable_unchecked<2>();
            for (py::ssize_t i = 0; i < static_cast<py::ssize_t>(links.size()); ++i) {
              buffer(i, 0) = links[static_cast<std::size_t>(i)].init;
              buffer(i, 1) = links[static_cast<std::size_t>(i)].term;
            }
            return result;
          },
          "(n_links, 2) array of [init, term] node ids")
      .def(
          "set_capacities",
          [](Network<Real>& network, const DoubleArray& capacities) {
            const auto values = capacities.unchecked<1>();
            if (values.shape(0) != network.number_of_links()) {
              throw std::invalid_argument("capacities must have one entry per link");
            }
            auto& links = network.mutable_links();
            for (py::ssize_t i = 0; i < values.shape(0); ++i) {
              links[static_cast<std::size_t>(i)].capacity = static_cast<Real>(values(i));
            }
          },
          py::arg("capacities"), "Bulk-set per-link capacities")
      .def("__repr__", [](const Network<Real>& network) {
        return "<Network '" + network.name() + "' links=" +
               std::to_string(network.number_of_links()) +
               " nodes=" + std::to_string(network.number_of_nodes()) +
               " od_pairs=" + std::to_string(network.number_of_od_pairs()) + ">";
      });

  py::class_<TrafficAssignmentApproach<Real>,
             std::shared_ptr<TrafficAssignmentApproach<Real>>>(
      m, "TrafficAssignmentApproach", "Abstract TAP solver (user equilibrium)")
      .def_property_readonly("approach_name",
                             &TrafficAssignmentApproach<Real>::GetApproachName)
      .def_property_readonly("network", &TrafficAssignmentApproach<Real>::network,
                             py::return_value_policy::reference_internal)
      .def(
          "compute_traffic_flows",
          [](TrafficAssignmentApproach<Real>& approach, bool statistics_recording) {
            py::gil_scoped_release release;
            approach.ComputeTrafficFlows(statistics_recording);
          },
          py::arg("statistics_recording") = false,
          "Solve user equilibrium (releases the GIL)")
      .def("reset", &TrafficAssignmentApproach<Real>::Reset,
           "Reset the network to its initial state (zero flows, no routes)")
      .def("set_output_root", &TrafficAssignmentApproach<Real>::SetOutputRoot,
           py::arg("path"), "Output directory for optional statistics recording")
      .def("set_route_search_threads",
           &TrafficAssignmentApproach<Real>::SetRouteSearchThreadCount,
           py::arg("thread_count"));

  py::class_<TapasApproach<Real>, TrafficAssignmentApproach<Real>,
             std::shared_ptr<TapasApproach<Real>>>(
      m, "TapasApproach", "Bush-based TAPAS solver (fastest, machine-precision RGAP)")
      .def(py::init<Network<Real>&, Real, int, Real, Real>(), py::arg("network"),
           py::arg("alpha") = static_cast<Real>(1e-14L), py::arg("max_iterations") = 200,
           py::arg("mu") = static_cast<Real>(0.5L), py::arg("v") = static_cast<Real>(0.25L),
           py::keep_alive<1, 2>());

  py::class_<RouteBasedApproach<Real>, TrafficAssignmentApproach<Real>,
             std::shared_ptr<RouteBasedApproach<Real>>>(
      m, "RouteBasedApproach", "Path-based solver with explicit route enumeration")
      .def(py::init<Network<Real>&, Real, std::string, std::size_t, int, int, int, Real>(),
           py::arg("network"), py::arg("alpha") = static_cast<Real>(1e-14L),
           py::arg("shift_method") = "NewtonStep",
           py::arg("route_search_threads") = static_cast<std::size_t>(1),
           py::arg("max_iterations") = 200, py::arg("full_iteration_count") = 3,
           py::arg("origin_iteration_count") = 1,
           py::arg("ema_alpha") = static_cast<Real>(0.7L), py::keep_alive<1, 2>());

  py::class_<DirectedLinkCapacityConstraint>(
      m, "LinkConstraint", "Per-link capacity bounds [lower, upper] for CNDP")
      .def(py::init<std::size_t, std::size_t, double, double, double>(),
           py::arg("init_node"), py::arg("term_node"), py::arg("lower_bound"),
           py::arg("upper_bound"), py::arg("investment_cost") = 1.0)
      .def_readwrite("init_node", &DirectedLinkCapacityConstraint::init_node)
      .def_readwrite("term_node", &DirectedLinkCapacityConstraint::term_node)
      .def_readwrite("lower_bound", &DirectedLinkCapacityConstraint::lower_bound)
      .def_readwrite("upper_bound", &DirectedLinkCapacityConstraint::upper_bound)
      .def_readwrite("investment_cost_param",
                     &DirectedLinkCapacityConstraint::investment_cost_param)
      .def("__repr__", [](const DirectedLinkCapacityConstraint& constraint) {
        return "<LinkConstraint " + constraint.ToString() + " [" +
               std::to_string(constraint.lower_bound) + ", " +
               std::to_string(constraint.upper_bound) + "]>";
      });

  py::class_<OptimizationStepConfig>(
      m, "StepConfig", "One CNDP pipeline step (equivalent of a [[pipeline]] TOML entry)")
      .def(py::init<>())
      .def_readwrite("type", &OptimizationStepConfig::type)
      .def_readwrite("name", &OptimizationStepConfig::name)
      .def_readwrite("max_iterations", &OptimizationStepConfig::max_iterations)
      .def_readwrite("tolerance", &OptimizationStepConfig::tolerance)
      .def_readwrite("algorithm", &OptimizationStepConfig::algorithm)
      .def_readwrite("local_algorithm", &OptimizationStepConfig::local_algorithm)
      .def_readwrite("local_max_iterations", &OptimizationStepConfig::local_max_iterations)
      .def_readwrite("local_tolerance", &OptimizationStepConfig::local_tolerance)
      .def_readwrite("population_size", &OptimizationStepConfig::population_size)
      .def_readwrite("step_size", &OptimizationStepConfig::step_size)
      .def_readwrite("fd_epsilon", &OptimizationStepConfig::fd_epsilon)
      .def_readwrite("gradient_method", &OptimizationStepConfig::gradient_method)
      .def_readwrite("stochastic_optimizer", &OptimizationStepConfig::stochastic_optimizer)
      .def("__repr__", [](const OptimizationStepConfig& config) {
        std::string text = "<StepConfig type='" + config.type + "'";
        if (!config.algorithm.empty()) text += " algorithm='" + config.algorithm + "'";
        text += " max_iterations=" + std::to_string(config.max_iterations) + ">";
        return text;
      });

  py::class_<CndMetricsConfig>(
      m, "MetricsConfig", "File-output options for CNDP runs ([metrics] TOML section)")
      .def(py::init<>())
      .def_readwrite("enable_trace", &CndMetricsConfig::enable_trace)
      .def_readwrite("enable_relative_gap", &CndMetricsConfig::enable_relative_gap)
      .def_readwrite("relative_gap_sample_period",
                     &CndMetricsConfig::relative_gap_sample_period)
      .def_readwrite("flush_every_n_points", &CndMetricsConfig::flush_every_n_points)
      .def_readwrite("write_metadata_json", &CndMetricsConfig::write_metadata_json)
      .def_readwrite("write_summary_csv", &CndMetricsConfig::write_summary_csv)
      .def_readwrite("append_dataset_subdir", &CndMetricsConfig::append_dataset_subdir)
      .def_readwrite("output_root", &CndMetricsConfig::output_root)
      .def_readwrite("run_id", &CndMetricsConfig::run_id)
      .def_readwrite("scenario_name", &CndMetricsConfig::scenario_name);

  py::class_<BilevelCND<Real>>(
      m, "BilevelCND",
      "Bilevel CNDP solver: pipeline of optimization steps over link capacities")
      .def(py::init<Network<Real>&, std::shared_ptr<TrafficAssignmentApproach<Real>>,
                    const std::vector<DirectedLinkCapacityConstraint>&,
                    const std::vector<OptimizationStepConfig>&, Real, Real, Real, double,
                    const CndMetricsConfig&, std::size_t, const std::string&>(),
           py::arg("network"), py::arg("approach"), py::arg("constraints"),
           py::arg("pipeline"),
           py::arg("link_capacity_selection_threshold") = static_cast<Real>(1e-3L),
           py::arg("budget_threshold") = static_cast<Real>(1e-1L),
           py::arg("budget_function_multiplier") = static_cast<Real>(5.0L),
           py::arg("budget_upper_bound") = 100000.0,
           py::arg("metrics") = CndMetricsConfig(),
           py::arg("route_search_threads") = static_cast<std::size_t>(1),
           py::arg("progress_format") = "none", py::keep_alive<1, 2>())
      .def("set_verbose", &BilevelCND<Real>::SetVerbose, py::arg("verbose"))
      .def("set_statistics_enabled", &BilevelCND<Real>::SetStatisticsEnabled,
           py::arg("enabled"))
      .def("set_final_diagnostics_enabled", &BilevelCND<Real>::SetFinalDiagnosticsEnabled,
           py::arg("enabled"))
      .def("set_route_search_threads", &BilevelCND<Real>::SetRouteSearchThreadCount,
           py::arg("thread_count"))
      .def(
          "compute_network_design",
          [](BilevelCND<Real>& solver) {
            py::gil_scoped_release release;
            solver.ComputeNetworkDesign();
          },
          "Run the optimization pipeline (releases the GIL)");

  m.def(
      "_build_network_from_dataset",
      [](const std::string& dataset, const std::string& data_root) {
        py::gil_scoped_release release;
        return BuildNetworkFromDataset(dataset, data_root);
      },
      py::arg("dataset"), py::arg("data_root"));

  m.def("_build_network_from_arrays", &BuildNetworkFromArrays, py::arg("name"),
        py::arg("n_nodes"), py::arg("n_zones"), py::arg("init_node"),
        py::arg("term_node"), py::arg("capacity"), py::arg("length"),
        py::arg("free_flow_time"), py::arg("b"), py::arg("power"), py::arg("speed"),
        py::arg("toll"), py::arg("link_type"), py::arg("demand"));

  m.def(
      "_load_constraints",
      [](const std::string& path, bool verbose) {
        DirectedConstraintLoader loader;
        loader.SetVerbose(verbose);
        return loader.LoadFromFile(path);
      },
      py::arg("path"), py::arg("verbose") = false);
}

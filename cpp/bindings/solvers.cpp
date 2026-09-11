#include "bindings.hpp"

void BindSolvers(py::module_& m) {
  py::class_<TrafficAssignmentApproach, std::shared_ptr<TrafficAssignmentApproach>>(
      m, "TrafficAssignmentApproach")
      .def_property_readonly("approach_name", &TrafficAssignmentApproach::GetApproachName)
      .def_property_readonly("network", &TrafficAssignmentApproach::network)
      .def("compute_traffic_flows", &TrafficAssignmentApproach::ComputeTrafficFlows,
           py::arg("statistics_recording") = false, py::call_guard<py::gil_scoped_release>())
      .def("reset", &TrafficAssignmentApproach::Reset)
      .def("set_output_root", &TrafficAssignmentApproach::SetOutputRoot, py::arg("path"))
      .def("set_route_search_threads", &TrafficAssignmentApproach::SetRouteSearchThreadCount,
           py::arg("thread_count"));

  py::class_<TapasApproach, TrafficAssignmentApproach, std::shared_ptr<TapasApproach>>(
      m, "TapasApproach")
      .def(py::init<std::shared_ptr<Network>, Real, int, Real, Real>(), py::arg("network"),
           py::arg("alpha") = 1e-14L, py::arg("max_iterations") = 200,
           py::arg("mu") = 0.5L, py::arg("v") = 0.25L);

  py::class_<RouteBasedApproach, TrafficAssignmentApproach, std::shared_ptr<RouteBasedApproach>>(
      m, "RouteBasedApproach")
      .def(py::init<std::shared_ptr<Network>, Real, const std::string&, std::size_t,
                    int, int, int, Real>(), py::arg("network"), py::arg("alpha") = 1e-14L,
           py::arg("shift_method") = "NewtonStep", py::arg("route_search_threads") = 1,
           py::arg("max_iterations") = 200, py::arg("full_iteration_count") = 3,
           py::arg("origin_iteration_count") = 1, py::arg("ema_alpha") = 0.7L);

  py::class_<BilevelCND>(m, "BilevelCND")
      .def(py::init<std::shared_ptr<Network>, std::shared_ptr<TrafficAssignmentApproach>,
                    const std::vector<LinkConstraint>&, const std::vector<StepConfig>&,
                    const CndpOptions&>(),
           py::arg("network"), py::arg("approach"), py::arg("constraints"),
           py::arg("pipeline"), py::arg("options"))
      // Retain the original constructor for existing object-level Python workflows.
      .def(py::init([](std::shared_ptr<Network> network,
                      std::shared_ptr<TrafficAssignmentApproach> approach,
                      const std::vector<LinkConstraint>& constraints,
                      const std::vector<StepConfig>& pipeline, Real link_threshold,
                      Real budget_threshold, Real theta, double budget,
                      const MetricsConfig& metrics, std::size_t threads,
                      const std::string& progress) {
        CndpOptions options;
        options.link_capacity_selection_threshold = link_threshold;
        options.budget_threshold = budget_threshold;
        options.budget_function_multiplier = theta;
        options.budget_upper_bound = budget;
        options.metrics = metrics;
        options.route_search_threads = threads;
        options.progress_format = progress;
        // Historically this constructor accepts already prepared constraints.
        options.exclude_insensitive_links = false;
        return std::make_unique<BilevelCND>(network, approach, constraints, pipeline, options);
      }), py::arg("network"), py::arg("approach"), py::arg("constraints"), py::arg("pipeline"),
          py::arg("link_capacity_selection_threshold") = 1e-3L,
          py::arg("budget_threshold") = 1e-1L, py::arg("budget_function_multiplier") = 5.0L,
          py::arg("budget_upper_bound") = 100000.0, py::arg("metrics") = MetricsConfig(),
          py::arg("route_search_threads") = 1, py::arg("progress_format") = "none")
      .def("set_verbose", &BilevelCND::SetVerbose, py::arg("verbose"))
      .def("set_statistics_enabled", &BilevelCND::SetStatisticsEnabled, py::arg("enabled"))
      .def("set_final_diagnostics_enabled", &BilevelCND::SetFinalDiagnosticsEnabled,
           py::arg("enabled"))
      .def("set_route_search_threads", &BilevelCND::SetRouteSearchThreadCount,
           py::arg("thread_count"))
      .def_property_readonly("active_constraint_count", &BilevelCND::active_constraint_count)
      .def("compute_network_design", &BilevelCND::ComputeNetworkDesign,
           py::call_guard<py::gil_scoped_release>());

  m.def("_make_approach", &MakeApproach, py::arg("network"), py::arg("options"));
  m.def("_solve_tap", py::overload_cast<TrafficAssignmentApproach&, bool, bool>(&SolveTap),
        py::arg("approach"), py::arg("reset") = false,
        py::arg("statistics_recording") = false, py::call_guard<py::gil_scoped_release>());
}

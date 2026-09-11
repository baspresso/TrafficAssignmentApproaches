#include "bindings.hpp"

void BindResults(py::module_& m) {
  py::class_<TapResult>(m, "NativeTapResult")
      .def_readonly("approach", &TapResult::approach)
      .def_readonly("network_name", &TapResult::network_name)
      .def_readonly("relative_gap", &TapResult::relative_gap)
      .def_readonly("total_travel_time", &TapResult::total_travel_time)
      .def_readonly("beckmann_objective", &TapResult::beckmann_objective)
      .def_readonly("solve_seconds", &TapResult::solve_seconds)
      .def_property_readonly("flows", [](const TapResult& r) { return ToArray(r.flows); })
      .def_property_readonly("link_costs", [](const TapResult& r) { return ToArray(r.link_costs); });

  py::class_<CndpResult>(m, "NativeCndpResult")
      .def_readonly("objective", &CndpResult::objective)
      .def_readonly("total_travel_time", &CndpResult::total_travel_time)
      .def_readonly("budget", &CndpResult::budget)
      .def_readonly("budget_upper_bound", &CndpResult::budget_upper_bound)
      .def_readonly("elapsed_seconds", &CndpResult::elapsed_seconds)
      .def_property_readonly("capacities", [](const CndpResult& r) { return ToArray(r.capacities); })
      .def_property_readonly("flows", [](const CndpResult& r) { return ToArray(r.flows); })
      .def_property_readonly("lower_bounds", [](const CndpResult& r) { return ToArray(r.lower_bounds); })
      .def_property_readonly("upper_bounds", [](const CndpResult& r) { return ToArray(r.upper_bounds); });
}

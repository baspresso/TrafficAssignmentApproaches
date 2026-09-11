#include "bindings.hpp"

void BindOptions(py::module_& m) {
  py::class_<LinkConstraint>(
      m, "LinkConstraint", "Per-link capacity bounds [lower, upper] for CNDP")
      .def(py::init<std::size_t, std::size_t, double, double, double>(),
           py::arg("init_node"), py::arg("term_node"), py::arg("lower_bound"),
           py::arg("upper_bound"), py::arg("investment_cost") = 1.0)
      .def_readwrite("init_node", &LinkConstraint::init_node)
      .def_readwrite("term_node", &LinkConstraint::term_node)
      .def_readwrite("lower_bound", &LinkConstraint::lower_bound)
      .def_readwrite("upper_bound", &LinkConstraint::upper_bound)
      .def_readwrite("investment_cost_param",
                     &LinkConstraint::investment_cost_param)
      .def("__repr__", [](const LinkConstraint& constraint) {
        return "<LinkConstraint " + constraint.ToString() + " [" +
               std::to_string(constraint.lower_bound) + ", " +
               std::to_string(constraint.upper_bound) + "]>";
      });

  py::class_<StepConfig>(
      m, "StepConfig", "One CNDP pipeline step (equivalent of a [[pipeline]] TOML entry)")
      .def(py::init<>())
      .def_readwrite("type", &StepConfig::type)
      .def_readwrite("name", &StepConfig::name)
      .def_readwrite("max_iterations", &StepConfig::max_iterations)
      .def_readwrite("tolerance", &StepConfig::tolerance)
      .def_readwrite("algorithm", &StepConfig::algorithm)
      .def_readwrite("local_algorithm", &StepConfig::local_algorithm)
      .def_readwrite("local_max_iterations", &StepConfig::local_max_iterations)
      .def_readwrite("local_tolerance", &StepConfig::local_tolerance)
      .def_readwrite("population_size", &StepConfig::population_size)
      .def_readwrite("step_size", &StepConfig::step_size)
      .def_readwrite("fd_epsilon", &StepConfig::fd_epsilon)
      .def_readwrite("gradient_method", &StepConfig::gradient_method)
      .def_readwrite("stochastic_optimizer", &StepConfig::stochastic_optimizer)
      .def("__repr__", [](const StepConfig& config) {
        std::string text = "<StepConfig type='" + config.type + "'";
        if (!config.algorithm.empty()) text += " algorithm='" + config.algorithm + "'";
        text += " max_iterations=" + std::to_string(config.max_iterations) + ">";
        return text;
      });

  py::class_<MetricsConfig>(
      m, "MetricsConfig", "File-output options for CNDP runs ([metrics] TOML section)")
      .def(py::init<>())
      .def_readwrite("enable_trace", &MetricsConfig::enable_trace)
      .def_readwrite("enable_relative_gap", &MetricsConfig::enable_relative_gap)
      .def_readwrite("relative_gap_sample_period",
                     &MetricsConfig::relative_gap_sample_period)
      .def_readwrite("flush_every_n_points", &MetricsConfig::flush_every_n_points)
      .def_readwrite("write_metadata_json", &MetricsConfig::write_metadata_json)
      .def_readwrite("write_summary_csv", &MetricsConfig::write_summary_csv)
      .def_readwrite("append_dataset_subdir", &MetricsConfig::append_dataset_subdir)
      .def_readwrite("output_root", &MetricsConfig::output_root)
      .def_readwrite("run_id", &MetricsConfig::run_id)
      .def_readwrite("scenario_name", &MetricsConfig::scenario_name);

  py::class_<TapOptions>(m, "TapOptions")
      .def(py::init<>())
      .def(py::init<const TapOptions&>())
      .def_readwrite("approach", &TapOptions::approach)
      .def_readwrite("relative_gap_tolerance", &TapOptions::relative_gap_tolerance)
      .def_readwrite("max_iterations", &TapOptions::max_iterations)
      .def_readwrite("mu", &TapOptions::mu)
      .def_readwrite("v", &TapOptions::v)
      .def_readwrite("shift_method", &TapOptions::shift_method)
      .def_readwrite("route_search_threads", &TapOptions::route_search_threads)
      .def_readwrite("full_iteration_count", &TapOptions::full_iteration_count)
      .def_readwrite("origin_iteration_count", &TapOptions::origin_iteration_count)
      .def_readwrite("ema_alpha", &TapOptions::ema_alpha);

  py::class_<CndpOptions>(m, "CndpOptions")
      .def(py::init<>())
      .def(py::init<const CndpOptions&>())
      .def_readwrite("link_capacity_selection_threshold", &CndpOptions::link_capacity_selection_threshold)
      .def_readwrite("budget_threshold", &CndpOptions::budget_threshold)
      .def_readwrite("budget_function_multiplier", &CndpOptions::budget_function_multiplier)
      .def_readwrite("budget_upper_bound", &CndpOptions::budget_upper_bound)
      .def_readwrite("route_search_threads", &CndpOptions::route_search_threads)
      .def_readwrite("exclude_insensitive_links", &CndpOptions::exclude_insensitive_links)
      .def_readwrite("statistics_enabled", &CndpOptions::statistics_enabled)
      .def_readwrite("final_diagnostics", &CndpOptions::final_diagnostics)
      .def_readwrite("verbose", &CndpOptions::verbose)
      .def_readwrite("progress_format", &CndpOptions::progress_format)
      .def_readwrite("metrics", &CndpOptions::metrics);

}

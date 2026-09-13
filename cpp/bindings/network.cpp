#include "bindings.hpp"
#include <limits>

namespace {
std::shared_ptr<Network> BuildNetworkFromArrays(
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
  for (const auto* values : {&capacity, &length, &free_flow_time, &b, &power, &speed, &toll}) {
    if (values->ndim() != 1) throw std::invalid_argument("per-link arrays must be one-dimensional");
  }
  if (link_type.ndim() != 1) throw std::invalid_argument("link_type must be one-dimensional");
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

  std::vector<LinkData> links;
  links.reserve(static_cast<std::size_t>(n_links));
  for (py::ssize_t i = 0; i < n_links; ++i) {
    if (init(i) < 0 || init(i) >= n_nodes || term(i) < 0 || term(i) >= n_nodes) {
      throw std::invalid_argument("link " + std::to_string(i) +
                                  " has a node id outside [0, n_nodes)");
    }
    if (type_view(i) < 0 || type_view(i) > std::numeric_limits<int>::max()) {
      throw std::invalid_argument("link_type must be an integer in [0, 2147483647]");
    }
    links.push_back(LinkData{static_cast<int>(init(i)), static_cast<int>(term(i)),
                       static_cast<Real>(cap(i)), static_cast<Real>(len(i)),
                       static_cast<Real>(fft(i)), static_cast<Real>(b_view(i)),
                       static_cast<Real>(power_view(i)),
                       static_cast<Real>(speed_view(i)),
                       static_cast<Real>(toll_view(i)),
                       static_cast<int>(type_view(i))});
  }

  std::vector<std::vector<Real>> trips(static_cast<std::size_t>(n_zones),
                                       std::vector<Real>(static_cast<std::size_t>(n_zones), Real(0)));
  for (py::ssize_t origin = 0; origin < n_zones; ++origin) {
    for (py::ssize_t dest = 0; dest < n_zones; ++dest) {
      trips[static_cast<std::size_t>(origin)][static_cast<std::size_t>(dest)] =
          static_cast<Real>(demand_view(origin, dest));
    }
  }

  return std::make_shared<Network>(NetworkData{
      name, n_nodes, n_zones, std::move(links), std::move(trips)});
}
}

void BindNetwork(py::module_& m) {
  py::class_<Link>(m, "Link")
      .def_property_readonly("init", [](const Link& l) { return l.data().init; })
      .def_property_readonly("term", [](const Link& l) { return l.data().term; })
      .def_property_readonly("type", [](const Link& l) { return l.data().type; })
      .def_property_readonly("length", [](const Link& l) { return l.data().length; })
      .def_property_readonly("free_flow_time", [](const Link& l) { return l.data().free_flow_time; })
      .def_property_readonly("b", [](const Link& l) { return l.data().b; })
      .def_property_readonly("power", [](const Link& l) { return l.data().power; })
      .def_property_readonly("speed", [](const Link& l) { return l.data().speed; })
      .def_property_readonly("toll", [](const Link& l) { return l.data().toll; })
      .def_property("capacity", &Link::capacity, &Link::set_capacity)
      .def_property("flow", &Link::flow, &Link::set_flow)
      .def("delay", &Link::delay, py::arg("flow") = -1.0L)
      .def("__repr__", [](const Link& l) {
        const auto d = l.data();
        return "<Link " + std::to_string(d.init) + "->" + std::to_string(d.term) +
               " capacity=" + std::to_string(static_cast<double>(l.capacity())) +
               " flow=" + std::to_string(static_cast<double>(l.flow())) + ">";
      });

  py::class_<Network, std::shared_ptr<Network>>(m, "Network")
      .def_property_readonly("name", &Network::name)
      .def_property_readonly("number_of_links", &Network::number_of_links)
      .def_property_readonly("number_of_nodes", &Network::number_of_nodes)
      .def_property_readonly("number_of_zones", &Network::number_of_zones)
      .def_property_readonly("number_of_od_pairs", &Network::number_of_od_pairs)
      .def("total_travel_time", &Network::total_travel_time)
      .def("beckmann_objective", &Network::beckmann_objective)
      .def("relative_gap", &Network::relative_gap, py::call_guard<py::gil_scoped_release>())
      .def("reset", &Network::reset)
      .def("link", &Network::link, py::arg("index"))
      .def("capacities", [](const Network& n) { return ToArray(n.capacities()); })
      .def("flows", [](const Network& n) { return ToArray(n.flows()); })
      .def("free_flow_times", [](const Network& n) { return ToArray(n.free_flow_times()); })
      .def("link_costs", [](const Network& n) { return ToArray(n.link_costs()); })
      .def("link_nodes", [](const Network& n) {
        const auto nodes = n.link_nodes();
        py::array_t<std::int64_t> result(
            {static_cast<py::ssize_t>(nodes.size()), static_cast<py::ssize_t>(2)});
        auto buffer = result.mutable_unchecked<2>();
        for (std::size_t i = 0; i < nodes.size(); ++i) {
          buffer(i, 0) = nodes[i][0];
          buffer(i, 1) = nodes[i][1];
        }
        return result;
      })
      .def("set_capacities", [](Network& n, const DoubleArray& values) {
        const auto view = values.unchecked<1>();
        std::vector<Real> capacities(static_cast<std::size_t>(view.shape(0)));
        for (py::ssize_t i = 0; i < view.shape(0); ++i) capacities[i] = view(i);
        n.set_capacities(capacities);
      }, py::arg("capacities"))
      .def("__repr__", [](const Network& n) {
        return "<Network '" + n.name() + "' links=" + std::to_string(n.number_of_links()) +
               " nodes=" + std::to_string(n.number_of_nodes()) +
               " od_pairs=" + std::to_string(n.number_of_od_pairs()) + ">";
      });

  m.def("_build_network_from_dataset", &LoadNetwork, py::arg("dataset"), py::arg("data_root"),
        py::call_guard<py::gil_scoped_release>());
  m.def("_load_constraints", &LoadConstraints, py::arg("path"), py::arg("verbose") = false,
        py::arg("node_index_base") = 1,
        py::call_guard<py::gil_scoped_release>());
  m.def("_build_network_from_arrays", &BuildNetworkFromArrays, py::arg("name"),
        py::arg("n_nodes"), py::arg("n_zones"), py::arg("init_node"),
        py::arg("term_node"), py::arg("capacity"), py::arg("length"),
        py::arg("free_flow_time"), py::arg("b"), py::arg("power"), py::arg("speed"),
        py::arg("toll"), py::arg("link_type"), py::arg("demand"));
}

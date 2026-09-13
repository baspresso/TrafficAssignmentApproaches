#include <traffic_assignment/network.hpp>

#include <stdexcept>
#include <limits>
#include <utility>

#include "detail/access.hpp"
#include "detail/input_validation.hpp"
#include "detail/tap/core/NetworkBuilder.h"
#include "detail/cnd/DirectedConstraintLoader.h"

namespace traffic_assignment {
namespace {
template <typename Getter>
std::vector<Real> LinkValues(const detail::NetworkImpl& impl, Getter getter) {
  std::vector<Real> values;
  values.reserve(impl.network->number_of_links());
  for (const auto& link : impl.network->links()) values.push_back(getter(link));
  return values;
}
}

Network::Network(const NetworkData& data) : impl_(std::make_shared<detail::NetworkImpl>()) {
  if (data.number_of_zones <= 0 || data.number_of_nodes < data.number_of_zones) {
    throw std::invalid_argument("need 0 < n_zones <= n_nodes (zone nodes are ids 0..n_zones-1)");
  }
  if (data.demand.size() != static_cast<std::size_t>(data.number_of_zones)) {
    throw std::invalid_argument("demand must have one row per zone");
  }
  for (std::size_t origin = 0; origin < data.demand.size(); ++origin) {
    const auto& row = data.demand[origin];
    if (row.size() != data.demand.size()) {
      throw std::invalid_argument("demand must be an n_zones x n_zones matrix");
    }
    for (std::size_t dest = 0; dest < row.size(); ++dest) {
      detail::ValidateNonnegative(row[dest], "demand[" + std::to_string(origin) +
                                  ", " + std::to_string(dest) + "]");
    }
  }
  if (data.links.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::invalid_argument("too many links for native indices");
  }
  std::vector<TrafficAssignment::Link<Real>> links;
  links.reserve(data.links.size());
  for (std::size_t i = 0; i < data.links.size(); ++i) {
    const auto& link = data.links[i];
    const auto field = "link[" + std::to_string(i) + "]";
    if (link.init < 0 || link.init >= data.number_of_nodes ||
        link.term < 0 || link.term >= data.number_of_nodes) {
      throw std::invalid_argument(field + " has a node id outside [0, n_nodes)");
    }
    detail::ValidateNonnegative(link.capacity, field + ".capacity", true);
    detail::ValidateNonnegative(link.free_flow_time, field + ".free_flow_time");
    detail::ValidateNonnegative(link.b, field + ".b");
    detail::ValidateNonnegative(link.power, field + ".power");
    detail::ValidateNonnegative(link.length, field + ".length");
    detail::ValidateNonnegative(link.speed, field + ".speed");
    detail::ValidateNonnegative(link.toll, field + ".toll");
    if (link.type < 0) throw std::invalid_argument(field + ".link_type must be nonnegative");
    links.emplace_back(link.init, link.term, link.capacity, link.length,
                       link.free_flow_time, link.b, link.power, link.speed,
                       link.toll, link.type);
  }
  TrafficAssignment::NetworkBuilder builder;
  auto [adjacency, reverse] = builder.BuildAdjacencyLists<Real>(links, data.number_of_nodes);
  impl_->network = std::make_unique<TrafficAssignment::Network<Real>>(
      data.name, data.number_of_nodes, data.number_of_zones, std::move(links),
      data.demand, std::move(adjacency), std::move(reverse));
}

Network::Network(std::shared_ptr<detail::NetworkImpl> impl) : impl_(std::move(impl)) {}
Network::~Network() = default;
std::string Network::name() const { return impl_->network->name(); }
int Network::number_of_nodes() const { return impl_->network->number_of_nodes(); }
int Network::number_of_zones() const { return impl_->network->number_of_zones(); }
int Network::number_of_links() const { return impl_->network->number_of_links(); }
int Network::number_of_od_pairs() const { return impl_->network->number_of_od_pairs(); }
Real Network::total_travel_time() const { return impl_->network->TotalTravelTime(); }
Real Network::beckmann_objective() const { return impl_->network->ObjectiveFunction(); }
Real Network::relative_gap() const { return impl_->network->RelativeGap(); }
void Network::reset() { impl_->network->Reset(); }

Link Network::link(int index) const {
  if (index < 0 || index >= number_of_links()) throw std::out_of_range("link index out of range");
  return Link(impl_, static_cast<std::size_t>(index));
}

std::vector<Real> Network::flows() const {
  return LinkValues(*impl_, [](const auto& link) { return link.flow; });
}
std::vector<Real> Network::capacities() const {
  return LinkValues(*impl_, [](const auto& link) { return link.capacity; });
}
std::vector<Real> Network::free_flow_times() const {
  return LinkValues(*impl_, [](const auto& link) { return link.free_flow_time; });
}
std::vector<Real> Network::link_costs() const {
  return LinkValues(*impl_, [](const auto& link) { return link.Delay(); });
}
std::vector<std::array<int, 2>> Network::link_nodes() const {
  std::vector<std::array<int, 2>> nodes;
  nodes.reserve(number_of_links());
  for (const auto& link : impl_->network->links()) nodes.push_back({link.init, link.term});
  return nodes;
}
void Network::set_capacities(const std::vector<Real>& values) {
  if (values.size() != static_cast<std::size_t>(number_of_links())) {
    throw std::invalid_argument("capacities must have one entry per link");
  }
  for (const auto value : values) detail::ValidateNonnegative(value, "capacity", true);
  for (std::size_t i = 0; i < values.size(); ++i) {
    impl_->network->mutable_links()[i].capacity = values[i];
  }
}

Link::Link(std::shared_ptr<detail::NetworkImpl> impl, std::size_t index)
    : impl_(std::move(impl)), index_(index) {}
LinkData Link::data() const {
  const auto& link = impl_->network->links()[index_];
  return {link.init, link.term, link.capacity, link.length, link.free_flow_time,
          link.b, link.power, link.speed, link.toll, link.type};
}
Real Link::capacity() const { return impl_->network->links()[index_].capacity; }
Real Link::flow() const { return impl_->network->links()[index_].flow; }
void Link::set_capacity(Real value) {
  detail::ValidateNonnegative(value, "capacity", true);
  impl_->network->mutable_links()[index_].capacity = value;
}
void Link::set_flow(Real value) { impl_->network->mutable_links()[index_].flow = value; }
Real Link::delay(Real flow) const { return impl_->network->links()[index_].Delay(flow); }

std::shared_ptr<Network> LoadNetwork(const std::string& dataset, const std::string& data_root) {
  TrafficAssignment::NetworkBuilder builder;
  auto [nodes, zones, links] = builder.LoadNetworkData(dataset, data_root);
  auto trips = builder.LoadTripData<Real>(dataset, zones, data_root);
  return std::make_shared<Network>(NetworkData{
      dataset, nodes, zones, std::move(links), std::move(trips)});
}

std::vector<LinkConstraint> LoadConstraints(const std::string& path, bool verbose, int node_index_base) {
  TrafficAssignment::DirectedConstraintLoader loader;
  loader.SetVerbose(verbose);
  return loader.LoadFromFile(path, node_index_base);
}

}  // namespace traffic_assignment

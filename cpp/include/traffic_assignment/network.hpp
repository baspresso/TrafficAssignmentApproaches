#pragma once

#include <traffic_assignment/options.hpp>

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace traffic_assignment {
namespace detail {
struct NetworkImpl;
struct Access;
}

/// One directed link in an input network. Node identifiers are zero-based.
struct LinkData {
  int init = 0;
  int term = 0;
  Real capacity = 0;
  Real length = 0;
  Real free_flow_time = 0;
  Real b = 0;
  Real power = 1;
  Real speed = 0;
  Real toll = 0;
  int type = 1;
};

struct NetworkData {
  std::string name;
  int number_of_nodes = 0;
  int number_of_zones = 0;
  std::vector<LinkData> links;
  std::vector<std::vector<Real>> demand;
};

/// A link view retains its network storage, including after Python GC.
class Link {
 public:
  LinkData data() const;
  Real capacity() const;
  Real flow() const;
  void set_capacity(Real value);
  void set_flow(Real value);
  Real delay(Real flow = -1) const;

 private:
  friend class Network;
  Link(std::shared_ptr<detail::NetworkImpl> impl, std::size_t index);
  std::shared_ptr<detail::NetworkImpl> impl_;
  std::size_t index_;
};

/// Network state with private, stable-address native storage.
/// Mutable operations on the same network must be serialized by the caller.
class Network {
 public:
  explicit Network(const NetworkData& data);
  ~Network();
  Network(const Network&) = delete;
  Network& operator=(const Network&) = delete;

  std::string name() const;
  int number_of_nodes() const;
  int number_of_zones() const;
  int number_of_links() const;
  int number_of_od_pairs() const;
  Real total_travel_time() const;
  Real beckmann_objective() const;
  Real relative_gap() const;
  void reset();
  Link link(int index) const;
  std::vector<Real> flows() const;
  std::vector<Real> capacities() const;
  std::vector<Real> free_flow_times() const;
  std::vector<Real> link_costs() const;
  std::vector<std::array<int, 2>> link_nodes() const;
  void set_capacities(const std::vector<Real>& values);

 private:
  friend struct detail::Access;
  explicit Network(std::shared_ptr<detail::NetworkImpl> impl);
  std::shared_ptr<detail::NetworkImpl> impl_;
};

/// Load the project's preprocessed *_net.csv and *_trips.csv files.
std::shared_ptr<Network> LoadNetwork(const std::string& dataset,
                                     const std::string& data_root);
std::vector<LinkConstraint> LoadConstraints(const std::string& path,
                                            bool verbose = false);

}  // namespace traffic_assignment

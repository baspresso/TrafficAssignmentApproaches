#pragma once

#include <traffic_assignment/solve.hpp>
#include "tap/core/Network.h"
#include "tap/algorithms/TrafficAssignmentApproach.h"

namespace traffic_assignment::detail {

struct NetworkImpl {
  std::unique_ptr<TrafficAssignment::Network<Real>> network;
};

// Only implementation files can access algorithm templates and solver state.
struct Access {
  static TrafficAssignment::Network<Real>& Get(Network& network) {
    return *network.impl_->network;
  }
  static std::shared_ptr<Network> Wrap(std::shared_ptr<NetworkImpl> impl) {
    return std::shared_ptr<Network>(new Network(std::move(impl)));
  }
  static std::shared_ptr<TrafficAssignment::TrafficAssignmentApproach<Real>> Get(
      TrafficAssignmentApproach& approach);
};

}  // namespace traffic_assignment::detail

#include <traffic_assignment/solve.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include "detail/access.hpp"
#include "detail/input_validation.hpp"
#include "detail/cnd/BilevelCND.h"

namespace traffic_assignment {

struct BilevelCND::Impl {
  std::shared_ptr<Network> network;
  std::shared_ptr<TrafficAssignmentApproach> approach;
  std::size_t active_constraints = 0;
  std::unique_ptr<TrafficAssignment::BilevelCND<Real>> solver;
};

BilevelCND::BilevelCND(std::shared_ptr<Network> network,
                       std::shared_ptr<TrafficAssignmentApproach> approach,
                       const std::vector<LinkConstraint>& constraints,
                       const std::vector<StepConfig>& pipeline,
                       const CndpOptions& options)
    : impl_(std::make_unique<Impl>()) {
  if (!network || !approach) throw std::invalid_argument("network and approach must not be null");
  if (approach->network() != network) {
    throw std::invalid_argument("approach must belong to the supplied network");
  }
  if (constraints.size() != static_cast<std::size_t>(network->number_of_links())) {
    throw std::invalid_argument("constraints must have one entry per link");
  }
  const auto nodes = network->link_nodes();
  for (std::size_t i = 0; i < constraints.size(); ++i) {
    const auto& constraint = constraints[i];
    const auto field = "constraint[" + std::to_string(i) + "]";
    detail::ValidateConstraint(constraint, field);
    if (constraint.init_node != static_cast<std::size_t>(nodes[i][0]) ||
        constraint.term_node != static_cast<std::size_t>(nodes[i][1])) {
      throw std::invalid_argument(field + " endpoints must match the network link in native order");
    }
  }
  // Own the prepared constraints: callers' inputs are never modified.
  auto prepared = constraints;
  auto& native_network = detail::Access::Get(*network);
  if (options.exclude_insensitive_links) {
    constexpr Real epsilon = 1e-10L;
    for (std::size_t i = 0; i < prepared.size(); ++i) {
      const auto& link = native_network.links()[i];
      if (std::abs(link.b) < epsilon || std::abs(link.power) < epsilon ||
          std::abs(link.free_flow_time) < epsilon) {
        prepared[i].upper_bound = prepared[i].lower_bound;
      }
    }
  }
  impl_->active_constraints = static_cast<std::size_t>(std::count_if(
      prepared.begin(), prepared.end(), [](const auto& c) { return c.upper_bound > c.lower_bound; }));
  impl_->network = std::move(network);
  impl_->approach = std::move(approach);
  impl_->solver = std::make_unique<TrafficAssignment::BilevelCND<Real>>(
      native_network, detail::Access::Get(*impl_->approach), prepared, pipeline,
      options.link_capacity_selection_threshold, options.budget_threshold,
      options.budget_function_multiplier, options.budget_upper_bound, options.metrics,
      options.route_search_threads, options.progress_format);
  SetVerbose(options.verbose);
  SetStatisticsEnabled(options.statistics_enabled);
  SetFinalDiagnosticsEnabled(options.final_diagnostics);
}

BilevelCND::~BilevelCND() = default;
void BilevelCND::SetVerbose(bool enabled) { impl_->solver->SetVerbose(enabled); }
void BilevelCND::SetStatisticsEnabled(bool enabled) { impl_->solver->SetStatisticsEnabled(enabled); }
void BilevelCND::SetFinalDiagnosticsEnabled(bool enabled) {
  impl_->solver->SetFinalDiagnosticsEnabled(enabled);
}
void BilevelCND::SetRouteSearchThreadCount(std::size_t count) {
  impl_->solver->SetRouteSearchThreadCount(count);
}
std::size_t BilevelCND::active_constraint_count() const { return impl_->active_constraints; }
CndpResult BilevelCND::ComputeNetworkDesign() { return impl_->solver->ComputeNetworkDesign(); }

}  // namespace traffic_assignment

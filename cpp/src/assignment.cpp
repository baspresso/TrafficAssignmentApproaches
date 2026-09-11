#include <traffic_assignment/solve.hpp>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <stdexcept>
#include <utility>

#include "detail/access.hpp"
#include "detail/tap/algorithms/tapas/TapasApproach.h"
#include "detail/tap/algorithms/route_based/RouteBasedApproach.h"

namespace traffic_assignment {
namespace {
std::string ApproachKey(std::string value) {
  value.erase(std::remove(value.begin(), value.end(), '_'), value.end());
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value;
}

TapOptions TapasOptions(Real alpha, int max_iterations, Real mu, Real v) {
  TapOptions options;
  options.relative_gap_tolerance = alpha;
  options.max_iterations = max_iterations;
  options.mu = mu;
  options.v = v;
  return options;
}

TapOptions RouteOptions(Real alpha, const std::string& shift_method, std::size_t threads,
                         int max_iterations, int full, int origin, Real ema) {
  TapOptions options;
  options.approach = "routebased";
  options.relative_gap_tolerance = alpha;
  options.shift_method = shift_method;
  options.route_search_threads = threads;
  options.max_iterations = max_iterations;
  options.full_iteration_count = full;
  options.origin_iteration_count = origin;
  options.ema_alpha = ema;
  return options;
}
}

struct TrafficAssignmentApproach::Impl {
  // Destruction order: destroy the solver before releasing its network.
  std::shared_ptr<Network> network;
  std::shared_ptr<TrafficAssignment::TrafficAssignmentApproach<Real>> approach;
};

TrafficAssignmentApproach::TrafficAssignmentApproach(std::shared_ptr<Network> network,
                                                       const TapOptions& options)
    : impl_(std::make_unique<Impl>()) {
  if (!network) throw std::invalid_argument("network must not be null");
  impl_->network = std::move(network);
  auto& native_network = detail::Access::Get(*impl_->network);
  if (ApproachKey(options.approach) == "routebased") {
    impl_->approach = std::make_shared<TrafficAssignment::RouteBasedApproach<Real>>(
        native_network, options.relative_gap_tolerance, options.shift_method,
        options.route_search_threads, options.max_iterations, options.full_iteration_count,
        options.origin_iteration_count, options.ema_alpha);
  } else {
    impl_->approach = std::make_shared<TrafficAssignment::TapasApproach<Real>>(
        native_network, options.relative_gap_tolerance, options.max_iterations,
        options.mu, options.v);
  }
}

TrafficAssignmentApproach::~TrafficAssignmentApproach() = default;
std::shared_ptr<Network> TrafficAssignmentApproach::network() const { return impl_->network; }
std::string TrafficAssignmentApproach::GetApproachName() const {
  return impl_->approach->GetApproachName();
}
void TrafficAssignmentApproach::ComputeTrafficFlows(bool recording) {
  impl_->approach->ComputeTrafficFlows(recording);
}
void TrafficAssignmentApproach::Reset() { impl_->approach->Reset(); }
void TrafficAssignmentApproach::SetOutputRoot(const std::string& path) {
  impl_->approach->SetOutputRoot(path);
}
void TrafficAssignmentApproach::SetRouteSearchThreadCount(std::size_t count) {
  impl_->approach->SetRouteSearchThreadCount(count);
}
std::size_t TrafficAssignmentApproach::GetRouteSearchThreadCount() const {
  return impl_->approach->GetRouteSearchThreadCount();
}

std::shared_ptr<TrafficAssignment::TrafficAssignmentApproach<Real>> detail::Access::Get(
    TrafficAssignmentApproach& approach) {
  return approach.impl_->approach;
}

TapasApproach::TapasApproach(std::shared_ptr<Network> network, Real alpha,
                             int max_iterations, Real mu, Real v)
    : TrafficAssignmentApproach(std::move(network), TapasOptions(alpha, max_iterations, mu, v)) {}

RouteBasedApproach::RouteBasedApproach(std::shared_ptr<Network> network, Real alpha,
                                       const std::string& shift_method, std::size_t threads,
                                       int max_iterations, int full, int origin, Real ema)
    : TrafficAssignmentApproach(std::move(network),
                               RouteOptions(alpha, shift_method, threads, max_iterations,
                                            full, origin, ema)) {}

std::shared_ptr<TrafficAssignmentApproach> MakeApproach(std::shared_ptr<Network> network,
                                                       const TapOptions& options) {
  const auto key = ApproachKey(options.approach);
  if (key == "tapas" || key == "task" || key == "tasktapas") {
    return std::make_shared<TapasApproach>(std::move(network), options.relative_gap_tolerance,
                                          options.max_iterations, options.mu, options.v);
  }
  if (key == "routebased") {
    return std::make_shared<RouteBasedApproach>(std::move(network), options.relative_gap_tolerance,
        options.shift_method, options.route_search_threads, options.max_iterations,
        options.full_iteration_count, options.origin_iteration_count, options.ema_alpha);
  }
  throw std::invalid_argument("Unsupported approach '" + options.approach +
                              "'. Supported: Tapas, RouteBased.");
}

TapResult SolveTap(TrafficAssignmentApproach& approach, bool reset, bool recording) {
  if (reset) approach.Reset();
  const auto start = std::chrono::steady_clock::now();
  approach.ComputeTrafficFlows(recording);
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  const auto network = approach.network();
  return {approach.GetApproachName(), network->name(), network->relative_gap(),
          network->total_travel_time(), network->beckmann_objective(),
          network->flows(), network->link_costs(), seconds};
}

TapResult SolveTap(std::shared_ptr<Network> network, const TapOptions& options) {
  auto approach = MakeApproach(std::move(network), options);
  return SolveTap(*approach);
}

}  // namespace traffic_assignment

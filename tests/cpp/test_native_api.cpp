#include <traffic_assignment/solve.hpp>

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace ta = traffic_assignment;

void Require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

ta::NetworkData TwoRoutesData() {
  ta::NetworkData data;
  data.name = "TwoRoutes";
  data.number_of_nodes = 4;
  data.number_of_zones = 2;
  data.links = {{0, 1, 10, 0, 2, 1, 1}, {0, 2, 10, 0, 1, 1, 1},
                {2, 1, 10, 0, 1, 1, 1}, {2, 3, 10, 0, 1, 0, 1}};
  data.demand = {{0, 10}, {0, 0}};
  return data;
}

std::shared_ptr<ta::Network> TwoRoutes() {
  return std::make_shared<ta::Network>(TwoRoutesData());
}

template <typename Function>
void RequireInvalid(Function run, const char* message) {
  try {
    run();
  } catch (const std::invalid_argument&) {
    return;
  }
  throw std::runtime_error(message);
}

void CheckInputValidation() {
  for (auto field : {&ta::LinkData::capacity, &ta::LinkData::length,
                     &ta::LinkData::free_flow_time, &ta::LinkData::b,
                     &ta::LinkData::power, &ta::LinkData::speed, &ta::LinkData::toll}) {
    for (auto value : {-1.0L, std::numeric_limits<ta::Real>::quiet_NaN(),
                      std::numeric_limits<ta::Real>::infinity()}) {
      auto data = TwoRoutesData();
      data.links[0].*field = value;
      RequireInvalid([&] { ta::Network invalid(data); }, "invalid native link value accepted");
    }
  }
  auto data = TwoRoutesData();
  data.links[0].capacity = 0;
  RequireInvalid([&] { ta::Network invalid(data); }, "zero capacity accepted");
  data = TwoRoutesData();
  data.links[0].term = data.number_of_nodes;
  RequireInvalid([&] { ta::Network invalid(data); }, "out-of-range node accepted");
  data = TwoRoutesData();
  data.demand[0][1] = -1;
  RequireInvalid([&] { ta::Network invalid(data); }, "negative demand accepted");
  data.demand[0][1] = std::numeric_limits<ta::Real>::quiet_NaN();
  RequireInvalid([&] { ta::Network invalid(data); }, "nonfinite demand accepted");
  data = TwoRoutesData();
  data.demand[0].pop_back();
  RequireInvalid([&] { ta::Network invalid(data); }, "ragged demand accepted");

  auto network = TwoRoutes();
  auto capacities = network->capacities();
  capacities[0] = 15;
  capacities[1] = 0;
  RequireInvalid([&] { network->set_capacities(capacities); }, "invalid capacity update accepted");
  Require(network->capacities()[0] == 10, "failed capacity update mutated the network");
  RequireInvalid([&] { network->link(0).set_capacity(-1); }, "invalid link capacity accepted");

  std::vector<ta::LinkConstraint> constraints;
  for (const auto& nodes : network->link_nodes()) constraints.emplace_back(nodes[0], nodes[1], 8, 20, 1);
  auto approach = ta::MakeApproach(network);
  ta::StepConfig step;
  step.type = "nlopt";
  step.algorithm = "LN_COBYLA";
  for (auto field : {&ta::LinkConstraint::lower_bound, &ta::LinkConstraint::upper_bound,
                     &ta::LinkConstraint::investment_cost_param}) {
    auto invalid = constraints;
    // The final link is insensitive: validation must precede its bound filtering.
    invalid.back().*field = std::numeric_limits<double>::quiet_NaN();
    RequireInvalid([&] { ta::BilevelCND solver(network, approach, invalid, {step}); },
                   "invalid native constraints accepted");
    Require(network->capacities()[0] == 10, "invalid constraints mutated capacities");
  }
  std::swap(constraints[0], constraints[1]);
  RequireInvalid([&] { ta::BilevelCND solver(network, approach, constraints, {step}); },
                 "misaligned native constraints accepted");
}

int main() {
  try {
    CheckInputValidation();
    // This translation unit includes only public headers and links the core.
    auto network = TwoRoutes();
    auto approach = ta::MakeApproach(network);
    const auto result = ta::SolveTap(*approach);
    Require(std::abs(result.relative_gap) < 1e-12L, "TAP equilibrium gap");
    Require(std::abs(result.flows[0] - 5) < 1e-6L, "TAP route split");
    Require(std::abs(result.total_travel_time - 30) < 1e-6L, "TAP total travel time");

    auto link = network->link(0);
    network.reset();
    Require(std::abs(link.flow() - 5) < 1e-6L, "network retained by solver and link view");
    approach->Reset();
    Require(link.flow() == 0, "link view refers to live network state");
    Require(std::abs(result.flows[0] - 5) < 1e-6L, "result remains a snapshot");
    approach.reset();
    Require(link.capacity() == 10, "link view survives solver destruction");

    network = TwoRoutes();
    ta::TapOptions tap_options;
    tap_options.approach = "route_based";
    const auto route_result = ta::SolveTap(network, tap_options);
    Require(std::abs(route_result.total_travel_time - 30) < 1e-6L, "route-based public API");

    std::vector<ta::LinkConstraint> constraints;
    for (const auto& nodes : network->link_nodes()) {
      constraints.emplace_back(nodes[0], nodes[1], 10, 20, 1);
    }
    ta::StepConfig step;
    step.type = "nlopt";
    step.algorithm = "LN_COBYLA";
    step.max_iterations = 3;
    ta::CndpOptions options;
    options.budget_upper_bound = 100;
    approach = ta::MakeApproach(network);
    ta::BilevelCND solver(network, approach, constraints, {step}, options);
    Require(solver.active_constraint_count() == 3, "native insensitive-link filtering");
    const auto design = solver.ComputeNetworkDesign();
    Require(std::isfinite(design.objective), "finite native CNDP objective");
    Require(std::abs(design.objective - design.total_travel_time - design.budget) < 1e-9,
             "native CNDP objective decomposition");
    Require(design.upper_bounds[3] == 10, "native prepared bounds returned");
    Require(constraints[3].upper_bound == 20, "caller constraints unchanged");
    Require(design.budget <= options.budget_upper_bound + 1e-6, "native budget bound");

    bool rejected = false;
    try {
      ta::BilevelCND mismatch(TwoRoutes(), approach, constraints, {step}, options);
    } catch (const std::invalid_argument&) {
      rejected = true;
    }
    Require(rejected, "mismatched network and approach rejected");
    std::cout << "Native API checks passed\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}

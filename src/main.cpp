#include <iostream>
#include <iomanip>
#include <memory>
#include "../include/tap/algorithms/route_based/RouteBasedApproach.h"
#include "../include/tap/algorithms/demand_based/PumpOutDemandBasedApproach.h"
#include "../include/tap/algorithms/tapas_based/TapasApproach.h"
#include "../include/tap/core/NetworkBuilder.h"
#include "../include/cnd/BilevelCND.h"
#include "../include/cnd/DirectedConstraintLoader.h"
#include <nlopt.hpp>
#include <Eigen/Dense>
#include <boost/multiprecision/cpp_dec_float.hpp>

int main() {
	TrafficAssignment::NetworkBuilder builder;
	auto network = builder.BuildFromDataset<long double>("SiouxFalls_SingleODPair_0");

	auto approach = std::make_shared<TrafficAssignment::RouteBasedApproach<long double>>(
		network,
		1e-14,
		"NewtonStep"
	);

  TrafficAssignment::DirectedConstraintLoader loader;
  loader.SetVerbose(true);
  auto constraints = loader.LoadFromFile("C:/Projects/TrafficAssignmentApproaches/data/TransportationNetworks/SiouxFalls/SiouxFalls_constraints.csv");

  std::cout << "Creating BilevelCND solver..." << std::endl;
  TrafficAssignment::BilevelCND  <long double> cnd(network, approach, constraints, 0, 40);
  std::cout << "      BilevelCND solver created" << std::endl;
  std::cout << "      Design variables: " << constraints.size() << std::endl;

  std::cout << "\nRunning bilevel optimization..." << std::endl;
  std::cout << "      (Requires complete Network and TA implementation)"
            << std::endl;

  auto [optimal_y, min_obj] = cnd.ComputeNetworkDesign();

  std::cout << "\nOptimization Results:" << std::endl;
  for (auto now : optimal_y) {
    std::cout << now << ' ';
  }

	return 0;
}
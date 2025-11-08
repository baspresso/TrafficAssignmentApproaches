#include <iostream>
#include <iomanip>
#include <memory>
#include "../include/traffic_assignment/algorithms/route_based/RouteBasedApproach.h"
#include "../include/traffic_assignment/algorithms/demand_based/PumpOutDemandBasedApproach.h"
#include "../include/traffic_assignment/algorithms/tapas_based/TapasApproach.h"
#include "../include/traffic_assignment/core/NetworkBuilder.h"

int main() {
	TrafficAssignment::NetworkBuilder builder;
	auto network = builder.BuildFromDataset<long double>("SiouxFalls");

	// TrafficAssignment::TapasApproach<long double> approach(
	// 	network,
	// 	1e-14,
	// 	"LineSearch"
	// );

	TrafficAssignment::RouteBasedApproach<long double> approach(
		network,
		1e-14,
		"NewtonStep"
	);

	approach.ComputeTrafficFlows();

	return 0;
}
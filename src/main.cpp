#include <iostream>
#include <iomanip>
#include "../include/traffic_assignment/algorithms/route_based/RouteBasedKrylatov2023Approach.h"
#include "../include/traffic_assignment/algorithms/demand_based/PumpOutDemandBasedApproach.h"
#include "../include/traffic_assignment/algorithms/tapas_based/TapasApproach.h"
#include "../include/traffic_assignment/core/NetworkBuilder.h"
#include <memory>

#include "../include/traffic_assignment/algorithms/tapas_based/TapasShiftMethodFactory.h"

int main() {
	TrafficAssignment::NetworkBuilder builder;
	auto network = builder.BuildFromDataset<long double>("SiouxFalls"); // Returns shared_ptr

	TrafficAssignment::TapasApproach<long double> approach(
		network,
		1e-14,
		"LineSearch"
	);

	approach.ComputeTrafficFlows();

	return 0;
}
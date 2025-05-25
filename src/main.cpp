#include <iostream>
#include <iomanip>
#include "../include/traffic_assignment/algorithms/route_based/RouteBasedKrylatov2023Approach.h"
#include "../include/traffic_assignment/algorithms/demand_based/PumpOutDemandBasedApproach.h"
#include "../include/traffic_assignment/core/NetworkBuilder.h"
#include <memory>

int main() {
	// 1. Use NetworkBuilder to create a shared_ptr<Network>
	TrafficAssignment::NetworkBuilder builder;
	auto network = builder.BuildFromDataset<long double>("SiouxFalls"); // Returns shared_ptr

	// 2. Pass the shared_ptr directly (no std::move needed)
	TrafficAssignment::RouteBasedKrylatov2023Approach<long double> approach(
		network,  // Shared ownership
		1e-14     // Alpha value (optional)
	);

	// 3. Execute traffic assignment
	approach.ComputeTrafficFlows();
	//TrafficAssignment::RouteBasedKrylatov2023Approach<long double>("SiouxFalls").ComputeTrafficFlows();
	//TrafficAssignment::TapasLineSearchApproach<long double>("SiouxFalls").ComputeTrafficFlows();
	//TrafficAssignment::PumpOutDemandBasedApproach<long double>("SiouxFalls").ComputeTrafficFlows();
	//TrafficAssignment::TapasNewtonStepApproach<long double>("SiouxFalls").ComputeTrafficFlows();
	//TrafficAssignment::TapasAdvancedGradientDescentAppoach<long double>("SiouxFalls").ComputeTrafficFlows();
	//TrafficAssignment::RouteBasedKrylatov2023Approach<boost::multiprecision::cpp_bin_float_100>("Anaheim").ComputeTrafficFlows();
	//TrafficAssignment::TapasLineSearchApproach<boost::multiprecision::cpp_bin_float_100>("Anaheim").ComputeTrafficFlows();
	//TrafficAssignment::TapasAdvancedGradientDescentAppoach<boost::multiprecision::cpp_bin_float_100>("SiouxFalls").ComputeTrafficFlows();
	return 0;
}
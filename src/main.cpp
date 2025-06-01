#include <iostream>
#include <iomanip>
#include "../include/traffic_assignment/algorithms/route_based/RouteBasedKrylatov2023Approach.h"
#include "../include/traffic_assignment/algorithms/demand_based/PumpOutDemandBasedApproach.h"
#include "../include/traffic_assignment/algorithms/tapas_based/TapasNewtonStepApproach.h"
#include "../include/traffic_assignment/algorithms/tapas_based/TapasLineSearchApproach.h"
#include "../include/traffic_assignment/core/NetworkBuilder.h"
#include <memory>

int main() {
	TrafficAssignment::NetworkBuilder builder;
	auto network = builder.BuildFromDataset<long double>("SiouxFalls"); // Returns shared_ptr

	TrafficAssignment::PumpOutDemandBasedApproach<long double> approach(
		network,
		1e-14
	);

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
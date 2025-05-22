#include <iostream>
#include <iomanip>
#include "../include/traffic_assignment/algorithms/route_based/RouteBasedKrylatov2023Approach.h"
#include "../include/traffic_assignment/algorithms/demand_based/PumpOutDemandBasedApproach.h"
#include "../include/traffic_assignment/core/NetworkBuilder.h"


int main() {
	//TrafficAssignment::NetworkBuilder builder;
	//auto network = builder.BuildFromDataset<double>("SiouxFalls");
	//std::cout << 1 << '\n';
	//TrafficAssignment::RouteBasedKrylatov2023Approach<long double>("SiouxFalls").ComputeTrafficFlows();
	//TrafficAssignment::TapasLineSearchApproach<long double>("SiouxFalls").ComputeTrafficFlows();
	TrafficAssignment::PumpOutDemandBasedApproach<long double>("SiouxFalls").ComputeTrafficFlows();
	//TrafficAssignment::TapasNewtonStepApproach<long double>("SiouxFalls").ComputeTrafficFlows();
	//TrafficAssignment::TapasAdvancedGradientDescentAppoach<long double>("SiouxFalls").ComputeTrafficFlows();
	//TrafficAssignment::RouteBasedKrylatov2023Approach<boost::multiprecision::cpp_bin_float_100>("Anaheim").ComputeTrafficFlows();
	//TrafficAssignment::TapasLineSearchApproach<boost::multiprecision::cpp_bin_float_100>("Anaheim").ComputeTrafficFlows();
	//TrafficAssignment::TapasAdvancedGradientDescentAppoach<boost::multiprecision::cpp_bin_float_100>("SiouxFalls").ComputeTrafficFlows();
	return 0;
}
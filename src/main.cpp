#include <iostream>
#include <iomanip>
//#include "../include/RouteBasedKrylatov2023Approach.h"
#include "../include/traffic_assignment/algorithms/RouteBasedKrylatov2023Approach.h"
//#include "../include/TapasLineSearchApproach.h"
//#include "../include/TapasAdvancedGradientDescentApproach.h"
#include "../include/traffic_assignment/algorithms/TapasNewtonStepApproach.h"
//#include "Eigen/Dense"
//#include "boost/multiprecision/cpp_bin_float.hpp"
//#include "../include/PumpOutDemandBasedApproach.h"


int main() {
	//TrafficAssignment::RouteBasedKrylatov2023Approach<long double>("SiouxFalls").ComputeTrafficFlows();
	//TrafficAssignment::TapasLineSearchApproach<long double>("SiouxFalls").ComputeTrafficFlows();
	//TrafficAssignment::PumpOutDemandBasedApproach<long double>("SiouxFalls").ComputeTrafficFlows();
	TrafficAssignment::TapasNewtonStepApproach<long double>("SiouxFalls").ComputeTrafficFlows();
	//TrafficAssignment::TapasAdvancedGradientDescentAppoach<long double>("SiouxFalls").ComputeTrafficFlows();
	//TrafficAssignment::RouteBasedKrylatov2023Approach<boost::multiprecision::cpp_bin_float_100>("Anaheim").ComputeTrafficFlows();
	//TrafficAssignment::TapasLineSearchApproach<boost::multiprecision::cpp_bin_float_100>("Anaheim").ComputeTrafficFlows();
	//TrafficAssignment::TapasAdvancedGradientDescentAppoach<boost::multiprecision::cpp_bin_float_100>("SiouxFalls").ComputeTrafficFlows();
	return 0;
}
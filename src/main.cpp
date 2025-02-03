#include <iostream>
#include <iomanip>
#include "../include/RouteBasedKrylatov2023Approach.h"
#include "../include/Link.h"
#include "../include/TapasLineSearchApproach.h"
#include "../include/TapasAdvancedGradientDescentApproach.h"
#include "../include/TrafficAssignmentApproach.h"
#include "../include/TapasNewtonStepApproach.h"
//#include "Eigen/Dense"
//#include "boost/multiprecision/cpp_bin_float.hpp"
#include "../include/DemandBasedApproach.h"

int main() {
	//TrafficAssignment::RouteBasedKrylatov2023Approach<long double>("SiouxFalls").ComputeTrafficFlows();
	//TrafficAssignment::TapasLineSearchApproach<long double>("SiouxFalls").ComputeTrafficFlows();
	TrafficAssignment::DemandBasedApproach<long double>("SiouxFalls").ComputeTrafficFlows();
	//TrafficAssignment::TapasNewtonStepApproach<long double>("SiouxFalls").ComputeTrafficFlows();
	//TrafficAssignment::TapasAdvancedGradientDescentAppoach<long double>("SiouxFalls").ComputeTrafficFlows();
	//TrafficAssignment::RouteBasedKrylatov2023Approach<boost::multiprecision::cpp_bin_float_100>("Anaheim").ComputeTrafficFlows();
	//TrafficAssignment::TapasLineSearchApproach<boost::multiprecision::cpp_bin_float_100>("Anaheim").ComputeTrafficFlows();
	//TrafficAssignment::TapasAdvancedGradientDescentAppoach<boost::multiprecision::cpp_bin_float_100>("SiouxFalls").ComputeTrafficFlows();
    double a = 1e-10;
    long double b = 1e-10;
    std::cout <<  std::setprecision(20) << (a /= 5234) << '\n';
    //std::cout <<  std::setprecision(20) << (b /= 5234) << '\n';
	return 0;
}
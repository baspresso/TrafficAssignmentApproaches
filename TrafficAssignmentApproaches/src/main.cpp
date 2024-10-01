#include <iostream>
#include "RouteBasedKrylatov2023Approach.h"
#include "TapasLineSearchApproach.h"

int main() {
	TrafficAssignment::RouteBasedKrylatov2023Approach<long double>("Anaheim").ComputeTrafficFlows();
	TrafficAssignment::TapasLineSearchApproach<long double>("Anaheim").ComputeTrafficFlows();
	return 0;
}

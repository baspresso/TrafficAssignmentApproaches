#include <iostream>
#include "RouteBasedKrylatov2023Approach.h"

int main() {
	TrafficAssignment::RouteBasedKrylatov2023Approach<long double>("SiouxFalls").ComputeTrafficFlows();
	return 0;
}

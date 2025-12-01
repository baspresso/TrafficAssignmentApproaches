#include <iostream>
#include <iomanip>
#include <memory>
#include "../include/traffic_assignment/algorithms/route_based/RouteBasedApproach.h"
#include "../include/traffic_assignment/algorithms/demand_based/PumpOutDemandBasedApproach.h"
#include "../include/traffic_assignment/algorithms/tapas_based/TapasApproach.h"
#include "../include/traffic_assignment/core/NetworkBuilder.h"
#include "../include/cnd/BilevelCND.h"
#include "../include/cnd/DirectedConstraintLoader.h"
#include <nlopt.hpp>


// double my_black_box_function(const std::vector<double> &x, std::vector<double> &grad, void *data)
// {
//     // Check if gradient computation is requested.
//     // For derivative-free algorithms, grad will always be empty/NULL.
//     if (!grad.empty()) {
//         // You could theoretically compute a finite-difference gradient here,
//         // but for true black-box functions, you can simply leave this empty.
//         // NLopt will not call this for derivative-free algorithms.
//     }

//     // Calculate and return the function value.
//     double value = x[0] * (x[0] + 1) + (x[1] + 12) * x[1] + 5; // Example: f(x) = x1² + x2²
//     return value;
// }

int main() {
	std::cout << "yo" << '\n';

	// // Create an optimizer object with the Nelder-Mead algorithm in 2 dimensions.
    // nlopt::opt opt(nlopt::LN_NELDERMEAD, 2);

    // // Define and set the lower bounds (highly recommended).
    // std::vector<double> lb = {-10.0, -10.0}; // Lower bounds
    // std::vector<double> ub = {10.0, 10.0};   // Upper bounds
    // opt.set_lower_bounds(lb);
    // opt.set_upper_bounds(ub);

    // // Set the objective function.
    // opt.set_min_objective(my_black_box_function, nullptr);

    // // Set a stopping tolerance.
    // opt.set_xtol_rel(1e-4);

    // // Initial guess.
    // std::vector<double> x = {1.5, 1.5};
    // double minf;

    // try {
    //     // Run the optimization.
    //     nlopt::result result = opt.optimize(x, minf);

    //     std::cout << "Found minimum at: (" << x[0] << ", " << x[1] << ")" << std::endl;
    //     std::cout << "Minimum value: " << minf << std::endl;
    //     std::cout << "NLopt result code: " << result << std::endl;

    // } catch (std::exception &e) {
    //     std::cout << "NLopt failed: " << e.what() << std::endl;
    // }

	TrafficAssignment::NetworkBuilder builder;
	auto network = builder.BuildFromDataset<long double>("SiouxFalls");

	auto approach = std::make_shared<TrafficAssignment::RouteBasedApproach<long double>>(
		network,
		1e-14,
		"NewtonStep"
	);

    TrafficAssignment::DirectedConstraintLoader loader;
    loader.SetVerbose(true);
    auto constraints = loader.LoadFromFile("C:/Projects/TrafficAssignmentApproaches/data/TransportationNetworks/SiouxFalls/SiouxFalls_constraints.csv");

    //TrafficAssignment::BilevelCND <long double> cnd(network, approach, constraints);

	// Step 4: Create BilevelCND solver with constraints
      std::cout << "\n[4/5] Creating BilevelCND solver..." << std::endl;
      TrafficAssignment::BilevelCND  <long double> cnd(network, approach, constraints);
      std::cout << "      BilevelCND solver created" << std::endl;
      std::cout << "      Design variables: " << constraints.size() << std::endl;

    //    TrafficAssignment::Configure solver
    //    cnd.SetObjectiveWeight(0.7);  // 70% travel time, 30% investment
    //    cnd.SetOptimizationAlgorithm(nlopt::LD_SLSQP);
    //    cnd.SetOptimizationTolerance(1e-4);
    //    cnd.SetMaxIterations(500);
    //    cnd.SetVerbose(true);

    //   std::cout << "\n      Solver Configuration:" << std::endl;
    //   std::cout << "        Objective weight: 0.70 (70% TT, 30% cost)" << std::endl;
    //   std::cout << "        Algorithm: LD_SLSQP" << std::endl;
    //   std::cout << "        Tolerance: 1e-4" << std::endl;
    //   std::cout << "        Max iterations: 500" << std::endl;

      // Step 5: Run optimization
      std::cout << "\n[5/5] Running bilevel optimization..." << std::endl;
      std::cout << "      (Requires complete Network and TA implementation)"
                << std::endl;

    auto [optimal_y, min_obj] = cnd.ComputeNetworkDesign();

       std::cout << "\nOptimization Results:" << std::endl;
      std::cout << "  Minimum objective value: " << std::fixed
                << std::setprecision(6) << min_obj << std::endl;
      for (auto now : optimal_y) {
        std::cout << now << ' ';
      }

	return 0;
}
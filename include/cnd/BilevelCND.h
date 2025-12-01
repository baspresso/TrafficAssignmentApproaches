#ifndef BILEVEL_CND_H
#define BILEVEL_CND_H

#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <memory>
#include <functional>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "../traffic_assignment/utils/DataProcessor.h"
#include "../traffic_assignment/algorithms/common/TrafficAssignmentApproach.h"
#include "../traffic_assignment/core/Network.h"
#include "DirectedConstraintLoader.h"
#include <nlopt.hpp>

namespace TrafficAssignment {

template <typename T>
class BilevelCND {
public:
  /**
   * @brief Constructor with constraints
   *
   * @param network Reference to the transportation network
   * @param approach Shared pointer to traffic assignment approach (lower level solver)
   * @param constraints Vector of directed link capacity constraints
   */
  BilevelCND(
    Network<T>& network,
    std::shared_ptr<TrafficAssignmentApproach<T>> approach,
    const std::vector<DirectedLinkCapacityConstraint>& constraints)
    : network_(network),
      approach_(approach),
      constraints_(constraints),
      objective_weight_(0.5),
      optimization_algorithm_(nlopt::LN_COBYLA),
      optimization_tolerance_(1e-4),
      max_iterations_(100),
      verbose_(true),
      objective_value_(std::numeric_limits<double>::max()),
      total_travel_time_(0.0),
      total_investment_cost_(0.0) {}

  /**
   * @brief Destructor
   */
  ~BilevelCND() = default;

  /**
   * @brief Compute the optimal network design
   * @return std::pair with optimal capacities and minimum objective value
   */
  std::pair<std::vector<double>, T> ComputeNetworkDesign() {

    if (verbose_) {
      std::cout << "\n" << std::string(70, '=') << std::endl;
      std::cout << "BILEVEL CONTINUOUS NETWORK DESIGN PROBLEM" << std::endl;
      std::cout << std::string(70, '=') << std::endl;
      std::cout << "\nProblem Configuration:" << std::endl;
      std::cout << "  Design variables: " << constraints_.size() << std::endl;
      std::cout << "  Objective weight (α): " << std::fixed << std::setprecision(4)
                << objective_weight_ << std::endl;
      std::cout << "  Algorithm: " << GetAlgorithmName(optimization_algorithm_)
                << std::endl;
      std::cout << "  Tolerance: " << optimization_tolerance_ << std::endl;
      std::cout << "  Max iterations: " << max_iterations_ << std::endl;
    }

    // Initialize initial guess if not provided
    // if (initial_guess_.empty()) {
    //   initial_guess_.resize(constraints_.size());
    //   for (std::size_t i = 0; i < constraints_.size(); ++i) {
    //     initial_guess_[i] = (constraints_[i].lower_bound + constraints_[i].upper_bound) / 2;
    //   }
    //   if (verbose_) {
    //     std::cout << "\n  Using lower bounds as initial guess" << std::endl;
    //   }
    // }

    initial_guess_.resize(constraints_.size());
      for (std::size_t i = 0; i < constraints_.size(); ++i) {
        initial_guess_[i] = (constraints_[i].lower_bound + constraints_[i].upper_bound) / 2;
      }

    int n_vars = static_cast<int>(constraints_.size());
    nlopt::opt optimizer(optimization_algorithm_, n_vars);

    optimizer.set_min_objective(&BilevelCND<T>::NLOptObjectiveFunctionWrapper, this);

    std::vector<double> lower_bounds(n_vars);
    std::vector<double> upper_bounds(n_vars);

    for (int i = 0; i < n_vars; ++i) {
      lower_bounds[i] = constraints_[i].lower_bound;
      upper_bounds[i] = constraints_[i].upper_bound; 
    }

    optimizer.set_maxeval(50);

    optimizer.set_lower_bounds(lower_bounds);
    optimizer.set_upper_bounds(upper_bounds);
    
    std::vector<double> x(initial_guess_.begin(), initial_guess_.end());
    double minf;
    nlopt::result result = optimizer.optimize(x, minf);

    return {x, network_.TotalTravelTime()};
  }

private:
  Network<T>& network_;
  std::shared_ptr<TrafficAssignmentApproach<T>> approach_;
  std::vector<DirectedLinkCapacityConstraint> constraints_;
  //std::map<std::pair<std::size_t, std::size_t>, std::size_t> arc_to_link_map_;

  // Optimization parameters
  T objective_weight_;
  nlopt::algorithm optimization_algorithm_;
  T optimization_tolerance_;
  int max_iterations_;
  bool verbose_;
  std::vector<T> initial_guess_;

  // Results storage
  T objective_value_;
  T total_travel_time_;
  T total_investment_cost_;


  private:
  class ObjectiveFunctor {
  public:
    ObjectiveFunctor(BilevelCND<T>* instance) : instance_(instance) {}
    
    double operator()(const std::vector<double>& x,
                      std::vector<double>& grad,
                      void* data) {
      return instance_->NLOptObjectiveFunction(x, grad);
    }
    
  private:
    BilevelCND<T>* instance_;
  };

  static double NLOptObjectiveFunctionWrapper(
    unsigned n,
    const double* x,
    double* grad,
    void* data)
{
    auto* instance = static_cast<BilevelCND<T>*>(data);
    return instance->NLOptObjectiveFunction(n, x, grad);
}

  double NLOptObjectiveFunction(
    unsigned n,
    const double* x,
    double* grad)
{
    for (unsigned i = 0; i < n; i++) {
        network_.mutable_links()[i].capacity = x[i];
    }

    approach_->ComputeTrafficFlows();

    double v = network_.TotalTravelTime();
    std::cout << v << "\n";
    return v;
}

  /**
   * @brief Get human-readable algorithm name
   */
  static std::string GetAlgorithmName(nlopt::algorithm algo) {
    switch (algo) {
      case nlopt::LD_SLSQP:
        return "LD_SLSQP (Sequential Least-Squares Programming)";
      case nlopt::LD_MMA:
        return "LD_MMA (Method of Moving Asymptotes)";
      case nlopt::LN_COBYLA:
        return "LN_COBYLA (Constrained Optimization BY Linear Approximations)";
      case nlopt::GN_ISRES:
        return "GN_ISRES (Improved Stochastic Ranking Evolution Strategy)";
      default:
        return "Custom Algorithm (" + std::to_string(algo) + ")";
    }
  }

  // Prevent copying
  BilevelCND(const BilevelCND&) = delete;
  BilevelCND& operator=(const BilevelCND&) = delete;
};

}  // namespace TrafficAssignment

#endif  // BILEVEL_CND_H
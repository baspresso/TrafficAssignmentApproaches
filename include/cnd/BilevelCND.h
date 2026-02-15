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
#include <chrono>
#include <sstream>
#include "../tap/utils/DataProcessor.h"
#include "../tap/algorithms/common/TrafficAssignmentApproach.h"
#include "../tap/core/Network.h"
#include "DirectedConstraintLoader.h"
#include <nlopt.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

namespace mp = boost::multiprecision;
using high_prec = mp::cpp_dec_float_50;

namespace TrafficAssignment {

template <typename T>
class BilevelCND {
  using MatrixXd = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
  using MatrixHightPrecision = Eigen::Matrix<high_prec, Eigen::Dynamic, Eigen::Dynamic>;
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
    const std::vector<DirectedLinkCapacityConstraint>& constraints,
    int max_standard_iterations = 1000,
    int max_optimality_condition_iterations = 40,
    T optimization_tolerance = 1e-4,
    T link_capacity_selection_threshold = 1e-3,
    T budget_threshold = 1e-1,
    T budget_function_multiplier = 5,
    double budget_upper_bound = 100000,
    nlopt::algorithm optimization_algorithm = nlopt::LN_COBYLA
  )
    : network_(network),
      approach_(approach),
      constraints_(constraints),
      max_standard_iterations_(max_standard_iterations),
      max_optimality_condition_iterations_(max_optimality_condition_iterations),
      link_capacity_selection_threshold_(link_capacity_selection_threshold),
      budget_threshold_(budget_threshold),
      optimization_tolerance_(optimization_tolerance),
      budget_function_multiplier_(budget_function_multiplier),
      budget_upper_bound_(budget_upper_bound),
      optimization_algorithm_(optimization_algorithm),
      verbose_(true),
      objective_value_(std::numeric_limits<double>::max()),
      total_travel_time_(0.0),
      total_investment_cost_(0.0),
      standard_objective_eval_count_(0),
      standard_total_eval_target_(0),
      standard_progress_active_(false),
      optimality_condition_iteration_count_(0),
      optimality_condition_total_iterations_(0),
      optimality_condition_progress_active_(false) {}

  /**
   * @brief Destructor
   */
  ~BilevelCND() = default;

  /**
   * @brief Compute the optimal network design
   * @return std::pair with optimal capacities and minimum objective value
   */
  void ComputeNetworkDesign() {
    if (verbose_) {
      std::cout << "\n" << std::string(70, '=') << std::endl;
      std::cout << "BILEVEL CONTINUOUS NETWORK DESIGN PROBLEM" << std::endl;
      std::cout << std::string(70, '=') << std::endl;
      std::cout << "\nProblem Configuration:" << std::endl;
      std::cout << "  Design variables: " << constraints_.size() << std::endl;
      std::cout << "  Algorithm: " << GetAlgorithmName(optimization_algorithm_)
                << std::endl;
      std::cout << "  Tolerance: " << optimization_tolerance_ << std::endl;
      std::cout << "  Max standard iterations: " << max_standard_iterations_ << std::endl;
      std::cout << "  Max optimality condition iterations: " << max_optimality_condition_iterations_ << std::endl;
    }
    const auto optimization_start_time = std::chrono::steady_clock::now();
    initial_guess_.resize(constraints_.size());
    for (std::size_t i = 0; i < constraints_.size(); ++i) {
      initial_guess_[i] = constraints_[i].lower_bound;
      network_.mutable_links()[i].capacity = initial_guess_[i];
    }

    std::cout << "StandardOptimization-start\n";
    StandardOptimization();
    approach_->Reset();
    approach_->ComputeTrafficFlows();

    std::cout << "OptimalityCondition-start\n";

    OptimalityConditionOptimization();

    //std::cout << "budget" << BudgetFunction(x.size(), x.data()) << '\n';

    //OptimalityCondition();
    approach_->Reset();
    approach_->ComputeTrafficFlows();

    std::cout << Budget() << '\n';
    const auto optimization_end_time = std::chrono::steady_clock::now();
    const double optimization_elapsed_seconds =
      std::chrono::duration<double>(optimization_end_time - optimization_start_time).count();
    std::ostringstream timing_line;
    double total_travel_time = network_.TotalTravelTime();
    double budget_function = static_cast<double>(Budget());
    double objective_function = total_travel_time + budget_function;
    timing_line << "optimization_time = "
                << std::fixed << std::setprecision(2)
                << optimization_elapsed_seconds << "s " << std::setprecision(10)
                << "objective_function=" << objective_function
                << " total_travel_time=" << total_travel_time
                << " budget_function=" << budget_function << "\n";
    std::cout << timing_line.str();

    std::cout << "RecordStatistics-Start\n";

    RecordStatistics();

    std::cout << "RecordStatistics-End\n";
  }

private:
  Network<T>& network_;
  std::shared_ptr<TrafficAssignmentApproach<T>> approach_;
  std::vector<DirectedLinkCapacityConstraint> constraints_;

  nlopt::algorithm optimization_algorithm_;
  T optimization_tolerance_;
  int max_standard_iterations_;
  int max_optimality_condition_iterations_;
  T link_capacity_selection_threshold_;
  T budget_threshold_;
  T budget_function_multiplier_;
  double budget_upper_bound_;
  bool verbose_;
  std::vector<T> initial_guess_;

  // Results storage
  T objective_value_;
  T total_travel_time_;
  T total_investment_cost_;
  int standard_objective_eval_count_;
  int standard_total_eval_target_;
  std::chrono::steady_clock::time_point standard_optimization_start_time_;
  bool standard_progress_active_;
  int optimality_condition_iteration_count_;
  int optimality_condition_total_iterations_;
  std::chrono::steady_clock::time_point optimality_condition_start_time_;
  bool optimality_condition_progress_active_;

  void StartStandardOptimizationProgress() {
    standard_objective_eval_count_ = 0;
    standard_total_eval_target_ = max_standard_iterations_;
    standard_optimization_start_time_ = std::chrono::steady_clock::now();
    standard_progress_active_ = standard_total_eval_target_ > 0;
  }

  void FinishStandardOptimizationProgress() {
    if (!standard_progress_active_) {
      return;
    }
    std::cout << '\n';
    standard_progress_active_ = false;
  }

  void PrintStandardOptimizationProgress(double objective_function,
                                         double total_travel_time,
                                         double budget_function) {
    if (!standard_progress_active_) {
      return;
    }

    ++standard_objective_eval_count_;
    int total = std::max(1, standard_total_eval_target_);
    int current = std::min(standard_objective_eval_count_, total);
    int left = std::max(0, total - current);

    constexpr int kBarWidth = 30;
    double ratio = static_cast<double>(current) / static_cast<double>(total);
    int filled = static_cast<int>(ratio * kBarWidth);

    std::string bar;
    bar.reserve(kBarWidth);
    for (int i = 0; i < kBarWidth; ++i) {
      if (i < filled) {
        bar.push_back('=');
      } else if (i == filled && filled < kBarWidth) {
        bar.push_back('>');
      } else {
        bar.push_back(' ');
      }
    }

    auto now = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(now - standard_optimization_start_time_).count();
    double it_per_sec = elapsed > 0.0 ? static_cast<double>(current) / elapsed : 0.0;

    std::ostringstream progress_line;
    progress_line << '\r'
                  << std::setw(3) << static_cast<int>(ratio * 100.0) << "%|"
                  << bar << "| "
                  << current << "/" << total
                  << " left:" << left
                  << " [" << std::fixed << std::setprecision(2) << it_per_sec << "it/s] "
                  << std::defaultfloat
                  << "objective_function=" << objective_function
                  << " total_travel_time=" << total_travel_time
                  << " budget_function=" << budget_function;
    std::cout << progress_line.str() << std::flush;
  }

  void StartOptimalityConditionProgress() {
    optimality_condition_iteration_count_ = 0;
    optimality_condition_total_iterations_ = max_optimality_condition_iterations_;
    optimality_condition_start_time_ = std::chrono::steady_clock::now();
    optimality_condition_progress_active_ = optimality_condition_total_iterations_ > 0;
  }

  void FinishOptimalityConditionProgress() {
    if (!optimality_condition_progress_active_) {
      return;
    }
    std::cout << '\n';
    optimality_condition_progress_active_ = false;
  }

  void PrintOptimalityConditionProgress(double objective_function,
                                        double total_travel_time,
                                        double budget_function) {
    if (!optimality_condition_progress_active_) {
      return;
    }

    ++optimality_condition_iteration_count_;
    int total = std::max(1, optimality_condition_total_iterations_);
    int current = std::min(optimality_condition_iteration_count_, total);
    int left = std::max(0, total - current);

    constexpr int kBarWidth = 30;
    double ratio = static_cast<double>(current) / static_cast<double>(total);
    int filled = static_cast<int>(ratio * kBarWidth);

    std::string bar;
    bar.reserve(kBarWidth);
    for (int i = 0; i < kBarWidth; ++i) {
      if (i < filled) {
        bar.push_back('=');
      } else if (i == filled && filled < kBarWidth) {
        bar.push_back('>');
      } else {
        bar.push_back(' ');
      }
    }

    auto now = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(now - optimality_condition_start_time_).count();
    double it_per_sec = elapsed > 0.0 ? static_cast<double>(current) / elapsed : 0.0;

    std::ostringstream progress_line;
    progress_line << '\r'
                  << std::setw(3) << static_cast<int>(ratio * 100.0) << "%|"
                  << bar << "| "
                  << current << "/" << total
                  << " left:" << left
                  << " [" << std::fixed << std::setprecision(2) << it_per_sec << "it/s] "
                  << std::defaultfloat
                  << "objective_function=" << objective_function
                  << " total_travel_time=" << total_travel_time
                  << " budget_function=" << budget_function;
    std::cout << progress_line.str() << std::flush;
  }

  void EnforceBudgetFeasibility(std::vector<double>& x,
                                const std::vector<double>& lower_bounds) {
    const int n = static_cast<int>(x.size());
    double budget = BudgetFunction(n, x.data());
    if (budget <= budget_upper_bound_) {
      return;
    }

    const double multiplier = static_cast<double>(budget_function_multiplier_);
    if (multiplier <= 0.0) {
      return;
    }

    double excess = (budget - budget_upper_bound_) / multiplier;
    if (excess <= 0.0) {
      return;
    }

    double total_slack = 0.0;
    for (int i = 0; i < n; ++i) {
      total_slack += std::max(0.0, x[i] - lower_bounds[i]);
    }
    if (total_slack <= 0.0) {
      return;
    }

    const double keep_ratio = std::max(0.0, (total_slack - excess) / total_slack);
    for (int i = 0; i < n; ++i) {
      const double slack = std::max(0.0, x[i] - lower_bounds[i]);
      x[i] = lower_bounds[i] + slack * keep_ratio;
    }

    budget = BudgetFunction(n, x.data());
    if (budget <= budget_upper_bound_) {
      return;
    }

    // Final numeric cleanup: trim residual violation directly on slack variables.
    excess = (budget - budget_upper_bound_) / multiplier;
    for (int i = 0; i < n && excess > 0.0; ++i) {
      const double reducible = std::max(0.0, x[i] - lower_bounds[i]);
      const double reduce_by = std::min(reducible, excess);
      x[i] -= reduce_by;
      excess -= reduce_by;
    }
  }


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
    void* data) {
    auto* instance = static_cast<BilevelCND<T>*>(data);
    return instance->NLOptObjectiveFunction(n, x, grad);
  }

  T Budget() {
    T result = 0;
    for (int link_index = 0; link_index < network_.number_of_links(); link_index++) {
      result +=  network_.mutable_links()[link_index].capacity - constraints_[link_index].lower_bound;
    }
    return result * budget_function_multiplier_;
  }

  double BudgetFunction(int n, const double* x) {
    double result = 0;
    for (int i = 0; i < n; i++) {
      result += x[i] - constraints_[i].lower_bound;
    }
    return result * budget_function_multiplier_;
  }

  static double BudgetConstraintWrapper(
    unsigned n,
    const double* x,
    double* grad,
    void* data) {
      auto* instance = static_cast<BilevelCND<T>*>(data);
      return instance->BudgetConstraint(n, x, grad);
  }
  
  double BudgetConstraint(unsigned n, const double* x, double* grad) {
    if (grad) {
      for (unsigned i = 0; i < n; i++) {
        grad[i] = budget_function_multiplier_;
      }
    }
    double budget = BudgetFunction(n, x);
    return budget - budget_upper_bound_;
  }

  double NLOptObjectiveFunction(
    unsigned n,
    const double* x,
    double* grad) {
    for (unsigned i = 0; i < n; i++) {
        network_.mutable_links()[i].capacity = x[i];
    }

    approach_->Reset();

    approach_->ComputeTrafficFlows();

    double total_travel_time = network_.TotalTravelTime();
    double budget_function = BudgetFunction(n, x);
    double objective_function = total_travel_time + budget_function;
    PrintStandardOptimizationProgress(objective_function, total_travel_time, budget_function);
    return objective_function;
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

  void writeVectorsToCSV(const std::vector<std::vector<T>>& vectors, 
                      const std::vector<std::string>& headers) {
    std::filesystem::path current_path = "C:/Projects/TrafficAssignmentApproaches";
    auto file_path = current_path / "data" / "TransportationNetworks" / network_.name() / (network_.name() + "_optimality.csv");
    std::ofstream file(file_path);
    // Записываем заголовки
    for (size_t i = 0; i < headers.size(); ++i) {
      file << headers[i];
      if (i != headers.size() - 1) file << ",";
    }
    file << "\n";
    
    // Записываем данные построчно
    size_t rows = vectors[0].size();
    for (size_t row = 0; row < rows; ++row) {
      for (size_t col = 0; col < vectors.size(); ++col) {
        file << vectors[col][row];
        if (col != vectors.size() - 1) file << ",";
      }
      file << "\n";
    }
    
    file.close();
}

  void RecordStatistics() {
    std::vector <T> link_condition_result = OptimalityCondition();
    std::vector <T> flow(network_.number_of_links(), 0);
    std::vector <T> lower_bound(network_.number_of_links(), 0);
    std::vector <T> upper_bound(network_.number_of_links(), 0);
    std::vector <T> capacity(network_.number_of_links(), 0);

    for (int i = 0; i < network_.number_of_links(); i++) {
      flow[i] = network_.mutable_links()[i].flow;
      lower_bound[i] = constraints_[i].lower_bound;
      upper_bound[i] = constraints_[i].upper_bound;
      capacity[i] = network_.mutable_links()[i].capacity;
    }
    std::vector <std::vector <T>> vectors = 
    {
      flow,
      lower_bound,
      upper_bound,
      capacity,
      link_condition_result
    };
    std::vector <std::string> headers =
    {
      "flow",
      "lower_bound",
      "upper_bound",
      "capacity",
      "condition_result"
    };
    writeVectorsToCSV(vectors, headers);
  }

  MatrixXd CapacityDerColumn(std::vector <std::vector <int>> routes, int link_index) {
    MatrixXd res = MatrixXd::Constant(routes.size(), 1, 0.0);
    for (int i = 0; i < routes.size(); i++) {
      if (std::find(routes[i].begin(), routes[i].end(), link_index) != routes[i].end())  {
        res(i, 0) = network_.links()[link_index].DelayCapacityDer();
      }
    }
    return res;
  }

  MatrixHightPrecision
  ConvertEigenMatrix(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& source) {
    MatrixHightPrecision target(
      source.rows(), source.cols());
    for (int i = 0; i < source.rows(); ++i) {
      for (int j = 0; j < source.cols(); ++j) {
        target(i, j) = static_cast<high_prec>(source(i, j));
      }
    }
    return target;
  }
  
  T OptimalityFunction(int link_index, bool flow_computatation = true) {
    if (flow_computatation) {
      approach_->Reset();
      approach_->ComputeTrafficFlows();
    }
    T result = 0;
    for (int od_pair_index = 0; od_pair_index < network_.number_of_od_pairs(); od_pair_index++) {
      const int routes_count = network_.mutable_od_pairs()[od_pair_index].GetRoutesCount();
      auto routes = network_.mutable_od_pairs()[od_pair_index].GetRoutes();
      auto e_column = MatrixHightPrecision::Constant(routes_count, 1, 1.0);
      auto jacobi = ConvertEigenMatrix(RoutesJacobiMatrix(routes, network_.mutable_links()));
      auto capacity_der_column = ConvertEigenMatrix(CapacityDerColumn(routes, link_index));
      result += network_.od_pairs()[od_pair_index].GetDemand() * 
                                            static_cast<T>(- (e_column.transpose() * jacobi.inverse() * capacity_der_column)(0, 0) / 
                                            (e_column.transpose() * jacobi.inverse() * e_column)(0, 0));
    }
    result /= budget_function_multiplier_;
    return result;
  }

  std::vector <int> MiddleBoundTypeIndeces() {
    std::vector <int> middle_bound_type_indeces;
    T bound_treshold = 1e-5;

    for (int link_index = 0; link_index < network_.number_of_links(); link_index++) {
      auto cur_capacity = network_.mutable_links()[link_index].capacity;
      if (std::abs(constraints_[link_index].upper_bound - cur_capacity) > bound_treshold &&
          std::abs(constraints_[link_index].lower_bound - cur_capacity) > bound_treshold) {
        middle_bound_type_indeces.push_back(link_index);
      }
    }
    return middle_bound_type_indeces;
  }

  bool TransferResult(T amount, int max_optimality_index, int min_optimality_index) {
    network_.mutable_links()[max_optimality_index].capacity += amount;
    network_.mutable_links()[min_optimality_index].capacity -= amount;
    bool result = OptimalityFunction(max_optimality_index) > OptimalityFunction(min_optimality_index);
    network_.mutable_links()[max_optimality_index].capacity -= amount;
    network_.mutable_links()[min_optimality_index].capacity += amount;
    return result;
  }

  T RespectedAmountConstranints(T amount, int max_optimality_index, int min_optimality_index) {
    auto max_optimality_capacity = network_.mutable_links()[max_optimality_index].capacity;
    auto min_optimality_capacity = network_.mutable_links()[min_optimality_index].capacity;
    amount = std::min(amount, constraints_[max_optimality_index].upper_bound - max_optimality_capacity);
    amount = std::min(amount, min_optimality_capacity - constraints_[min_optimality_index].lower_bound);
    return amount;
  }

  T OptimalityConditionTransferAmount(int max_optimality_index, int min_optimality_index) {
    // auto max_optimality_capacity = network_.mutable_links()[max_index].capacity;
    // auto min_optimality_capacity = network_.mutable_links()[min_index].capacity;

    T amount = 1;
    amount = RespectedAmountConstranints(amount, max_optimality_index, min_optimality_index);
    T next_amount = RespectedAmountConstranints(amount * 2, max_optimality_index, min_optimality_index);

    while (TransferResult(next_amount, max_optimality_index, min_optimality_index)) {
      amount = next_amount;
      next_amount = RespectedAmountConstranints(amount * 2, max_optimality_index, min_optimality_index);
      if (std::abs(next_amount - amount) < link_capacity_selection_threshold_ / 10) {
        break;
      }
    }
    while (!TransferResult(amount, max_optimality_index, min_optimality_index)) {
      amount /= 2;
    }
    return amount;
  }

  void StandardOptimization() {
    if (max_standard_iterations_ == 0) {
      return;
    }
    initial_guess_.resize(constraints_.size());
    for (std::size_t i = 0; i < constraints_.size(); ++i) {
      initial_guess_[i] = constraints_[i].lower_bound;
      network_.mutable_links()[i].capacity = initial_guess_[i];
    }
    std::vector<double> x(initial_guess_.begin(), initial_guess_.end());
    int n_vars = static_cast<int>(constraints_.size());
    nlopt::opt optimizer(optimization_algorithm_, n_vars);

    optimizer.set_min_objective(&BilevelCND<T>::NLOptObjectiveFunctionWrapper, this);

    std::vector<double> lower_bounds(n_vars);
    std::vector<double> upper_bounds(n_vars);

    for (int i = 0; i < n_vars; ++i) {
      lower_bounds[i] = constraints_[i].lower_bound;
      upper_bounds[i] = constraints_[i].upper_bound; 
    }

    optimizer.set_maxeval(max_standard_iterations_);

    optimizer.set_lower_bounds(lower_bounds);
    optimizer.set_upper_bounds(upper_bounds);

    double constraint_tol = 1e-4;
    
    optimizer.add_inequality_constraint(
      &BilevelCND<T>::BudgetConstraintWrapper,
      this,
      constraint_tol
    );

    double minf;
    
    StartStandardOptimizationProgress();
    try {
      optimizer.optimize(x, minf);
    } catch (...) {
      FinishStandardOptimizationProgress();
      throw;
    }
    FinishStandardOptimizationProgress();

    EnforceBudgetFeasibility(x, lower_bounds);

    for (unsigned i = 0; i < network_.number_of_links(); i++) {
      network_.mutable_links()[i].capacity = x[i];
    }
  }

  void OptimalityConditionOptimization() {
    if (max_optimality_condition_iterations_ == 0) {
      return;
    }

    StartOptimalityConditionProgress();
    try {
      for (int cnt = 0; cnt < max_optimality_condition_iterations_; cnt++) {
        (void)cnt;
        if (Budget() < budget_upper_bound_ - budget_threshold_) {
          NonUpperBoundBudgetOptimalityIteration();
        }
        else {
          UpperBoundBudgetOptimalityIteration();
        }

        approach_->Reset();
        approach_->ComputeTrafficFlows();
        double total_travel_time = network_.TotalTravelTime();
        double budget_function = static_cast<double>(Budget());
        double objective_function = total_travel_time + budget_function;
        PrintOptimalityConditionProgress(
          objective_function, total_travel_time, budget_function);
      }
    } catch (...) {
      FinishOptimalityConditionProgress();
      throw;
    }
    FinishOptimalityConditionProgress();
  }

  std::vector <T> OptimalityCondition() {
    //std::cout << std::fixed << std::setprecision(20);
    std::vector <T> link_condition_result(network_.number_of_links());
    for (int link_index = 0; link_index < network_.number_of_links(); link_index++) {
      link_condition_result[link_index] = OptimalityFunction(link_index);
      //std::cout << link_condition_result[link_index] << ' ';
    } 

    //std::cout << std::fixed << std::setprecision(10);
    return link_condition_result;
  }

  void UpperBoundBudgetOptimalityIteration() {
    int max_index = -1, min_index = -1;
    
    for (int link_index = 0; link_index < network_.number_of_links(); link_index++) {
      if (network_.mutable_links()[link_index].flow < 1) {
        continue;
      }
      auto cur_capacity = network_.mutable_links()[link_index].capacity;
      if (std::abs(constraints_[link_index].upper_bound - cur_capacity) > link_capacity_selection_threshold_ &&
          std::abs(constraints_[link_index].lower_bound - cur_capacity) > link_capacity_selection_threshold_) {
        max_index = link_index;
        min_index = link_index;
        break;
      }
    }
    if (max_index == -1) {
      return;
    }

    approach_->Reset();
    approach_->ComputeTrafficFlows();

    T max_value = OptimalityFunction(max_index, false);
    T min_value = max_value; 

    for (int link_index = 0; link_index < network_.number_of_links(); link_index++) {
      if (network_.mutable_links()[link_index].flow < 1) {
        continue;
      }
      auto cur_capacity = network_.mutable_links()[link_index].capacity;
      auto cur_value = OptimalityFunction(link_index, false);
      if (cur_value > max_value &&
          std::abs(constraints_[link_index].upper_bound - cur_capacity) > link_capacity_selection_threshold_) {
        max_index = link_index;
        max_value = cur_value;
      }
      if (cur_value < min_value &&
          std::abs(constraints_[link_index].lower_bound - cur_capacity) > link_capacity_selection_threshold_) {
        min_index = link_index;
        min_value = cur_value;
      }
    }

    T amount = OptimalityConditionTransferAmount(max_index, min_index);
    if (!optimality_condition_progress_active_) {
      std::cout << "max_index " << max_index << ' ' << "min_index " << min_index << " amount " << amount << '\n';
    }

    network_.mutable_links()[max_index].capacity += amount;
    network_.mutable_links()[min_index].capacity -= amount;
  }

  bool ProcessingCheck(int link_index, T middle_bound_value) {
    if (network_.mutable_links()[link_index].flow < 1e-3) {
      return false;
    }
    auto cur_capacity = network_.mutable_links()[link_index].capacity;
    auto cur_value = OptimalityFunction(link_index, false);
    if (std::abs(constraints_[link_index].lower_bound - cur_capacity) < link_capacity_selection_threshold_ &&
        cur_value <= middle_bound_value) {
      return false;
    }
    if (std::abs(constraints_[link_index].upper_bound - cur_capacity) < link_capacity_selection_threshold_ &&
        cur_value >= middle_bound_value) {
      return false;
    }
    return true;
  }

  bool CheckBudgetSatisfaction(T amount, int link_index) {
    network_.mutable_links()[link_index].capacity += amount;
    bool result = budget_upper_bound_ - Budget() > -budget_threshold_;
    network_.mutable_links()[link_index].capacity -= amount;
    return result;
  }

  T NonUpperBoundRespectedAmountConstranints(T amount, int max_optimality_index) {
    auto max_optimality_capacity = network_.mutable_links()[max_optimality_index].capacity;
    if (amount < 0) {
      return std::max(amount, constraints_[max_optimality_index].lower_bound - max_optimality_capacity);
    }
    else {
      amount = std::min(amount, constraints_[max_optimality_index].upper_bound - max_optimality_capacity);
      amount = std::min(amount, (budget_upper_bound_ - Budget()) / budget_function_multiplier_);
    }
    return amount;
  }

  bool NonUpperBoundTransferResult(T amount, int max_optimality_index, T middle_bound_value) {
    network_.mutable_links()[max_optimality_index].capacity += amount;
    bool result = false;
    if (amount < 0) {
      result = OptimalityFunction(max_optimality_index) < middle_bound_value;
    }
    else {
      result = OptimalityFunction(max_optimality_index) > middle_bound_value;
    }
    network_.mutable_links()[max_optimality_index].capacity -= amount;
    return result;
  }

  T NonUpperBoundOptimalityConditionTransferAmount(int max_optimality_index, T middle_bound_value) {
    T amount = 1;
    T next_amount = 0;
    auto max_optimality_capacity = network_.mutable_links()[max_optimality_index].capacity;
    amount = NonUpperBoundRespectedAmountConstranints(amount, max_optimality_index);
    auto cur_value = OptimalityFunction(max_optimality_index, false);
    if (cur_value < middle_bound_value) {
      next_amount = constraints_[max_optimality_index].lower_bound - max_optimality_capacity;
      if (NonUpperBoundTransferResult(next_amount, max_optimality_index, middle_bound_value)) {
        return next_amount;
      }
      amount = -1;
      amount = NonUpperBoundRespectedAmountConstranints(amount, max_optimality_index);
      next_amount = NonUpperBoundRespectedAmountConstranints(amount * -2, max_optimality_index);
      while (NonUpperBoundTransferResult(next_amount, max_optimality_index, middle_bound_value)) {
        amount = next_amount;
        next_amount = NonUpperBoundRespectedAmountConstranints(amount * -2, max_optimality_index);
        if (std::abs(next_amount - amount) < link_capacity_selection_threshold_ / 10) {
          break;
        }
      }
    }
    else {
      next_amount = NonUpperBoundRespectedAmountConstranints(constraints_[max_optimality_index].upper_bound - max_optimality_capacity, max_optimality_index);
      while (!NonUpperBoundTransferResult(next_amount, max_optimality_index, middle_bound_value)) {
        next_amount /= 2;
      }
      amount = next_amount;
    }
    return amount;
  }

  void NonUpperBoundBudgetOptimalityIteration() {
    T middle_bound_value = 1.0;
    int max_optimality_index = -1;
    
    approach_->Reset();
    approach_->ComputeTrafficFlows();

    for (int link_index = 0; link_index < network_.number_of_links(); link_index++) {
      if (ProcessingCheck(link_index, middle_bound_value)) {
        max_optimality_index = link_index;
        break;
      }
    }
    //std::cout << "max_optimality_index = " << max_optimality_index << '\n';
    if (max_optimality_index == -1) {
      return;
    }

    T max_delta = std::abs(middle_bound_value - OptimalityFunction(max_optimality_index, false));
    for (int link_index = 0; link_index < network_.number_of_links(); link_index++) {
      if (ProcessingCheck(link_index, middle_bound_value)) {
        auto cur_delta = std::abs(middle_bound_value - OptimalityFunction(link_index, false));
        if (cur_delta > max_delta) {
          max_optimality_index = link_index;
          max_delta = cur_delta;
        }
      }
    }
    if (!optimality_condition_progress_active_) {
      std::cout << "max_optimality_index = " << max_optimality_index << '\n';
      std::cout << "max_delta " << max_delta << '\n';
    }
    auto amount = NonUpperBoundOptimalityConditionTransferAmount(max_optimality_index, middle_bound_value);
    if (!optimality_condition_progress_active_) {
      std::cout << "amount " << amount << '\n';
    }
    network_.mutable_links()[max_optimality_index].capacity += amount;
  }

  BilevelCND(const BilevelCND&) = delete;
  BilevelCND& operator=(const BilevelCND&) = delete;
};

}  // namespace TrafficAssignment

#endif  // BILEVEL_CND_H

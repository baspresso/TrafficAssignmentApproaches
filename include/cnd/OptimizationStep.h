#ifndef OPTIMIZATION_STEP_H
#define OPTIMIZATION_STEP_H

#include <string>

namespace TrafficAssignment {

/**
 * @brief Configuration for a single optimization step in the CNDP pipeline.
 *
 * Each [step.N] INI section maps to one OptimizationStepConfig. The pipeline
 * executes steps sequentially, passing the best capacities forward as warm start.
 */
struct OptimizationStepConfig {
  std::string type;           ///< Step type: "nlopt", "optimality_condition", or "optimlib".
  std::string name;           ///< Optional display name; auto-generated from type+algorithm if empty.
  int max_iterations = 100;   ///< Maximum iterations/evaluations for this step.
  double tolerance = 1e-4;    ///< Convergence tolerance (relative objective change).

  // NLopt-specific
  std::string algorithm;          ///< NLopt algorithm name (e.g., "LN_COBYLA", "GN_ISRES").
  std::string local_algorithm;    ///< Local optimizer for composite algorithms (AUGLAG, global+local).
  int local_max_iterations = 0;   ///< Max evaluations for local optimizer (0 = max_iterations/10).
  double local_tolerance = 0.0;   ///< Convergence tolerance for local optimizer (0 = use main tolerance).

  // Population-based (OptimLib) specific
  int population_size = 0;  ///< Population size for DE/PSO (0 = auto: max(200, 2*n_vars)).

  // Gradient descent specific
  double step_size = 1.0;    ///< Initial step size for gradient descent (Armijo starting alpha).
  double fd_epsilon = 1e-4;  ///< Finite difference perturbation size for gradient estimation.
  std::string gradient_method;  ///< Gradient estimator: "finite_difference" (default), "spsa", "sensitivity".
  std::string stochastic_optimizer;  ///< Optimizer for stochastic gradient: "sgd" (default), "momentum", "adam".
};

/**
 * @brief Result returned by a single optimization step execution.
 */
struct StepResult {
  double objective = 0.0;         ///< Best objective value F(y) = TSTT(x(y)) + theta * BudgetCost(y).
  double total_travel_time = 0.0; ///< Total system travel time at the best solution.
  double budget = 0.0;            ///< Investment cost theta * sum(y_a - lb_a) at the best solution.
  bool success = true;            ///< Whether the step completed without error.
  int evaluations = 0;            ///< Number of objective function evaluations performed.
};

template <typename T>
class CndOptimizationContext;

/**
 * @brief Abstract interface for a single optimization step in the CNDP pipeline.
 *
 * Each step implements a different optimization strategy (derivative-free, sensitivity-based,
 * or population-based metaheuristic) that operates on the shared CndOptimizationContext
 * to improve link capacity decisions.
 *
 * @tparam T Numeric type for flow/capacity computations.
 */
template <typename T>
class OptimizationStep {
public:
  virtual ~OptimizationStep() = default;

  /// @brief Executes the optimization step, modifying network capacities via ctx.
  /// @return StepResult with the best objective, evaluations used, and success status.
  virtual StepResult Execute(CndOptimizationContext<T>& ctx) = 0;

  /// @brief Returns a human-readable name for this step (used in logging and traces).
  virtual std::string GetName() const = 0;
};

}  // namespace TrafficAssignment

#endif  // OPTIMIZATION_STEP_H

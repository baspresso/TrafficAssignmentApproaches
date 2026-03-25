#ifndef GRADIENT_ESTIMATOR_H
#define GRADIENT_ESTIMATOR_H

#include <functional>
#include <string>
#include <vector>

#include "../CndOptimizationContext.h"

namespace TrafficAssignment {

/**
 * @brief Abstract interface for gradient estimation strategies.
 *
 * Used by GradientDescentStep to compute gradients of the bilevel CNDP
 * objective. Implementations differ in accuracy vs. cost:
 * - FiniteDifference: n_vars TAP solves per gradient (exact but slow)
 * - SPSA: 2 TAP solves per gradient (noisy but dimension-independent)
 * - Sensitivity: 1 TAP solve + matrix ops (analytical, requires routes)
 *
 * @tparam T Numeric type for flow/capacity computations.
 */
template <typename T>
class GradientEstimator {
public:
  virtual ~GradientEstimator() = default;

  /// @brief Human-readable name for logging.
  virtual std::string GetName() const = 0;

  /**
   * @brief Computes the gradient of the objective at point x.
   *
   * @param ctx Shared optimization context (network, TAP solver, constraints).
   * @param x Current capacity vector.
   * @param f_x Objective value at x (precomputed).
   * @param lower_bounds Per-variable lower bounds.
   * @param upper_bounds Per-variable upper bounds.
   * @param n_vars Number of design variables.
   * @param[out] grad Output gradient vector (must be pre-allocated to n_vars).
   * @param eval_fn Black-box evaluation callback wrapping EvaluateObjectiveQuiet.
   *               Takes a capacity vector, returns objective value.
   */
  virtual void ComputeGradient(
    CndOptimizationContext<T>& ctx,
    const std::vector<double>& x, double f_x,
    const std::vector<double>& lower_bounds,
    const std::vector<double>& upper_bounds,
    int n_vars, std::vector<long double>& grad,
    std::function<double(const std::vector<double>&)> eval_fn
  ) = 0;

  /// @brief Number of TAP evaluations used in the last ComputeGradient call.
  virtual int LastEvalCount() const = 0;

  /// @brief Expected number of TAP evaluations per gradient for progress bar sizing.
  virtual int ExpectedEvalsPerGradient(int n_vars) const = 0;

  /// @brief Whether this estimator produces deterministic (low-noise) gradients.
  /// Deterministic estimators (FD, Sensitivity) work well with Armijo line search.
  /// Stochastic estimators (SPSA) need a decaying step schedule instead.
  virtual bool IsDeterministic() const { return true; }
};

}  // namespace TrafficAssignment

#endif  // GRADIENT_ESTIMATOR_H

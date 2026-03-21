#ifndef SENSITIVITY_ESTIMATOR_H
#define SENSITIVITY_ESTIMATOR_H

#include <cmath>
#include <vector>
#include <string>
#include <functional>
#include <iostream>
#include <memory>

#include "GradientEstimator.h"
#include "FiniteDifferenceEstimator.h"
#include "../../sensitivity/OptimalitySensitivity.h"

namespace TrafficAssignment {

/**
 * @brief Analytical gradient estimator using sensitivity analysis.
 *
 * Uses the Jacobian-based optimality function (Chiou 2005) to compute the
 * analytical gradient of the bilevel objective:
 *   grad[a] = theta * (phi_a + 1.0)
 * where phi_a is the optimality function for link a.
 *
 * Requires route information to be available (either from RouteBased approach
 * or post-solve extraction from TAPAS bushes). Falls back to finite differences
 * if no routes are available.
 *
 * Cost: 1 TAP solve + O(n_od * R^2) matrix operations.
 *
 * @tparam T Numeric type for flow/capacity computations.
 */
template <typename T>
class SensitivityEstimator : public GradientEstimator<T> {
public:
  explicit SensitivityEstimator(double fd_epsilon)
    : fd_fallback_(std::make_unique<FiniteDifferenceEstimator<T>>(fd_epsilon)) {}

  std::string GetName() const override { return "sensitivity"; }

  void ComputeGradient(
    CndOptimizationContext<T>& ctx,
    const std::vector<double>& x, double f_x,
    const std::vector<double>& lower_bounds,
    const std::vector<double>& upper_bounds,
    int n_vars, std::vector<long double>& grad,
    std::function<double(const std::vector<double>&)> eval_fn
  ) override {
    last_eval_count_ = 0;

    // Ensure routes are available (TAPAS: extract from bush; RouteBased: already there)
    ctx.approach->PopulateRoutes();

    // Check if any OD pair has routes
    if (!sensitivity_.HasRoutes(ctx)) {
      if (!fallback_warned_) {
        std::cout << "\n  [WARN] SensitivityEstimator: no routes available, falling back to FD\n";
        fallback_warned_ = true;
      }
      fd_fallback_->ComputeGradient(ctx, x, f_x, lower_bounds, upper_bounds,
                                     n_vars, grad, eval_fn);
      last_eval_count_ = fd_fallback_->LastEvalCount();
      return;
    }

    // Build OD cache (Jacobian inverses, denominators)
    auto od_cache = sensitivity_.BuildOdCache(ctx);

    const long double theta = static_cast<long double>(ctx.budget_function_multiplier);

    // Compute analytical gradient: grad[a] = theta * (1 - phi_a)
    // phi_a = -dTSTT/dy_a / theta (positive when increasing capacity reduces TSTT).
    // OptimalityFunctionFromCache applies a negation internally (see local_value = -num/den).
    // The full objective is F = TSTT + theta * sum(y_a - lb_a) + penalty
    // dF/dy_a = dTSTT/dy_a + theta = -theta*phi_a + theta = theta * (1 - phi_a)
    for (int a = 0; a < n_vars; ++a) {
      T phi_a = sensitivity_.OptimalityFunctionFromCache(ctx, a, od_cache);
      if (!ctx.IsFiniteScalar(phi_a)) {
        grad[a] = 0.0L;
        continue;
      }
      grad[a] = theta * (1.0L - static_cast<long double>(phi_a));
    }

    // Add soft budget penalty derivative if budget is violated
    double budget_function = ctx.BudgetFunction(n_vars, x.data());
    long double budget_violation = static_cast<long double>(budget_function) - static_cast<long double>(ctx.budget_upper_bound);
    if (budget_violation > 0.0L) {
      // d/dy_a [penalty_factor * (B - B_max)^2]
      //   = penalty_factor * 2 * (B - B_max) * dB/dy_a
      //   = penalty_factor * 2 * (B - B_max) * theta
      long double penalty_grad = static_cast<long double>(CndOptimizationContext<T>::kDefaultBudgetViolationPenaltyFactor)
                                 * 2.0L * budget_violation * theta;
      for (int a = 0; a < n_vars; ++a) {
        grad[a] += penalty_grad;
      }
    }
  }

  int LastEvalCount() const override { return last_eval_count_; }
  int ExpectedEvalsPerGradient(int /*n_vars*/) const override { return 1; }

private:
  OptimalitySensitivity<T> sensitivity_;
  std::unique_ptr<FiniteDifferenceEstimator<T>> fd_fallback_;
  bool fallback_warned_ = false;
  int last_eval_count_ = 0;
};

}  // namespace TrafficAssignment

#endif  // SENSITIVITY_ESTIMATOR_H

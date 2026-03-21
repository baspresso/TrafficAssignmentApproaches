#ifndef FINITE_DIFFERENCE_ESTIMATOR_H
#define FINITE_DIFFERENCE_ESTIMATOR_H

#include <cmath>
#include <vector>
#include <string>
#include <functional>

#include "GradientEstimator.h"

namespace TrafficAssignment {

/**
 * @brief Forward finite-difference gradient estimator.
 *
 * For each variable i, perturbs x[i] by h_i = max(epsilon, epsilon * |x[i]|).
 * If x[i] is at the upper bound, uses backward difference instead.
 * Cost: n_vars TAP solves per gradient.
 *
 * @tparam T Numeric type for flow/capacity computations.
 */
template <typename T>
class FiniteDifferenceEstimator : public GradientEstimator<T> {
public:
  explicit FiniteDifferenceEstimator(double fd_epsilon)
    : fd_epsilon_(fd_epsilon) {}

  std::string GetName() const override { return "finite_difference"; }

  void ComputeGradient(
    CndOptimizationContext<T>& ctx,
    const std::vector<double>& x, double f_x,
    const std::vector<double>& lower_bounds,
    const std::vector<double>& upper_bounds,
    int n_vars, std::vector<long double>& grad,
    std::function<double(const std::vector<double>&)> eval_fn
  ) override {
    last_eval_count_ = 0;

    for (int i = 0; i < n_vars; ++i) {
      long double h = std::max(static_cast<long double>(fd_epsilon_),
                               static_cast<long double>(fd_epsilon_) * std::abs(static_cast<long double>(x[i])));

      bool backward = (x[i] + static_cast<double>(h) > upper_bounds[i]);

      std::vector<double> x_pert = x;
      if (backward) {
        x_pert[i] = x[i] - static_cast<double>(h);
        x_pert[i] = std::max(lower_bounds[i], x_pert[i]);
        x_pert[i] = std::max(kMinCapacity, x_pert[i]);
        h = static_cast<long double>(x[i] - x_pert[i]);
        if (h < 1e-15L) {
          grad[i] = 0.0L;
          continue;
        }
        long double f_pert = static_cast<long double>(eval_fn(x_pert));
        ++last_eval_count_;
        grad[i] = (static_cast<long double>(f_x) - f_pert) / h;
      } else {
        x_pert[i] = x[i] + static_cast<double>(h);
        x_pert[i] = std::min(upper_bounds[i], x_pert[i]);
        h = static_cast<long double>(x_pert[i] - x[i]);
        if (h < 1e-15L) {
          grad[i] = 0.0L;
          continue;
        }
        long double f_pert = static_cast<long double>(eval_fn(x_pert));
        ++last_eval_count_;
        grad[i] = (f_pert - static_cast<long double>(f_x)) / h;
      }
    }
  }

  int LastEvalCount() const override { return last_eval_count_; }
  int ExpectedEvalsPerGradient(int n_vars) const override { return n_vars; }

private:
  static constexpr double kMinCapacity = 1e-6;

  double fd_epsilon_;
  int last_eval_count_ = 0;
};

}  // namespace TrafficAssignment

#endif  // FINITE_DIFFERENCE_ESTIMATOR_H

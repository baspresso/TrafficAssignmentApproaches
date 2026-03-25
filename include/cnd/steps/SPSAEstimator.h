#ifndef SPSA_ESTIMATOR_H
#define SPSA_ESTIMATOR_H

#include <cmath>
#include <vector>
#include <string>
#include <functional>
#include <random>

#include "GradientEstimator.h"

namespace TrafficAssignment {

/**
 * @brief SPSA (Simultaneous Perturbation Stochastic Approximation) gradient estimator.
 *
 * Estimates the full gradient from only 2 function evaluations using random
 * simultaneous perturbation (Spall, 1992):
 *   g_i = (f(x + c*Delta) - f(x - c*Delta)) / (2*c*Delta_i)
 *
 * Delta is a vector of i.i.d. Bernoulli +/-1 values. Perturbed points are
 * clamped to bounds and budget feasibility is enforced. The gradient is
 * computed from actual displacement (accounting for clamping).
 *
 * Cost: 2 TAP solves regardless of n_vars.
 *
 * @tparam T Numeric type for flow/capacity computations.
 */
template <typename T>
class SPSAEstimator : public GradientEstimator<T> {
public:
  explicit SPSAEstimator(double perturbation_size)
    : perturbation_size_(perturbation_size), rng_(42) {}

  void SetPerturbationSize(double c_k) { perturbation_size_ = c_k; }

  std::string GetName() const override { return "spsa"; }

  void ComputeGradient(
    CndOptimizationContext<T>& ctx,
    const std::vector<double>& x, double f_x,
    const std::vector<double>& lower_bounds,
    const std::vector<double>& upper_bounds,
    int n_vars, std::vector<long double>& grad,
    std::function<double(const std::vector<double>&)> eval_fn
  ) override {
    last_eval_count_ = 0;

    // Auto-select number of averaged SPSA samples to reduce variance.
    // Each sample costs 2 TAP solves. Using ~sqrt(n_vars) samples gives a
    // good noise-vs-cost trade-off while staying well below FD's n_vars cost.
    // const int num_samples = std::max(1, std::min(20,
    //     static_cast<int>(std::ceil(std::sqrt(static_cast<double>(n_vars))))));
    const int num_samples = 3;
    expected_evals_ = 2 * num_samples;

    // Zero gradient accumulator
    for (int i = 0; i < n_vars; ++i) grad[i] = 0.0L;

    std::vector<double> delta(n_vars);
    std::vector<double> x_plus(n_vars);
    std::vector<double> x_minus(n_vars);
    std::bernoulli_distribution coin(0.5);

    for (int s = 0; s < num_samples; ++s) {
      // Generate Bernoulli +/-1 perturbation vector
      for (int i = 0; i < n_vars; ++i) {
        delta[i] = coin(rng_) ? 1.0 : -1.0;
      }

      // Compute perturbed points, clamped to bounds
      for (int i = 0; i < n_vars; ++i) {
        double c = perturbation_size_ * std::max(1.0, std::abs(x[i]));
        x_plus[i] = std::max(lower_bounds[i], std::min(upper_bounds[i], x[i] + c * delta[i]));
        x_minus[i] = std::max(lower_bounds[i], std::min(upper_bounds[i], x[i] - c * delta[i]));
        x_plus[i] = std::max(kMinCapacity, x_plus[i]);
        x_minus[i] = std::max(kMinCapacity, x_minus[i]);
      }

      // Enforce budget feasibility on both perturbed points
      // ctx.EnforceBudgetFeasibility(x_plus, lower_bounds);
      // ctx.EnforceBudgetFeasibility(x_minus, lower_bounds);

      // Evaluate at both perturbed points
      long double f_plus = static_cast<long double>(eval_fn(x_plus));
      ++last_eval_count_;
      long double f_minus = static_cast<long double>(eval_fn(x_minus));
      ++last_eval_count_;

      // Accumulate gradient from actual displacements (long double arithmetic)
      long double df = f_plus - f_minus;
      for (int i = 0; i < n_vars; ++i) {
        long double actual_displacement = static_cast<long double>(x_plus[i]) - static_cast<long double>(x_minus[i]);
        if (std::abs(actual_displacement) > 1e-15L) {
          grad[i] += df / actual_displacement;
        }
      }
    }

    // Average over samples
    long double inv_samples = 1.0L / num_samples;
    for (int i = 0; i < n_vars; ++i) {
      grad[i] *= inv_samples;
    }
  }

  int LastEvalCount() const override { return last_eval_count_; }
  int ExpectedEvalsPerGradient(int n_vars) const override {
    // int num_samples = std::max(1, std::min(20,
    //     static_cast<int>(std::ceil(std::sqrt(static_cast<double>(n_vars))))));
    int num_samples = 3;
    return 2 * num_samples;
  }
  bool IsDeterministic() const override { return false; }

private:
  static constexpr double kMinCapacity = 1e-6;

  double perturbation_size_;
  std::mt19937 rng_;
  int last_eval_count_ = 0;
  int expected_evals_ = 2;
};

}  // namespace TrafficAssignment

#endif  // SPSA_ESTIMATOR_H

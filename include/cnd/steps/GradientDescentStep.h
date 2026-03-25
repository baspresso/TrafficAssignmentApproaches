#ifndef GRADIENT_DESCENT_STEP_H
#define GRADIENT_DESCENT_STEP_H

#include <vector>
#include <string>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <iostream>
#include <memory>
#include <functional>

#include "../OptimizationStep.h"
#include "../CndOptimizationContext.h"
#include "GradientEstimator.h"
#include "FiniteDifferenceEstimator.h"
#include "SPSAEstimator.h"
#include "SensitivityEstimator.h"

namespace TrafficAssignment {

/**
 * @brief Projected gradient descent optimization step with finite-difference gradients.
 *
 * Computes gradients via forward finite differences (one TAP solve per variable per iteration),
 * then performs projected gradient descent with Armijo backtracking line search.
 * Box constraints [lb, ub] are enforced by projection, and budget feasibility is maintained
 * via the same proportional projection used by other steps.
 *
 * Best suited for small-medium networks or short refinement runs after another optimizer.
 * Each iteration costs ~(n_vars + line_search_evals) TAP solves.
 *
 * @tparam T Numeric type for flow/capacity computations.
 */
template <typename T>
class GradientDescentStep : public OptimizationStep<T> {
public:
  explicit GradientDescentStep(const OptimizationStepConfig& config)
    : config_(config) {
    // Factory: create gradient estimator based on config
    std::string method = config.gradient_method;
    std::transform(method.begin(), method.end(), method.begin(),
      [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (method == "spsa") {
      estimator_ = std::make_unique<SPSAEstimator<T>>(config.fd_epsilon);
    } else if (method == "sensitivity" || method == "analytical") {
      estimator_ = std::make_unique<SensitivityEstimator<T>>(config.fd_epsilon);
    } else {
      // Default: finite_difference (also handles "fd", "finite_difference", "")
      estimator_ = std::make_unique<FiniteDifferenceEstimator<T>>(config.fd_epsilon);
    }
  }

  std::string GetName() const override {
    if (!config_.name.empty()) return config_.name;
    std::string name = "gradient_descent(" + estimator_->GetName();
    if (!config_.stochastic_optimizer.empty() && config_.stochastic_optimizer != "sgd") {
      name += "+" + config_.stochastic_optimizer;
    }
    name += ")";
    return name;
  }

  StepResult Execute(CndOptimizationContext<T>& ctx) override {
    StepResult result;
    if (config_.max_iterations == 0) {
      result.success = true;
      return result;
    }

    const int n_vars = static_cast<int>(ctx.constraints.size());

    // Read current capacities
    std::vector<double> x(n_vars);
    std::vector<double> lower_bounds(n_vars);
    std::vector<double> upper_bounds(n_vars);
    for (int i = 0; i < n_vars; ++i) {
      x[i] = static_cast<double>(ctx.network.mutable_links()[i].capacity);
      lower_bounds[i] = ctx.constraints[i].lower_bound;
      upper_bounds[i] = ctx.constraints[i].upper_bound;
    }

    // Evaluate initial objective
    double f_x = EvaluateObjective(ctx, x, lower_bounds, n_vars);

    // Track best solution
    std::vector<double> best_x = x;
    double best_f = f_x;

    // Progress bar: each iteration costs gradient_evals + a few line search evals
    const int evals_per_grad = estimator_->ExpectedEvalsPerGradient(n_vars);
    const bool stochastic = !estimator_->IsDeterministic();
    const int ls_evals = stochastic ? 1 : 3;  // stochastic: no line search, just 1 step eval
    const int total_evals_estimate = config_.max_iterations * (evals_per_grad + ls_evals);
    progress_.Start(total_evals_estimate);

    std::vector<long double> grad(n_vars);

    // Create evaluation callback for gradient estimators
    auto eval_fn = [this, &ctx, &lower_bounds, n_vars](const std::vector<double>& x_eval) -> double {
      return EvaluateObjectiveQuiet(ctx, x_eval, lower_bounds, n_vars);
    };

    if (stochastic) {
      ExecuteStochastic(ctx, x, f_x, lower_bounds, upper_bounds, n_vars,
                         best_x, best_f, grad, eval_fn);
    } else {
      ExecuteDeterministic(ctx, x, f_x, lower_bounds, upper_bounds, n_vars,
                            best_x, best_f, grad, eval_fn);
    }

    progress_.Finish();

    // Write best solution back to network
    for (int i = 0; i < n_vars; ++i) {
      ctx.network.mutable_links()[i].capacity = static_cast<T>(best_x[i]);
    }

    result.objective = best_f;
    result.total_travel_time = ctx.network.TotalTravelTime();
    result.budget = static_cast<double>(ctx.Budget());
    result.evaluations = eval_count_;
    result.success = true;

    return result;
  }

private:
  static constexpr double kArmijoC1 = 1e-4;        ///< Sufficient decrease parameter.
  static constexpr double kBacktrackFactor = 0.5;   ///< Step size reduction factor.
  static constexpr int kMaxBacktracks = 25;         ///< Max line search attempts.
  static constexpr double kMinGradientNorm = 1e-12; ///< Convergence threshold on gradient norm.
  static constexpr double kMinCapacity = 1e-6;      ///< Minimum capacity clamp.
  static constexpr double kMinStepSize = 1e-8;      ///< Minimum step size before giving up.
  static constexpr int kMaxLineSearchFailures = 3;  ///< Max consecutive line search failures.
  // Spall (1998) SPSA constants for step/perturbation decay schedules.
  static constexpr double kSpsaAlpha = 0.602;          ///< Step size decay exponent: a_k = a0 / (k+1+A)^alpha.
  static constexpr double kSpsaGamma = 0.101;          ///< Perturbation decay exponent: c_k = c0 / (k+1)^gamma.
  static constexpr double kSpsaStabilityFraction = 0.1; ///< Stability constant A = fraction * max_iterations.
  static constexpr double kGradClipMultiplier = 10.0;  ///< Clip gradient components exceeding this * mean(|g|).
  // Adam optimizer constants (Kingma & Ba 2015).
  static constexpr double kAdamBeta1 = 0.9;     ///< First moment decay rate.
  static constexpr double kAdamBeta2 = 0.999;    ///< Second moment decay rate.
  static constexpr double kAdamEpsilon = 1e-8;   ///< Numerical stability constant for Adam denominator.
  // Momentum constant.
  static constexpr double kMomentumBeta = 0.9;   ///< Momentum decay rate for gradient EMA.

  OptimizationStepConfig config_;
  std::unique_ptr<GradientEstimator<T>> estimator_;
  typename CndOptimizationContext<T>::ProgressState progress_;
  int eval_count_ = 0;
  int consecutive_ls_failures_ = 0;

  /**
   * @brief Deterministic gradient descent with Armijo line search.
   * Used for FD and Sensitivity estimators that produce reliable gradients.
   */
  void ExecuteDeterministic(CndOptimizationContext<T>& ctx,
                             std::vector<double>& x, double f_x,
                             const std::vector<double>& lower_bounds,
                             const std::vector<double>& upper_bounds,
                             int n_vars,
                             std::vector<double>& best_x, double& best_f,
                             std::vector<long double>& grad,
                             const std::function<double(const std::vector<double>&)>& eval_fn) {
    double current_step_size = config_.step_size;

    for (int iter = 0; iter < config_.max_iterations; ++iter) {
      estimator_->ComputeGradient(ctx, x, f_x, lower_bounds, upper_bounds,
                                   n_vars, grad, eval_fn);

      // Re-evaluate f_x at the base point after gradient computation.
      // The gradient evaluations ran TAP solves from perturbed capacities,
      // so the TAP warm-start state has drifted.
      f_x = EvaluateObjectiveQuiet(ctx, x, lower_bounds, n_vars);
      if (f_x < best_f) { best_f = f_x; best_x = x; }

      long double grad_norm_sq = 0.0L;
      for (int i = 0; i < n_vars; ++i) grad_norm_sq += grad[i] * grad[i];
      if (grad_norm_sq < static_cast<long double>(kMinGradientNorm) * static_cast<long double>(kMinGradientNorm)) {
        if (ctx.verbose) {
          std::cout << "\n  GD converged: gradient norm " << std::sqrt(grad_norm_sq)
                    << " < " << kMinGradientNorm << std::endl;
        }
        break;
      }

      long double grad_norm = std::sqrt(grad_norm_sq);
      for (int i = 0; i < n_vars; ++i) grad[i] /= grad_norm;

      double alpha = ArmijoLineSearch(ctx, x, f_x, grad,
                                       lower_bounds, upper_bounds, n_vars,
                                       current_step_size);

      if (alpha <= 0.0) {
        ++consecutive_ls_failures_;
        current_step_size *= 0.1;
        if (consecutive_ls_failures_ >= kMaxLineSearchFailures ||
            current_step_size < kMinStepSize) {
          if (ctx.verbose) {
            std::cout << "\n  GD: line search failed " << consecutive_ls_failures_
                      << " times, stopping at iteration " << iter << std::endl;
          }
          break;
        }
        if (ctx.verbose) {
          std::cout << "\n  GD: line search failed at iteration " << iter
                    << ", shrinking step to " << current_step_size << std::endl;
        }
        continue;
      }
      consecutive_ls_failures_ = 0;

      if (alpha >= current_step_size * 0.5) {
        current_step_size = std::min(current_step_size * 1.5, config_.step_size * 10.0);
      } else {
        current_step_size = std::max(alpha, config_.step_size * 0.01);
      }

      for (int i = 0; i < n_vars; ++i) {
        double x_new_i = x[i] - alpha * static_cast<double>(grad[i]);
        x_new_i = std::max(lower_bounds[i], std::min(upper_bounds[i], x_new_i));
        x_new_i = std::max(kMinCapacity, x_new_i);
        x[i] = x_new_i;
      }
      ctx.EnforceBudgetFeasibility(x, lower_bounds);

      double f_new = EvaluateObjective(ctx, x, lower_bounds, n_vars);
      if (f_new < best_f) { best_f = f_new; best_x = x; }

      double rel_change = std::abs(f_new - f_x) / std::max(1.0, std::abs(f_x));
      if (rel_change < config_.tolerance && f_new <= f_x) {
        if (ctx.verbose) {
          std::cout << "\n  GD converged: relative change " << rel_change
                    << " < " << config_.tolerance << std::endl;
        }
        break;
      }
      f_x = f_new;
    }
  }

  /**
   * @brief Stochastic gradient descent with Spall (1998) SPSA schedule.
   *
   * Uses standard SPSA decay schedules with convergence guarantees:
   *   a_k = a0 / (k+1+A)^alpha    (step size, alpha=0.602)
   *   c_k = c0 / (k+1)^gamma      (perturbation size, gamma=0.101)
   *
   * Supports three optimizer modes (selected via config_.stochastic_optimizer):
   *   - "sgd" (default): Raw gradient step with Spall schedule.
   *   - "momentum": Gradient EMA (heavy ball) for noise reduction.
   *   - "adam": Per-variable adaptive learning rates (Kingma & Ba 2015).
   *
   * Gradient clipping tames noisy outlier components.
   * Polyak-Ruppert iterate averaging after burn-in provides a smoothed solution.
   */
  void ExecuteStochastic(CndOptimizationContext<T>& ctx,
                          std::vector<double>& x, double f_x,
                          const std::vector<double>& lower_bounds,
                          const std::vector<double>& upper_bounds,
                          int n_vars,
                          std::vector<double>& best_x, double& best_f,
                          std::vector<long double>& grad,
                          const std::function<double(const std::vector<double>&)>& eval_fn) {
    static constexpr int kStagnationWindow = 150;
    int iters_since_improvement = 0;

    // Determine optimizer mode
    enum class OptimizerMode { kSgd, kMomentum, kAdam };
    OptimizerMode mode = OptimizerMode::kSgd;
    {
      std::string opt = config_.stochastic_optimizer;
      std::transform(opt.begin(), opt.end(), opt.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
      if (opt == "momentum") mode = OptimizerMode::kMomentum;
      else if (opt == "adam") mode = OptimizerMode::kAdam;
    }

    // Spall schedule parameters
    const long double a0 = static_cast<long double>(config_.step_size);
    const long double c0 = static_cast<long double>(config_.fd_epsilon);
    const long double A = static_cast<long double>(kSpsaStabilityFraction) * config_.max_iterations;

    // Polyak-Ruppert iterate averaging after burn-in
    const int burn_in = static_cast<int>(0.2 * config_.max_iterations);
    std::vector<double> x_avg(n_vars, 0.0);
    int avg_count = 0;

    // Momentum state: gradient EMA
    std::vector<long double> g_ema(n_vars, 0.0L);

    // Adam state: first and second moment estimates
    std::vector<long double> adam_m(n_vars, 0.0L);
    std::vector<long double> adam_v(n_vars, 0.0L);

    // Get SPSA estimator pointer for perturbation decay
    auto* spsa_est = dynamic_cast<SPSAEstimator<T>*>(estimator_.get());

    if (ctx.verbose && mode != OptimizerMode::kSgd) {
      std::cout << "  GD(SPSA): using "
                << (mode == OptimizerMode::kMomentum ? "momentum" : "adam")
                << " optimizer" << std::endl;
    }

    for (int iter = 0; iter < config_.max_iterations; ++iter) {
      // Decaying perturbation: c_k = c0 / (k+1)^gamma
      long double c_k = c0 / std::pow(static_cast<long double>(iter + 1), static_cast<long double>(kSpsaGamma));
      if (spsa_est) {
        spsa_est->SetPerturbationSize(static_cast<double>(c_k));
      }

      estimator_->ComputeGradient(ctx, x, f_x, lower_bounds, upper_bounds,
                                   n_vars, grad, eval_fn);

      // Skip near-zero gradients
      long double grad_norm_sq = 0.0L;
      for (int i = 0; i < n_vars; ++i) grad_norm_sq += grad[i] * grad[i];
      if (grad_norm_sq < static_cast<long double>(kMinGradientNorm) * static_cast<long double>(kMinGradientNorm)) continue;

      // Gradient clipping: clip components exceeding kGradClipMultiplier * mean(|g|)
      long double mean_abs_grad = 0.0L;
      for (int i = 0; i < n_vars; ++i) mean_abs_grad += std::abs(grad[i]);
      mean_abs_grad /= n_vars;
      long double clip_threshold = static_cast<long double>(kGradClipMultiplier) * mean_abs_grad;
      if (clip_threshold > 0.0L) {
        for (int i = 0; i < n_vars; ++i) {
          grad[i] = std::max(-clip_threshold, std::min(clip_threshold, grad[i]));
        }
      }

      // Decaying step size: a_k = a0 / (k+1+A)^alpha
      long double a_k = a0 / std::pow(static_cast<long double>(iter + 1) + A, static_cast<long double>(kSpsaAlpha));

      // Compute update direction based on optimizer mode
      if (mode == OptimizerMode::kAdam) {
        // Adam: per-variable adaptive learning rates
        const int t = iter + 1;  // 1-indexed for bias correction
        for (int i = 0; i < n_vars; ++i) {
          adam_m[i] = static_cast<long double>(kAdamBeta1) * adam_m[i]
                    + (1.0L - static_cast<long double>(kAdamBeta1)) * grad[i];
          adam_v[i] = static_cast<long double>(kAdamBeta2) * adam_v[i]
                    + (1.0L - static_cast<long double>(kAdamBeta2)) * grad[i] * grad[i];
        }
        // Bias-corrected estimates
        long double bc1 = 1.0L - std::pow(static_cast<long double>(kAdamBeta1), static_cast<long double>(t));
        long double bc2 = 1.0L - std::pow(static_cast<long double>(kAdamBeta2), static_cast<long double>(t));
        for (int i = 0; i < n_vars; ++i) {
          long double m_hat = adam_m[i] / bc1;
          long double v_hat = adam_v[i] / bc2;
          double x_new_i = x[i] - static_cast<double>(a_k * m_hat / (std::sqrt(v_hat) + static_cast<long double>(kAdamEpsilon)));
          x_new_i = std::max(lower_bounds[i], std::min(upper_bounds[i], x_new_i));
          x_new_i = std::max(kMinCapacity, x_new_i);
          x[i] = x_new_i;
        }
      } else if (mode == OptimizerMode::kMomentum) {
        // Momentum: gradient EMA (heavy ball)
        for (int i = 0; i < n_vars; ++i) {
          g_ema[i] = static_cast<long double>(kMomentumBeta) * g_ema[i]
                   + (1.0L - static_cast<long double>(kMomentumBeta)) * grad[i];
          double x_new_i = x[i] - static_cast<double>(a_k * g_ema[i]);
          x_new_i = std::max(lower_bounds[i], std::min(upper_bounds[i], x_new_i));
          x_new_i = std::max(kMinCapacity, x_new_i);
          x[i] = x_new_i;
        }
      } else {
        // SGD: raw gradient step (original behavior)
        for (int i = 0; i < n_vars; ++i) {
          double x_new_i = x[i] - static_cast<double>(a_k * grad[i]);
          x_new_i = std::max(lower_bounds[i], std::min(upper_bounds[i], x_new_i));
          x_new_i = std::max(kMinCapacity, x_new_i);
          x[i] = x_new_i;
        }
      }
      ctx.EnforceBudgetFeasibility(x, lower_bounds);

      f_x = EvaluateObjective(ctx, x, lower_bounds, n_vars);

      if (f_x < best_f) {
        best_f = f_x;
        best_x = x;
        iters_since_improvement = 0;
      } else {
        ++iters_since_improvement;
        if (iters_since_improvement >= kStagnationWindow) {
          if (ctx.verbose) {
            std::cout << "\n  GD(SPSA): no improvement in " << kStagnationWindow
                      << " iterations, stopping at iteration " << iter << std::endl;
          }
          break;
        }
      }

      // Accumulate iterate average after burn-in
      if (iter >= burn_in) {
        for (int i = 0; i < n_vars; ++i) {
          x_avg[i] += x[i];
        }
        ++avg_count;
      }
    }

    // Evaluate averaged iterate and use if better than tracked best
    if (avg_count > 0) {
      for (int i = 0; i < n_vars; ++i) {
        x_avg[i] /= avg_count;
        x_avg[i] = std::max(lower_bounds[i], std::min(upper_bounds[i], x_avg[i]));
        x_avg[i] = std::max(kMinCapacity, x_avg[i]);
      }
      ctx.EnforceBudgetFeasibility(x_avg, lower_bounds);
      double f_avg = EvaluateObjective(ctx, x_avg, lower_bounds, n_vars);
      if (f_avg < best_f) {
        best_f = f_avg;
        best_x = x_avg;
        if (ctx.verbose) {
          std::cout << "\n  GD(SPSA): iterate average improved objective: "
                    << f_avg << std::endl;
        }
      }
    }

    // Restore best solution
    x = best_x;
    f_x = best_f;
  }

  /**
   * @brief Armijo backtracking line search with projected directional derivative.
   *
   * Uses the correct projected gradient descent Armijo condition:
   *   f(x_trial) <= f(x) + c1 * grad^T * (x_trial - x)
   * where x_trial = project(x - alpha * grad) accounts for box + budget projection.
   *
   * Falls back to simple decrease (f_trial < f_x) if Armijo fails, to handle
   * noisy finite-difference gradients on the non-smooth bilevel objective.
   */
  double ArmijoLineSearch(CndOptimizationContext<T>& ctx,
                           const std::vector<double>& x,
                           double f_x,
                           const std::vector<long double>& grad,
                           const std::vector<double>& lower_bounds,
                           const std::vector<double>& upper_bounds,
                           int n_vars,
                           double initial_alpha) {
    double alpha = initial_alpha;

    std::vector<double> x_trial(n_vars);
    double best_decrease_alpha = 0.0;
    double best_decrease_f = f_x;

    for (int k = 0; k < kMaxBacktracks; ++k) {
      // Compute trial point with projection
      for (int i = 0; i < n_vars; ++i) {
        double trial = x[i] - alpha * static_cast<double>(grad[i]);
        trial = std::max(lower_bounds[i], std::min(upper_bounds[i], trial));
        trial = std::max(kMinCapacity, trial);
        x_trial[i] = trial;
      }
      ctx.EnforceBudgetFeasibility(x_trial, lower_bounds);

      // Compute projected directional derivative: grad^T * (x_trial - x)
      // This is negative when x_trial moves in a descent direction.
      long double directional_deriv = 0.0L;
      for (int i = 0; i < n_vars; ++i) {
        directional_deriv += grad[i] * static_cast<long double>(x_trial[i] - x[i]);
      }

      double f_trial = EvaluateObjectiveQuiet(ctx, x_trial, lower_bounds, n_vars);

      // Track best simple decrease for fallback
      if (f_trial < best_decrease_f) {
        best_decrease_f = f_trial;
        best_decrease_alpha = alpha;
      }

      // Armijo sufficient decrease condition using projected directional derivative.
      // If directional_deriv >= 0 the projected step is not descent — keep backtracking.
      if (directional_deriv < 0.0L &&
          static_cast<long double>(f_trial) <= static_cast<long double>(f_x) + static_cast<long double>(kArmijoC1) * directional_deriv) {
        return alpha;
      }

      alpha *= kBacktrackFactor;
    }

    // Fallback: accept any step that achieved simple decrease.
    // This handles noisy bilevel objectives where Armijo is too strict.
    if (best_decrease_alpha > 0.0) {
      return best_decrease_alpha;
    }

    return 0.0;  // No decrease found at all
  }

  /**
   * @brief Evaluates objective with logging and progress bar updates.
   *
   * Sets capacities, solves TAP, computes TSTT + budget cost + soft penalty.
   */
  double EvaluateObjective(CndOptimizationContext<T>& ctx,
                            const std::vector<double>& x,
                            const std::vector<double>& lower_bounds,
                            int n_vars) {
    double obj = EvaluateObjectiveCore(ctx, x, lower_bounds, n_vars, true);
    return obj;
  }

  /**
   * @brief Evaluates objective without progress bar updates (for gradient/line search).
   */
  double EvaluateObjectiveQuiet(CndOptimizationContext<T>& ctx,
                                 const std::vector<double>& x,
                                 const std::vector<double>& lower_bounds,
                                 int n_vars) {
    double obj = EvaluateObjectiveCore(ctx, x, lower_bounds, n_vars, false);
    return obj;
  }

  /**
   * @brief Core objective evaluation: set capacities, solve TAP, compute penalized objective.
   */
  double EvaluateObjectiveCore(CndOptimizationContext<T>& ctx,
                                const std::vector<double>& x,
                                const std::vector<double>& lower_bounds,
                                int n_vars,
                                bool show_progress) {
    // Set capacities
    for (int i = 0; i < n_vars; ++i) {
      ctx.network.mutable_links()[i].capacity = static_cast<T>(
        std::max(kMinCapacity, x[i])
      );
    }

    // Solve TAP
    double ta_compute_seconds = 0.0;
    const bool ta_valid = ctx.ComputeTrafficFlowsSafely(ta_compute_seconds, "gd_eval");

    double total_travel_time = ctx.network.TotalTravelTime();
    double budget_function = ctx.BudgetFunction(n_vars, x.data());
    const double budget_violation = std::max(0.0, budget_function - ctx.budget_upper_bound);
    double soft_budget_penalty = 0.0;
    if (budget_violation > 0.0) {
      soft_budget_penalty =
        CndOptimizationContext<T>::kDefaultBudgetViolationPenaltyFactor *
        budget_violation * budget_violation;
      if (!CndOptimizationContext<T>::IsFiniteNumber(soft_budget_penalty)) {
        soft_budget_penalty = CndOptimizationContext<T>::kInvalidObjectivePenalty;
      }
    }
    double objective_function = total_travel_time + budget_function + soft_budget_penalty;

    double log_total_travel_time = total_travel_time;
    double log_budget_function = budget_function;
    double log_objective_function = objective_function;
    double progress_total_travel_time = total_travel_time;
    double progress_budget_function = budget_function;
    double progress_objective_function = objective_function;

    bool objective_valid = ta_valid &&
      CndOptimizationContext<T>::IsFiniteNumber(total_travel_time) &&
      CndOptimizationContext<T>::IsFiniteNumber(budget_function) &&
      CndOptimizationContext<T>::IsFiniteNumber(objective_function);
    if (!objective_valid) {
      ctx.PrintInvalidStateWarning("gd_eval", "objective_invalid");
      log_total_travel_time = std::numeric_limits<double>::quiet_NaN();
      log_budget_function = std::numeric_limits<double>::quiet_NaN();
      log_objective_function = std::numeric_limits<double>::quiet_NaN();
      progress_total_travel_time = CndOptimizationContext<T>::kInvalidObjectivePenalty;
      progress_budget_function =
        CndOptimizationContext<T>::IsFiniteNumber(budget_function) ? budget_function : 0.0;
      progress_objective_function = CndOptimizationContext<T>::kInvalidObjectivePenalty;
      objective_function = CndOptimizationContext<T>::kInvalidObjectivePenalty;
    }

    ++ctx.counters.standard_trace_step;
    ++eval_count_;
    ctx.UpdateRuntimeAndBudgetMetrics(ta_compute_seconds, progress_budget_function);
    ctx.LogQualityTimePoint("gd_eval",
                            ctx.counters.standard_trace_step,
                            log_objective_function,
                            log_total_travel_time,
                            log_budget_function,
                            ta_compute_seconds);
    if (show_progress) {
      CndOptimizationContext<T>::PrintProgressBarWithMetrics(
        progress_, "eval/s",
        progress_objective_function,
        progress_total_travel_time,
        progress_budget_function
      );
    } else {
      CndOptimizationContext<T>::PrintProgressBar(progress_, " eval/s");
    }

    return objective_function;
  }
};

}  // namespace TrafficAssignment

#endif  // GRADIENT_DESCENT_STEP_H

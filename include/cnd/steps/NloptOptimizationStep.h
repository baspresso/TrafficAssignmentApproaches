#ifndef NLOPT_OPTIMIZATION_STEP_H
#define NLOPT_OPTIMIZATION_STEP_H

#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <nlopt.hpp>
#include "../OptimizationStep.h"
#include "../CndOptimizationContext.h"

namespace TrafficAssignment {

// Forward declaration for ParseNloptAlgorithm — defined at bottom of file
inline nlopt::algorithm ParseNloptAlgorithmForStep(const std::string& value);

template <typename T>
class NloptOptimizationStep : public OptimizationStep<T> {
public:
  explicit NloptOptimizationStep(const OptimizationStepConfig& config)
    : config_(config),
      ctx_(nullptr) {}

  std::string GetName() const override {
    if (!config_.name.empty()) return config_.name;
    return "nlopt(" + config_.algorithm + ")";
  }

  StepResult Execute(CndOptimizationContext<T>& ctx) override {
    StepResult result;
    if (config_.max_iterations == 0) {
      result.success = true;
      return result;
    }

    ctx_ = &ctx;

    const nlopt::algorithm algo = ParseNloptAlgorithmForStep(config_.algorithm);
    const bool hard_budget = CndOptimizationContext<T>::SupportsHardBudgetConstraint(algo);
    const double tolerance = config_.tolerance;

    const int n_vars = static_cast<int>(ctx.constraints.size());

    // Build initial x from current network capacities
    std::vector<double> x(n_vars);
    for (int i = 0; i < n_vars; ++i) {
      x[i] = static_cast<double>(ctx.network.mutable_links()[i].capacity);
    }

    nlopt::opt optimizer(algo, n_vars);
    optimizer.set_min_objective(&NloptOptimizationStep<T>::ObjectiveWrapper, this);

    // Composite algorithm support (AUGLAG or global+local)
    if (algo == nlopt::AUGLAG || algo == nlopt::AUGLAG_EQ) {
      std::string local_algo_str = config_.local_algorithm;
      if (local_algo_str.empty()) {
        local_algo_str = "LN_COBYLA";  // preserve existing default
      }
      const nlopt::algorithm local_algo = ParseNloptAlgorithmForStep(local_algo_str);
      nlopt::opt local_optimizer(local_algo, n_vars);
      local_optimizer.set_xtol_rel(
        config_.local_tolerance > 0.0 ? config_.local_tolerance : tolerance
      );
      int local_max = config_.local_max_iterations;
      if (local_max <= 0) {
        local_max = std::max(10, config_.max_iterations / 10);
      }
      local_optimizer.set_maxeval(local_max);
      optimizer.set_local_optimizer(local_optimizer);
    } else if (!config_.local_algorithm.empty()) {
      // Non-AUGLAG with local algorithm → set as local optimizer for global methods
      const nlopt::algorithm local_algo = ParseNloptAlgorithmForStep(config_.local_algorithm);
      nlopt::opt local_optimizer(local_algo, n_vars);
      local_optimizer.set_xtol_rel(
        config_.local_tolerance > 0.0 ? config_.local_tolerance : tolerance
      );
      int local_max = config_.local_max_iterations;
      if (local_max <= 0) {
        local_max = std::max(10, config_.max_iterations / 10);
      }
      local_optimizer.set_maxeval(local_max);
      optimizer.set_local_optimizer(local_optimizer);
    }

    std::vector<double> lower_bounds(n_vars);
    std::vector<double> upper_bounds(n_vars);
    for (int i = 0; i < n_vars; ++i) {
      lower_bounds[i] = ctx.constraints[i].lower_bound;
      upper_bounds[i] = ctx.constraints[i].upper_bound;
    }

    optimizer.set_maxeval(config_.max_iterations);
    optimizer.set_lower_bounds(lower_bounds);
    optimizer.set_upper_bounds(upper_bounds);

    if (hard_budget) {
      optimizer.add_inequality_constraint(
        &NloptOptimizationStep<T>::BudgetConstraintWrapper,
        this,
        CndOptimizationContext<T>::kConstraintTolerance
      );
    } else if (ctx.verbose) {
      std::cout << "  Selected algorithm does not support nonlinear constraints; "
                   "using soft budget penalty fallback.\n";
    }

    double minf;

    progress_.Start(config_.max_iterations);
    eval_count_ = 0;
    hard_budget_enabled_ = hard_budget;

    try {
      optimizer.optimize(x, minf);
    } catch (...) {
      progress_.Finish();
      ctx_ = nullptr;
      throw;
    }
    progress_.Finish();

    ctx.EnforceBudgetFeasibility(x, lower_bounds);

    for (int i = 0; i < n_vars; ++i) {
      ctx.network.mutable_links()[i].capacity = static_cast<T>(x[i]);
    }

    result.evaluations = eval_count_;
    result.objective = minf;
    result.total_travel_time = ctx.network.TotalTravelTime();
    result.budget = static_cast<double>(ctx.Budget());
    result.success = true;

    ctx_ = nullptr;
    return result;
  }

private:
  OptimizationStepConfig config_;
  CndOptimizationContext<T>* ctx_;
  typename CndOptimizationContext<T>::ProgressState progress_;
  int eval_count_ = 0;
  bool hard_budget_enabled_ = false;

  static double ObjectiveWrapper(unsigned n, const double* x, double* grad, void* data) {
    auto* self = static_cast<NloptOptimizationStep<T>*>(data);
    return self->ObjectiveFunction(n, x, grad);
  }

  double ObjectiveFunction(unsigned n, const double* x, double* grad) {
    (void)grad;
    auto& ctx = *ctx_;

    for (unsigned i = 0; i < n; i++) {
      ctx.network.mutable_links()[i].capacity = x[i];
    }

    double ta_compute_seconds = 0.0;
    const bool ta_valid = ctx.ComputeTrafficFlowsSafely(ta_compute_seconds, "standard_eval");

    double total_travel_time = ctx.network.TotalTravelTime();
    double budget_function = ctx.BudgetFunction(n, x);
    const double budget_violation = std::max(0.0, budget_function - ctx.budget_upper_bound);
    double soft_budget_penalty = 0.0;
    if (!hard_budget_enabled_ && budget_violation > 0.0) {
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
      ctx.PrintInvalidStateWarning("standard_eval", "objective_invalid");
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
    ctx.LogQualityTimePoint("standard_eval",
                            ctx.counters.standard_trace_step,
                            log_objective_function,
                            log_total_travel_time,
                            log_budget_function,
                            ta_compute_seconds);
    CndOptimizationContext<T>::PrintProgressBarWithMetrics(
      progress_, "it/s",
      progress_objective_function,
      progress_total_travel_time,
      progress_budget_function
    );
    return objective_function;
  }

  static double BudgetConstraintWrapper(unsigned n, const double* x, double* grad, void* data) {
    auto* self = static_cast<NloptOptimizationStep<T>*>(data);
    return self->BudgetConstraint(n, x, grad);
  }

  double BudgetConstraint(unsigned n, const double* x, double* grad) {
    auto& ctx = *ctx_;
    if (grad) {
      for (unsigned i = 0; i < n; i++) {
        grad[i] = ctx.budget_function_multiplier;
      }
    }
    double budget = ctx.BudgetFunction(n, x);
    return budget - ctx.budget_upper_bound;
  }
};

inline nlopt::algorithm ParseNloptAlgorithmForStep(const std::string& raw_value) {
  // Convert to uppercase for comparison
  std::string value = raw_value;
  std::transform(value.begin(), value.end(), value.begin(),
    [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

  if (value == "LN_COBYLA") return nlopt::LN_COBYLA;
  if (value == "LN_BOBYQA") return nlopt::LN_BOBYQA;
  if (value == "LN_NELDERMEAD") return nlopt::LN_NELDERMEAD;
  if (value == "LN_SBPLX") return nlopt::LN_SBPLX;
  if (value == "LN_PRAXIS") return nlopt::LN_PRAXIS;
  if (value == "LN_NEWUOA") return nlopt::LN_NEWUOA;
  if (value == "LN_NEWUOA_BOUND") return nlopt::LN_NEWUOA_BOUND;
  if (value == "GN_ISRES") return nlopt::GN_ISRES;
  if (value == "GN_AGS") return nlopt::GN_AGS;
  if (value == "GN_ESCH") return nlopt::GN_ESCH;
  if (value == "GN_CRS2_LM") return nlopt::GN_CRS2_LM;
  if (value == "AUGLAG") return nlopt::AUGLAG;
  if (value == "AUGLAG_EQ") return nlopt::AUGLAG_EQ;
  throw std::runtime_error(
    "Unsupported nlopt algorithm '" + raw_value +
    "'. Supported: LN_COBYLA, LN_BOBYQA, LN_NELDERMEAD, LN_SBPLX, LN_PRAXIS, "
    "LN_NEWUOA, LN_NEWUOA_BOUND, GN_ISRES, GN_AGS, GN_ESCH, GN_CRS2_LM, AUGLAG, AUGLAG_EQ."
  );
}

}  // namespace TrafficAssignment

#endif  // NLOPT_OPTIMIZATION_STEP_H

#ifndef OPTIMLIB_OPTIMIZATION_STEP_H
#define OPTIMLIB_OPTIMIZATION_STEP_H

#include <vector>
#include <string>
#include <algorithm>
#include <chrono>
#include <csignal>
#include <csetjmp>
#include <limits>
#include <stdexcept>
#include <functional>

#ifndef OPTIM_ENABLE_EIGEN_WRAPPERS
#define OPTIM_ENABLE_EIGEN_WRAPPERS
#endif
#include <Eigen/Dense>
#include <optim.hpp>

#include "../OptimizationStep.h"
#include "../CndOptimizationContext.h"

namespace TrafficAssignment {

template <typename T>
class OptimlibOptimizationStep : public OptimizationStep<T> {
public:
  explicit OptimlibOptimizationStep(const OptimizationStepConfig& config)
    : config_(config),
      ctx_(nullptr) {}

  std::string GetName() const override {
    if (!config_.name.empty()) return config_.name;
    return "optimlib(" + config_.algorithm + ")";
  }

  StepResult Execute(CndOptimizationContext<T>& ctx) override {
    StepResult result;
    if (config_.max_iterations == 0) {
      result.success = true;
      return result;
    }

    ctx_ = &ctx;

    const int n_vars = static_cast<int>(ctx.constraints.size());

    // Build initial x from current network capacities
    Eigen::VectorXd x(n_vars);
    for (int i = 0; i < n_vars; ++i) {
      x[i] = static_cast<double>(ctx.network.mutable_links()[i].capacity);
    }

    // Build bounds
    Eigen::VectorXd lower_bounds(n_vars);
    Eigen::VectorXd upper_bounds(n_vars);
    for (int i = 0; i < n_vars; ++i) {
      lower_bounds[i] = ctx.constraints[i].lower_bound;
      upper_bounds[i] = ctx.constraints[i].upper_bound;
    }
    lower_bounds_ = lower_bounds;

    // Configure settings
    optim::algo_settings_t settings;
    settings.vals_bound = true;
    settings.lower_bounds = lower_bounds;
    settings.upper_bounds = upper_bounds;
    settings.iter_max = static_cast<size_t>(config_.max_iterations);
    settings.rel_objfn_change_tol = config_.tolerance;
    settings.print_level = 0;

    const int pop_size = config_.population_size > 0
      ? config_.population_size
      : std::max(200, 2 * n_vars);

    // Narrow initial population range around current capacities to avoid
    // pathological capacity combos that crash the TA solver.
    // Use [0.8*x, 1.2*x] clamped to constraint bounds.
    Eigen::VectorXd init_lb(n_vars), init_ub(n_vars);
    for (int i = 0; i < n_vars; ++i) {
      double xi = x[i];
      double margin = std::max(0.2 * std::abs(xi), 1.0);
      init_lb[i] = std::max(lower_bounds[i], xi - margin);
      init_ub[i] = std::min(upper_bounds[i], xi + margin);
    }

    // Algorithm-specific settings
    settings.de_settings.n_pop = static_cast<size_t>(pop_size);
    settings.de_settings.n_gen = static_cast<size_t>(config_.max_iterations);
    settings.de_settings.initial_lb = init_lb;
    settings.de_settings.initial_ub = init_ub;

    settings.pso_settings.n_pop = static_cast<size_t>(pop_size);
    settings.pso_settings.n_gen = static_cast<size_t>(config_.max_iterations);
    settings.pso_settings.initial_lb = init_lb;
    settings.pso_settings.initial_ub = init_ub;

    settings.nm_settings.adaptive_pars = true;

    // Total expected evaluations for progress bar:
    // DE/PSO: pop_size * max_iterations; NM: just max_iterations
    const int total_evals = IsPopulationBased()
      ? pop_size * config_.max_iterations
      : config_.max_iterations;

    progress_.Start(total_evals);
    eval_count_ = 0;

    auto obj_fn = [](const Eigen::VectorXd& vals_inp,
                     Eigen::VectorXd* grad_out,
                     void* opt_data) -> double {
      auto* self = static_cast<OptimlibOptimizationStep<T>*>(opt_data);
      return self->ObjectiveFunction(vals_inp, grad_out);
    };

    std::string algo = config_.algorithm;
    std::transform(algo.begin(), algo.end(), algo.begin(),
      [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

    bool success = false;
    try {
      if (algo == "DE" || algo == "DIFFERENTIAL_EVOLUTION") {
        success = optim::de(x, obj_fn, this, settings);
      } else if (algo == "DE_PRMM") {
        success = optim::de_prmm(x, obj_fn, this, settings);
      } else if (algo == "PSO" || algo == "PARTICLE_SWARM") {
        success = optim::pso(x, obj_fn, this, settings);
      } else if (algo == "PSO_DV") {
        success = optim::pso_dv(x, obj_fn, this, settings);
      } else if (algo == "NM" || algo == "NELDER_MEAD") {
        success = optim::nm(x, obj_fn, this, settings);
      } else if (algo == "GD" || algo == "GRADIENT_DESCENT") {
        success = optim::gd(x, obj_fn, this, settings);
      } else {
        progress_.Finish();
        ctx_ = nullptr;
        throw std::runtime_error(
          "Unsupported optimlib algorithm '" + config_.algorithm +
          "'. Supported: DE, DE_PRMM, PSO, PSO_DV, NM, GD."
        );
      }
    } catch (...) {
      progress_.Finish();
      ctx_ = nullptr;
      throw;
    }
    progress_.Finish();

    // Post-optimization: enforce budget feasibility and write back
    std::vector<double> x_vec(n_vars);
    std::vector<double> lb_vec(n_vars);
    for (int i = 0; i < n_vars; ++i) {
      x_vec[i] = x[i];
      lb_vec[i] = lower_bounds[i];
    }
    ctx.EnforceBudgetFeasibility(x_vec, lb_vec);

    for (int i = 0; i < n_vars; ++i) {
      ctx.network.mutable_links()[i].capacity = static_cast<T>(x_vec[i]);
    }

    result.evaluations = eval_count_;
    result.objective = settings.opt_fn_value;
    result.total_travel_time = ctx.network.TotalTravelTime();
    result.budget = static_cast<double>(ctx.Budget());
    result.success = success;

    ctx_ = nullptr;
    return result;
  }

private:
  OptimizationStepConfig config_;
  CndOptimizationContext<T>* ctx_;
  typename CndOptimizationContext<T>::ProgressState progress_;
  int eval_count_ = 0;
  Eigen::VectorXd lower_bounds_;
  bool needs_reset_ = false;

  // SIGSEGV recovery for TA solver crashes
  static inline std::jmp_buf s_jump_buf_;
  static inline bool s_jump_valid_ = false;

  static void SignalHandler(int sig) {
    if (sig == SIGSEGV && s_jump_valid_) {
      s_jump_valid_ = false;
      std::longjmp(s_jump_buf_, 1);
    }
    // Re-raise if we can't handle it
    std::signal(sig, SIG_DFL);
    std::raise(sig);
  }

  bool IsPopulationBased() const {
    std::string algo = config_.algorithm;
    std::transform(algo.begin(), algo.end(), algo.begin(),
      [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return algo == "DE" || algo == "DE_PRMM" ||
           algo == "DIFFERENTIAL_EVOLUTION" ||
           algo == "PSO" || algo == "PSO_DV" ||
           algo == "PARTICLE_SWARM";
  }

  double ObjectiveFunction(const Eigen::VectorXd& vals_inp,
                           Eigen::VectorXd* /*grad_out*/) {
    auto& ctx = *ctx_;
    const int n = static_cast<int>(vals_inp.size());

    // Clamp values to bounds — population-based optimizers may explore
    // outside bounds, and extreme capacities can crash the TA solver.
    static constexpr double kMinCapacity = 1e-6;
    for (int i = 0; i < n; i++) {
      double val = vals_inp[i];
      val = std::max(val, ctx.constraints[i].lower_bound);
      val = std::min(val, ctx.constraints[i].upper_bound);
      val = std::max(val, kMinCapacity);
      ctx.network.mutable_links()[i].capacity = val;
    }

    double ta_compute_seconds = 0.0;
    // Population-based methods explore far from the feasible region.
    // The TA solver (especially Tapas) can segfault with pathological
    // capacity combinations. Use signal handling to catch SIGSEGV.
    bool ta_valid = false;
    if (needs_reset_) {
      for (int i = 0; i < ctx.network.number_of_links(); ++i) {
        ctx.network.mutable_links()[i].flow = T(0);
      }
      try { ctx.approach->Reset(); } catch (...) {}
      needs_reset_ = false;
    }

    // Install SIGSEGV handler around TA computation
    s_jump_valid_ = true;
    auto prev_handler = std::signal(SIGSEGV, &OptimlibOptimizationStep::SignalHandler);
    if (setjmp(s_jump_buf_) == 0) {
      const auto ta_start = std::chrono::steady_clock::now();
      ctx.approach->ComputeTrafficFlows();
      ta_compute_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - ta_start).count();
      ta_valid = !ctx.HasInvalidTrafficAssignmentState();
    } else {
      // Recovered from SIGSEGV
      ta_valid = false;
    }
    s_jump_valid_ = false;
    std::signal(SIGSEGV, prev_handler);

    if (!ta_valid) {
      ++ctx.counters.invalid_ta_state_count;
      needs_reset_ = true;
    }

    double total_travel_time = ctx.network.TotalTravelTime();
    double budget_function = ctx.BudgetFunction(n, vals_inp.data());
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
      progress_, "eval/s",
      progress_objective_function,
      progress_total_travel_time,
      progress_budget_function
    );
    return objective_function;
  }
};

}  // namespace TrafficAssignment

#endif  // OPTIMLIB_OPTIMIZATION_STEP_H

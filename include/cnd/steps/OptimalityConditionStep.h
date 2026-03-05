#ifndef OPTIMALITY_CONDITION_STEP_H
#define OPTIMALITY_CONDITION_STEP_H

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include <Eigen/Dense>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "../OptimizationStep.h"
#include "../CndOptimizationContext.h"
#include "../../tap/data/Link.h"

namespace TrafficAssignment {

namespace mp_internal = boost::multiprecision;
using high_prec_internal = mp_internal::cpp_dec_float_50;

template <typename T>
class OptimalityConditionStep : public OptimizationStep<T> {
  using MatrixXd = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
  using MatrixHightPrecision = Eigen::Matrix<high_prec_internal, Eigen::Dynamic, Eigen::Dynamic>;

public:
  explicit OptimalityConditionStep(const OptimizationStepConfig& config)
    : config_(config) {}

  std::string GetName() const override {
    if (!config_.name.empty()) return config_.name;
    return "optimality_condition";
  }

  StepResult Execute(CndOptimizationContext<T>& ctx) override {
    StepResult result;
    if (config_.max_iterations == 0) {
      result.success = true;
      return result;
    }

    typename CndOptimizationContext<T>::ProgressState progress;
    progress.Start(config_.max_iterations);

    try {
      for (int cnt = 0; cnt < config_.max_iterations; cnt++) {
        if (ctx.Budget() < ctx.budget_upper_bound - ctx.budget_threshold) {
          NonUpperBoundBudgetOptimalityIteration(ctx);
        } else {
          UpperBoundBudgetOptimalityIteration(ctx);
        }

        double ta_compute_seconds = 0.0;
        if (!ctx.ComputeTrafficFlowsSafely(ta_compute_seconds, "optcond_iter")) {
          throw std::runtime_error("TA computation failed during optimality condition iteration.");
        }
        double total_travel_time = ctx.network.TotalTravelTime();
        double budget_function = static_cast<double>(ctx.Budget());
        double objective_function = total_travel_time + budget_function;
        double log_total_travel_time = total_travel_time;
        double log_budget_function = budget_function;
        double log_objective_function = objective_function;
        double progress_total_travel_time = total_travel_time;
        double progress_budget_function = budget_function;
        double progress_objective_function = objective_function;
        if (!CndOptimizationContext<T>::IsFiniteNumber(total_travel_time) ||
            !CndOptimizationContext<T>::IsFiniteNumber(budget_function) ||
            !CndOptimizationContext<T>::IsFiniteNumber(objective_function)) {
          log_total_travel_time = std::numeric_limits<double>::quiet_NaN();
          log_budget_function = std::numeric_limits<double>::quiet_NaN();
          log_objective_function = std::numeric_limits<double>::quiet_NaN();
          progress_total_travel_time = CndOptimizationContext<T>::kInvalidObjectivePenalty;
          progress_budget_function =
            CndOptimizationContext<T>::IsFiniteNumber(budget_function) ? budget_function : 0.0;
          progress_objective_function = CndOptimizationContext<T>::kInvalidObjectivePenalty;
          ctx.PrintInvalidStateWarning("optcond_iter", "objective_invalid");
        }
        ctx.counters.optcond_trace_step = cnt + 1;
        ctx.UpdateRuntimeAndBudgetMetrics(ta_compute_seconds, progress_budget_function);
        ctx.LogQualityTimePoint("optcond_iter",
                                ctx.counters.optcond_trace_step,
                                log_objective_function,
                                log_total_travel_time,
                                log_budget_function,
                                ta_compute_seconds);
        CndOptimizationContext<T>::PrintProgressBarWithMetrics(
          progress, "it/s",
          progress_objective_function,
          progress_total_travel_time,
          progress_budget_function
        );
      }
    } catch (...) {
      progress.Finish();
      throw;
    }
    progress.Finish();

    result.total_travel_time = ctx.network.TotalTravelTime();
    result.budget = static_cast<double>(ctx.Budget());
    result.objective = result.total_travel_time + result.budget;
    result.evaluations = config_.max_iterations;
    result.success = true;
    return result;
  }

private:
  OptimizationStepConfig config_;

  // --- Optimality math ---

  struct OdCache {
    std::vector<std::vector<int>> routes;
    MatrixHightPrecision          jacobi_inverse;
    MatrixHightPrecision          e_col;
    high_prec_internal            denominator;
    T                             demand;
    bool                          valid = false;
  };

  MatrixXd CapacityDerColumn(CndOptimizationContext<T>& ctx,
                             const std::vector<std::vector<int>>& routes,
                             int link_index) {
    MatrixXd res = MatrixXd::Constant(routes.size(), 1, 0.0);
    for (std::size_t i = 0; i < routes.size(); i++) {
      if (std::find(routes[i].begin(), routes[i].end(), link_index) != routes[i].end()) {
        res(i, 0) = ctx.network.links()[link_index].DelayCapacityDer();
      }
    }
    return res;
  }

  MatrixHightPrecision ConvertEigenMatrix(const MatrixXd& source) {
    MatrixHightPrecision target(source.rows(), source.cols());
    for (int i = 0; i < source.rows(); ++i) {
      for (int j = 0; j < source.cols(); ++j) {
        target(i, j) = static_cast<high_prec_internal>(source(i, j));
      }
    }
    return target;
  }

  T OptimalityFunction(CndOptimizationContext<T>& ctx,
                        int link_index,
                        bool flow_computation = true) {
    if (flow_computation) {
      double ta_compute_seconds = 0.0;
      if (!ctx.ComputeTrafficFlowsSafely(ta_compute_seconds, "optimality_function")) {
        return ctx.NaNValue();
      }
    }
    T result = 0;
    for (int od_pair_index = 0; od_pair_index < ctx.network.number_of_od_pairs(); od_pair_index++) {
      const int routes_count = ctx.network.mutable_od_pairs()[od_pair_index].GetRoutesCount();
      if (routes_count <= 0) continue;
      auto routes = ctx.network.mutable_od_pairs()[od_pair_index].GetRoutes();
      auto e_column = MatrixHightPrecision::Constant(routes_count, 1, 1.0);
      auto jacobi = ConvertEigenMatrix(RoutesJacobiMatrix(routes, ctx.network.mutable_links()));
      auto capacity_der_column = ConvertEigenMatrix(CapacityDerColumn(ctx, routes, link_index));
      auto jacobi_inverse = jacobi.inverse();
      const auto numerator = (e_column.transpose() * jacobi_inverse * capacity_der_column)(0, 0);
      const auto denominator = (e_column.transpose() * jacobi_inverse * e_column)(0, 0);
      if (mp_internal::abs(denominator) <
          static_cast<high_prec_internal>(CndOptimizationContext<T>::kDenominatorZeroGuard)) {
        return ctx.NaNValue();
      }
      const T local_value = static_cast<T>(-numerator / denominator);
      if (!ctx.IsFiniteScalar(local_value)) return ctx.NaNValue();
      result += ctx.network.od_pairs()[od_pair_index].GetDemand() * local_value;
    }
    result /= ctx.budget_function_multiplier;
    if (!ctx.IsFiniteScalar(result)) return ctx.NaNValue();
    return result;
  }

  std::vector<OdCache> BuildOdCache(CndOptimizationContext<T>& ctx) {
    const int n_od = ctx.network.number_of_od_pairs();
    std::vector<OdCache> od_cache(n_od);
    for (int od = 0; od < n_od; ++od) {
      const int rc = ctx.network.mutable_od_pairs()[od].GetRoutesCount();
      if (rc <= 0) continue;
      OdCache& c = od_cache[od];
      c.routes  = ctx.network.mutable_od_pairs()[od].GetRoutes();
      c.demand  = ctx.network.od_pairs()[od].GetDemand();
      c.e_col   = MatrixHightPrecision::Constant(rc, 1, high_prec_internal(1));
      auto jacobi_hp = ConvertEigenMatrix(
          RoutesJacobiMatrix(c.routes, ctx.network.mutable_links()));
      c.jacobi_inverse = jacobi_hp.inverse();
      c.denominator    = (c.e_col.transpose() * c.jacobi_inverse * c.e_col)(0, 0);
      if (mp_internal::abs(c.denominator) <
          static_cast<high_prec_internal>(CndOptimizationContext<T>::kDenominatorZeroGuard)) {
        continue;
      }
      c.valid = true;
    }
    return od_cache;
  }

  T OptimalityFunctionFromCache(CndOptimizationContext<T>& ctx,
                                 int link_index,
                                 const std::vector<OdCache>& od_cache) {
    T result = T(0);
    for (const OdCache& c : od_cache) {
      if (!c.valid) continue;
      auto cap_der = ConvertEigenMatrix(CapacityDerColumn(ctx, c.routes, link_index));
      const high_prec_internal numerator =
          (c.e_col.transpose() * c.jacobi_inverse * cap_der)(0, 0);
      const T local_value = static_cast<T>(-numerator / c.denominator);
      if (!ctx.IsFiniteScalar(local_value)) return ctx.NaNValue();
      result += c.demand * local_value;
    }
    result /= ctx.budget_function_multiplier;
    return ctx.IsFiniteScalar(result) ? result : ctx.NaNValue();
  }

  // --- Upper-bound optimality iteration ---

  bool TransferResult(CndOptimizationContext<T>& ctx,
                      T amount, int max_optimality_index, int min_optimality_index) {
    ctx.network.mutable_links()[max_optimality_index].capacity += amount;
    ctx.network.mutable_links()[min_optimality_index].capacity -= amount;
    const T max_value = OptimalityFunction(ctx, max_optimality_index);
    const T min_value = OptimalityFunction(ctx, min_optimality_index);
    bool result = ctx.IsFiniteScalar(max_value) && ctx.IsFiniteScalar(min_value) && max_value > min_value;
    ctx.network.mutable_links()[max_optimality_index].capacity -= amount;
    ctx.network.mutable_links()[min_optimality_index].capacity += amount;
    return result;
  }

  T RespectedAmountConstraints(CndOptimizationContext<T>& ctx,
                                T amount, int max_optimality_index, int min_optimality_index) {
    auto max_optimality_capacity = ctx.network.mutable_links()[max_optimality_index].capacity;
    auto min_optimality_capacity = ctx.network.mutable_links()[min_optimality_index].capacity;
    amount = std::min(amount, static_cast<T>(ctx.constraints[max_optimality_index].upper_bound) - max_optimality_capacity);
    amount = std::min(amount, min_optimality_capacity - static_cast<T>(ctx.constraints[min_optimality_index].lower_bound));
    return amount;
  }

  T OptimalityConditionTransferAmount(CndOptimizationContext<T>& ctx,
                                       int max_optimality_index, int min_optimality_index) {
    T amount = static_cast<T>(CndOptimizationContext<T>::kInitialTransferAmount);
    amount = RespectedAmountConstraints(ctx, amount, max_optimality_index, min_optimality_index);
    T next_amount = RespectedAmountConstraints(
      ctx, amount * CndOptimizationContext<T>::kBinarySearchDivisor,
      max_optimality_index, min_optimality_index);

    while (TransferResult(ctx, next_amount, max_optimality_index, min_optimality_index)) {
      amount = next_amount;
      next_amount = RespectedAmountConstraints(
        ctx, amount * CndOptimizationContext<T>::kBinarySearchDivisor,
        max_optimality_index, min_optimality_index);
      if (std::abs(next_amount - amount) < ctx.link_capacity_selection_threshold / 10) {
        break;
      }
    }
    while (!TransferResult(ctx, amount, max_optimality_index, min_optimality_index)) {
      amount /= CndOptimizationContext<T>::kBinarySearchDivisor;
    }
    return amount;
  }

  void UpperBoundBudgetOptimalityIteration(CndOptimizationContext<T>& ctx) {
    int max_index = -1, min_index = -1;

    for (int link_index = 0; link_index < ctx.network.number_of_links(); link_index++) {
      if (ctx.network.mutable_links()[link_index].flow <
          CndOptimizationContext<T>::kInitialTransferAmount) continue;
      auto cur_capacity = ctx.network.mutable_links()[link_index].capacity;
      if (std::abs(ctx.constraints[link_index].upper_bound - static_cast<double>(cur_capacity)) >
              static_cast<double>(ctx.link_capacity_selection_threshold) &&
          std::abs(ctx.constraints[link_index].lower_bound - static_cast<double>(cur_capacity)) >
              static_cast<double>(ctx.link_capacity_selection_threshold)) {
        max_index = link_index;
        min_index = link_index;
        break;
      }
    }
    if (max_index == -1) return;

    double ta_compute_seconds = 0.0;
    if (!ctx.ComputeTrafficFlowsSafely(ta_compute_seconds, "upper_bound_opt_iteration")) return;

    const auto od_cache = BuildOdCache(ctx);

    T max_value = OptimalityFunctionFromCache(ctx, max_index, od_cache);
    if (!ctx.IsFiniteScalar(max_value)) return;
    T min_value = max_value;

    for (int link_index = 0; link_index < ctx.network.number_of_links(); link_index++) {
      if (ctx.network.mutable_links()[link_index].flow <
          CndOptimizationContext<T>::kInitialTransferAmount) continue;
      auto cur_capacity = ctx.network.mutable_links()[link_index].capacity;
      auto cur_value = OptimalityFunctionFromCache(ctx, link_index, od_cache);
      if (!ctx.IsFiniteScalar(cur_value)) continue;
      if (cur_value > max_value &&
          std::abs(ctx.constraints[link_index].upper_bound - static_cast<double>(cur_capacity)) >
              static_cast<double>(ctx.link_capacity_selection_threshold)) {
        max_index = link_index;
        max_value = cur_value;
      }
      if (cur_value < min_value &&
          std::abs(ctx.constraints[link_index].lower_bound - static_cast<double>(cur_capacity)) >
              static_cast<double>(ctx.link_capacity_selection_threshold)) {
        min_index = link_index;
        min_value = cur_value;
      }
    }

    T amount = OptimalityConditionTransferAmount(ctx, max_index, min_index);

    ctx.network.mutable_links()[max_index].capacity += amount;
    ctx.network.mutable_links()[min_index].capacity -= amount;
  }

  // --- Non-upper-bound optimality iteration ---

  bool ProcessingCheck(CndOptimizationContext<T>& ctx,
                       int link_index, T middle_bound_value,
                       const std::vector<OdCache>& od_cache) {
    if (ctx.network.mutable_links()[link_index].flow <
        CndOptimizationContext<T>::kFlowSignificanceThreshold) return false;
    auto cur_capacity = ctx.network.mutable_links()[link_index].capacity;
    auto cur_value = OptimalityFunctionFromCache(ctx, link_index, od_cache);
    if (!ctx.IsFiniteScalar(cur_value)) return false;
    if (std::abs(ctx.constraints[link_index].lower_bound - static_cast<double>(cur_capacity))
            < static_cast<double>(ctx.link_capacity_selection_threshold) &&
        cur_value <= middle_bound_value)
      return false;
    if (std::abs(ctx.constraints[link_index].upper_bound - static_cast<double>(cur_capacity))
            < static_cast<double>(ctx.link_capacity_selection_threshold) &&
        cur_value >= middle_bound_value)
      return false;
    return true;
  }

  T NonUpperBoundRespectedAmountConstraints(CndOptimizationContext<T>& ctx,
                                             T amount, int max_optimality_index) {
    auto max_optimality_capacity = ctx.network.mutable_links()[max_optimality_index].capacity;
    if (amount < 0) {
      return std::max(amount,
        static_cast<T>(ctx.constraints[max_optimality_index].lower_bound) - max_optimality_capacity);
    } else {
      amount = std::min(amount,
        static_cast<T>(ctx.constraints[max_optimality_index].upper_bound) - max_optimality_capacity);
      amount = std::min(amount,
        static_cast<T>((ctx.budget_upper_bound - static_cast<double>(ctx.Budget())) /
                        static_cast<double>(ctx.budget_function_multiplier)));
    }
    return amount;
  }

  bool NonUpperBoundTransferResult(CndOptimizationContext<T>& ctx,
                                    T amount, int max_optimality_index, T middle_bound_value) {
    ctx.network.mutable_links()[max_optimality_index].capacity += amount;
    const T optimality_value = OptimalityFunction(ctx, max_optimality_index);
    if (!ctx.IsFiniteScalar(optimality_value)) {
      ctx.network.mutable_links()[max_optimality_index].capacity -= amount;
      return false;
    }
    bool result = false;
    if (amount < 0) {
      result = optimality_value < middle_bound_value;
    } else {
      result = optimality_value > middle_bound_value;
    }
    ctx.network.mutable_links()[max_optimality_index].capacity -= amount;
    return result;
  }

  T NonUpperBoundOptimalityConditionTransferAmount(CndOptimizationContext<T>& ctx,
                                                    int max_optimality_index,
                                                    T middle_bound_value) {
    T amount = static_cast<T>(CndOptimizationContext<T>::kInitialTransferAmount);
    T next_amount = 0;
    auto max_optimality_capacity = ctx.network.mutable_links()[max_optimality_index].capacity;
    amount = NonUpperBoundRespectedAmountConstraints(ctx, amount, max_optimality_index);
    auto cur_value = OptimalityFunction(ctx, max_optimality_index, false);
    if (!ctx.IsFiniteScalar(cur_value)) return 0;

    if (cur_value < middle_bound_value) {
      next_amount = static_cast<T>(ctx.constraints[max_optimality_index].lower_bound) -
                    max_optimality_capacity;
      if (NonUpperBoundTransferResult(ctx, next_amount, max_optimality_index, middle_bound_value)) {
        return next_amount;
      }
      amount = static_cast<T>(-CndOptimizationContext<T>::kInitialTransferAmount);
      amount = NonUpperBoundRespectedAmountConstraints(ctx, amount, max_optimality_index);
      next_amount = NonUpperBoundRespectedAmountConstraints(
        ctx, amount * -CndOptimizationContext<T>::kBinarySearchDivisor, max_optimality_index);
      while (NonUpperBoundTransferResult(ctx, next_amount, max_optimality_index, middle_bound_value)) {
        amount = next_amount;
        next_amount = NonUpperBoundRespectedAmountConstraints(
          ctx, amount * -CndOptimizationContext<T>::kBinarySearchDivisor, max_optimality_index);
        if (std::abs(next_amount - amount) < ctx.link_capacity_selection_threshold / 10) {
          break;
        }
      }
    } else {
      next_amount = NonUpperBoundRespectedAmountConstraints(
        ctx,
        static_cast<T>(ctx.constraints[max_optimality_index].upper_bound) - max_optimality_capacity,
        max_optimality_index);
      while (!NonUpperBoundTransferResult(ctx, next_amount, max_optimality_index, middle_bound_value)) {
        next_amount /= CndOptimizationContext<T>::kBinarySearchDivisor;
      }
      amount = next_amount;
    }
    return amount;
  }

  void NonUpperBoundBudgetOptimalityIteration(CndOptimizationContext<T>& ctx) {
    T middle_bound_value = 1.0;
    int max_optimality_index = -1;

    double ta_compute_seconds = 0.0;
    if (!ctx.ComputeTrafficFlowsSafely(ta_compute_seconds, "non_upper_bound_opt_iteration")) return;

    const auto od_cache = BuildOdCache(ctx);

    for (int link_index = 0; link_index < ctx.network.number_of_links(); link_index++) {
      if (ProcessingCheck(ctx, link_index, middle_bound_value, od_cache)) {
        max_optimality_index = link_index;
        break;
      }
    }
    if (max_optimality_index == -1) return;

    const T base_value = OptimalityFunctionFromCache(ctx, max_optimality_index, od_cache);
    if (!ctx.IsFiniteScalar(base_value)) return;
    T max_delta = std::abs(middle_bound_value - base_value);
    for (int link_index = 0; link_index < ctx.network.number_of_links(); link_index++) {
      if (ProcessingCheck(ctx, link_index, middle_bound_value, od_cache)) {
        const T optimality_value = OptimalityFunctionFromCache(ctx, link_index, od_cache);
        if (!ctx.IsFiniteScalar(optimality_value)) continue;
        auto cur_delta = std::abs(middle_bound_value - optimality_value);
        if (cur_delta > max_delta) {
          max_optimality_index = link_index;
          max_delta = cur_delta;
        }
      }
    }
    auto amount = NonUpperBoundOptimalityConditionTransferAmount(
      ctx, max_optimality_index, middle_bound_value);
    ctx.network.mutable_links()[max_optimality_index].capacity += amount;
  }
};

}  // namespace TrafficAssignment

#endif  // OPTIMALITY_CONDITION_STEP_H

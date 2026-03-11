#ifndef OPTIMALITY_CONDITION_STEP_H
#define OPTIMALITY_CONDITION_STEP_H

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include <unordered_set>
#include <Eigen/Dense>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "../OptimizationStep.h"
#include "../CndOptimizationContext.h"
#include "../../tap/data/Link.h"

namespace TrafficAssignment {

namespace mp_internal = boost::multiprecision;
using high_prec_internal = mp_internal::cpp_dec_float_50;

/**
 * @brief Sensitivity-analysis-based optimization step using first-order optimality conditions.
 *
 * Implements the iterative capacity update procedure from Chiou (2005) Section 3.
 * Each iteration:
 * 1. Solves the lower-level TAP for current capacities y to get equilibrium flows x(y).
 * 2. Builds a Jacobian-based OD cache for efficient sensitivity computation.
 * 3. Computes the optimality function phi_a for each design link (Chiou 2005, Eq. 9):
 *      phi_a = dF/dy_a = sum_w d_w * (eT J^{-1} dc_a/dy_a) / (eT J^{-1} e) / theta
 *    where J is the route cost Jacobian, d_w is OD demand, and e is the all-ones vector.
 * 4. Transfers capacity between links with max/min phi_a to equalize optimality values.
 *
 * Two iteration modes based on budget usage:
 * - Budget exhausted (B(y) >= B_max - threshold): redistributes capacity within fixed budget
 *   by transferring from max-phi link to min-phi link (budget-neutral).
 * - Budget remaining (B(y) < B_max - threshold): expands the link with largest |phi_a - 1|,
 *   moving it toward the common optimality value of 1.0.
 *
 * Uses Boost.Multiprecision (cpp_dec_float_50) for Jacobian inversion to avoid
 * numerical issues with near-singular route cost matrices.
 *
 * @tparam T Numeric type for flow/capacity computations.
 */
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

  /// @brief Runs iterative capacity updates using optimality conditions. Chiou (2005) Section 3.
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

  /**
   * @brief Cached per-OD-pair data for efficient optimality function evaluation.
   *
   * Precomputes the route cost Jacobian inverse, the all-ones column eT*J^{-1}*e denominator,
   * and route-link membership sets. This avoids redundant matrix inversions when evaluating
   * phi_a across all design links.
   */
  struct OdCache {
    std::vector<std::vector<int>> routes;            ///< Route definitions (list of link IDs per route).
    MatrixHightPrecision          jacobi_inverse;    ///< J^{-1}: inverse of route cost Jacobian (high precision).
    MatrixHightPrecision          e_col;             ///< e: all-ones column vector (R x 1).
    MatrixHightPrecision          inv_transpose_e;   ///< J^{-1}^T * e: precomputed for O(R) dot product.
    high_prec_internal            denominator;       ///< eT * J^{-1} * e: shared denominator for all links.
    T                             demand;            ///< OD pair demand d_w.
    std::vector<std::unordered_set<int>> route_link_sets; ///< route_link_sets[r] = {link IDs in route r} for O(1) lookup.
    bool                          valid = false;     ///< False if Jacobian is near-singular.
  };

  /// @brief Builds the dc_a/dy_a column: partial derivative of route costs w.r.t. capacity of link_index.
  ///        Non-zero only for routes containing the link. Chiou (2005) Eq. 6.
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

  /// @brief Fast dc_a/dy_a column using precomputed route-link membership sets from OdCache.
  MatrixXd CapacityDerColumnCached(CndOptimizationContext<T>& ctx,
                                    const OdCache& cache,
                                    int link_index) {
    MatrixXd res = MatrixXd::Constant(cache.routes.size(), 1, 0.0);
    T delay_cap_der = ctx.network.links()[link_index].DelayCapacityDer();
    for (std::size_t i = 0; i < cache.routes.size(); i++) {
      if (cache.route_link_sets[i].count(link_index)) {
        res(i, 0) = delay_cap_der;
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

  /**
   * @brief Computes the optimality function phi_a for a single link. Chiou (2005) Eq. 9.
   *
   * phi_a = (1/theta) * sum_w d_w * (eT J_w^{-1} dc_a) / (eT J_w^{-1} e)
   *
   * At optimality, phi_a should be equal across all active design links (KKT condition).
   * This is the non-cached version that rebuilds the Jacobian from scratch.
   */
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

  /// @brief Builds the OD cache for all OD pairs: precomputes Jacobian inverses and denominators.
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

      // Build route-link membership sets for O(1) lookup
      c.route_link_sets.resize(rc);
      for (int r = 0; r < rc; ++r) {
        c.route_link_sets[r].insert(c.routes[r].begin(), c.routes[r].end());
      }

      auto jacobi_hp = ConvertEigenMatrix(
          RoutesJacobiMatrix(c.routes, ctx.network.mutable_links()));
      c.jacobi_inverse = jacobi_hp.inverse();
      c.denominator    = (c.e_col.transpose() * c.jacobi_inverse * c.e_col)(0, 0);
      if (mp_internal::abs(c.denominator) <
          static_cast<high_prec_internal>(CndOptimizationContext<T>::kDenominatorZeroGuard)) {
        continue;
      }
      // Precompute J⁻¹ᵀ * e for O(R) dot product in OptimalityFunctionFromCache
      c.inv_transpose_e = c.jacobi_inverse.transpose() * c.e_col;
      c.valid = true;
    }
    return od_cache;
  }

  /// @brief Computes phi_a using precomputed OdCache (avoids redundant Jacobian inversions).
  T OptimalityFunctionFromCache(CndOptimizationContext<T>& ctx,
                                 int link_index,
                                 const std::vector<OdCache>& od_cache) {
    T result = T(0);
    for (const OdCache& c : od_cache) {
      if (!c.valid) continue;
      auto cap_der = ConvertEigenMatrix(CapacityDerColumnCached(ctx, c, link_index));
      // Use precomputed inv_transpose_e for O(R) dot product instead of O(R²) matrix-vector multiply
      // numerator = eᵀ * J⁻¹ * cap_der = (J⁻¹ᵀ * e)ᵀ * cap_der = inv_transpose_eᵀ * cap_der
      const high_prec_internal numerator =
          (c.inv_transpose_e.transpose() * cap_der)(0, 0);
      const T local_value = static_cast<T>(-numerator / c.denominator);
      if (!ctx.IsFiniteScalar(local_value)) return ctx.NaNValue();
      result += c.demand * local_value;
    }
    result /= ctx.budget_function_multiplier;
    return ctx.IsFiniteScalar(result) ? result : ctx.NaNValue();
  }

  // --- Upper-bound budget optimality iteration (budget exhausted) ---

  /// @brief Clamps transfer amount to respect both links' capacity box constraints [lb, ub].
  T RespectedAmountConstraints(CndOptimizationContext<T>& ctx,
                                T amount, int max_optimality_index, int min_optimality_index) {
    auto max_optimality_capacity = ctx.network.mutable_links()[max_optimality_index].capacity;
    auto min_optimality_capacity = ctx.network.mutable_links()[min_optimality_index].capacity;
    amount = std::min(amount, static_cast<T>(ctx.constraints[max_optimality_index].upper_bound) - max_optimality_capacity);
    amount = std::min(amount, min_optimality_capacity - static_cast<T>(ctx.constraints[min_optimality_index].lower_bound));
    return amount;
  }

  /**
   * @brief Estimates capacity transfer amount using first-order sensitivity.
   *
   * Uses the gap between max/min optimality function values to scale the transfer
   * conservatively, rather than binary search with repeated TA solves.
   */
  T ApproximateTransferAmount(CndOptimizationContext<T>& ctx,
                               int max_optimality_index, int min_optimality_index,
                               T max_value, T min_value) {
    // The gap between optimality function values indicates the gradient magnitude.
    // Use a fraction of the constraint-bounded maximum as the transfer amount.
    // Start with initial amount and respect constraints.
    T max_amount = static_cast<T>(CndOptimizationContext<T>::kInitialTransferAmount);
    max_amount = RespectedAmountConstraints(ctx, max_amount, max_optimality_index, min_optimality_index);

    if (max_amount <= 0) return T(0);

    // Use the optimality gap to scale: larger gap → more confident in larger transfer
    T gap = max_value - min_value;
    if (gap <= 0) return max_amount;

    // Scale factor based on gap relative to absolute values
    T avg_abs = (std::abs(max_value) + std::abs(min_value)) / 2;
    T scale = (avg_abs > 0) ? std::min(gap / avg_abs, T(1.0)) : T(1.0);

    // Conservative step: use scale * max_amount, bounded below by threshold
    T amount = scale * max_amount;
    amount = std::max(amount, ctx.link_capacity_selection_threshold);
    amount = RespectedAmountConstraints(ctx, amount, max_optimality_index, min_optimality_index);

    return std::max(amount, T(0));
  }

  /// @brief Validates transfer with one TA solve; halves amount on failure (up to kMaxValidationRetries).
  T ValidatedTransferAmount(CndOptimizationContext<T>& ctx,
                              int max_optimality_index, int min_optimality_index,
                              T initial_amount) {
    static constexpr int kMaxValidationRetries = 5;
    T amount = initial_amount;

    for (int retry = 0; retry < kMaxValidationRetries && amount > ctx.link_capacity_selection_threshold / 10; ++retry) {
      ctx.network.mutable_links()[max_optimality_index].capacity += amount;
      ctx.network.mutable_links()[min_optimality_index].capacity -= amount;

      double ta_seconds = 0.0;
      bool ta_ok = ctx.ComputeTrafficFlowsSafely(ta_seconds, "transfer_validate");

      if (ta_ok) {
        // Build cache for the modified state and verify ordering still holds
        const auto validate_cache = BuildOdCache(ctx);
        T new_max = OptimalityFunctionFromCache(ctx, max_optimality_index, validate_cache);
        T new_min = OptimalityFunctionFromCache(ctx, min_optimality_index, validate_cache);

        ctx.network.mutable_links()[max_optimality_index].capacity -= amount;
        ctx.network.mutable_links()[min_optimality_index].capacity += amount;

        if (ctx.IsFiniteScalar(new_max) && ctx.IsFiniteScalar(new_min) && new_max > new_min) {
          return amount;
        }
      } else {
        ctx.network.mutable_links()[max_optimality_index].capacity -= amount;
        ctx.network.mutable_links()[min_optimality_index].capacity += amount;
      }

      // Halve and retry
      amount /= CndOptimizationContext<T>::kBinarySearchDivisor;
    }
    return amount;
  }

  /**
   * @brief Budget-exhausted iteration: transfers capacity from max-phi to min-phi link.
   *
   * Finds the links with maximum and minimum optimality function values among active links
   * (flow > threshold, capacity not at bounds), then transfers capacity to equalize phi.
   * Budget-neutral since capacity is redistributed between two links.
   */
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

    T amount = ApproximateTransferAmount(ctx, max_index, min_index, max_value, min_value);
    amount = ValidatedTransferAmount(ctx, max_index, min_index, amount);

    ctx.network.mutable_links()[max_index].capacity += amount;
    ctx.network.mutable_links()[min_index].capacity -= amount;
  }

  // --- Non-upper-bound budget optimality iteration (budget remaining) ---

  /// @brief Checks if a link is eligible for capacity adjustment (active flow, not at bound in wrong direction).
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

  /// @brief Clamps transfer amount to box constraints and remaining budget.
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

  /// @brief Estimates single-link capacity change to move phi_a toward the target value (1.0).
  T ApproximateNonUpperBoundTransferAmount(CndOptimizationContext<T>& ctx,
                                            int max_optimality_index,
                                            T cur_value, T middle_bound_value) {
    auto max_optimality_capacity = ctx.network.mutable_links()[max_optimality_index].capacity;

    if (cur_value < middle_bound_value) {
      // Need to decrease capacity (move toward lower bound)
      T amount = static_cast<T>(-CndOptimizationContext<T>::kInitialTransferAmount);
      amount = NonUpperBoundRespectedAmountConstraints(ctx, amount, max_optimality_index);

      // Scale by the relative gap
      T gap = middle_bound_value - cur_value;
      T avg_abs = (std::abs(cur_value) + std::abs(middle_bound_value)) / 2;
      T scale = (avg_abs > 0) ? std::min(gap / avg_abs, T(1.0)) : T(1.0);
      amount = scale * amount;  // amount is negative, scale makes it less negative (conservative)

      amount = NonUpperBoundRespectedAmountConstraints(ctx, amount, max_optimality_index);
      return amount;
    } else {
      // Need to increase capacity (move toward upper bound)
      T max_amount = static_cast<T>(ctx.constraints[max_optimality_index].upper_bound) -
                     max_optimality_capacity;
      max_amount = NonUpperBoundRespectedAmountConstraints(ctx, max_amount, max_optimality_index);

      T gap = cur_value - middle_bound_value;
      T avg_abs = (std::abs(cur_value) + std::abs(middle_bound_value)) / 2;
      T scale = (avg_abs > 0) ? std::min(gap / avg_abs, T(1.0)) : T(1.0);

      T amount = scale * max_amount;
      amount = std::max(amount, ctx.link_capacity_selection_threshold);
      amount = NonUpperBoundRespectedAmountConstraints(ctx, amount, max_optimality_index);
      return amount;
    }
  }

  /// @brief Validates single-link transfer with one TA solve; halves on failure.
  T ValidatedNonUpperBoundTransferAmount(CndOptimizationContext<T>& ctx,
                                          int max_optimality_index,
                                          T middle_bound_value,
                                          T initial_amount) {
    static constexpr int kMaxValidationRetries = 5;
    T amount = initial_amount;

    for (int retry = 0; retry < kMaxValidationRetries &&
         std::abs(amount) > ctx.link_capacity_selection_threshold / 10; ++retry) {
      ctx.network.mutable_links()[max_optimality_index].capacity += amount;

      double ta_seconds = 0.0;
      bool ta_ok = ctx.ComputeTrafficFlowsSafely(ta_seconds, "non_ub_transfer_validate");

      if (ta_ok) {
        const auto validate_cache = BuildOdCache(ctx);
        T optimality_value = OptimalityFunctionFromCache(ctx, max_optimality_index, validate_cache);

        ctx.network.mutable_links()[max_optimality_index].capacity -= amount;

        if (ctx.IsFiniteScalar(optimality_value)) {
          bool valid = (amount < 0) ? (optimality_value < middle_bound_value)
                                    : (optimality_value > middle_bound_value);
          if (valid) return amount;
        }
      } else {
        ctx.network.mutable_links()[max_optimality_index].capacity -= amount;
      }

      // Halve and retry
      amount /= CndOptimizationContext<T>::kBinarySearchDivisor;
    }
    return amount;
  }

  /**
   * @brief Budget-remaining iteration: expands the link with largest deviation from target phi = 1.0.
   *
   * When budget is not yet exhausted, adjusts a single link's capacity to move its optimality
   * function value toward the target of 1.0 (the KKT multiplier for the budget constraint
   * when budget is slack).
   */
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
    T best_value = base_value;
    for (int link_index = 0; link_index < ctx.network.number_of_links(); link_index++) {
      if (ProcessingCheck(ctx, link_index, middle_bound_value, od_cache)) {
        const T optimality_value = OptimalityFunctionFromCache(ctx, link_index, od_cache);
        if (!ctx.IsFiniteScalar(optimality_value)) continue;
        auto cur_delta = std::abs(middle_bound_value - optimality_value);
        if (cur_delta > max_delta) {
          max_optimality_index = link_index;
          max_delta = cur_delta;
          best_value = optimality_value;
        }
      }
    }

    T amount = ApproximateNonUpperBoundTransferAmount(
      ctx, max_optimality_index, best_value, middle_bound_value);
    amount = ValidatedNonUpperBoundTransferAmount(
      ctx, max_optimality_index, middle_bound_value, amount);
    ctx.network.mutable_links()[max_optimality_index].capacity += amount;
  }
};

}  // namespace TrafficAssignment

#endif  // OPTIMALITY_CONDITION_STEP_H

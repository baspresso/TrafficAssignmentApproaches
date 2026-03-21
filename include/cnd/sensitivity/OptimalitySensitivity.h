#ifndef OPTIMALITY_SENSITIVITY_H
#define OPTIMALITY_SENSITIVITY_H

#include <vector>
#include <cmath>
#include <unordered_set>
#include <Eigen/Dense>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "../CndOptimizationContext.h"
#include "../../tap/data/Link.h"

namespace TrafficAssignment {

namespace mp_sensitivity = boost::multiprecision;
using high_prec_sensitivity = long double;

/**
 * @brief Shared Jacobian-based sensitivity math for optimality condition computations.
 *
 * Extracts the reusable sensitivity analysis from OptimalityConditionStep, including:
 * - OD cache construction (Jacobian inverse, denominators, route-link sets)
 * - Optimality function phi_a computation (Chiou 2005, Eq. 9)
 * - Capacity derivative column construction
 *
 * Used by both OptimalityConditionStep (for iterative capacity updates) and
 * SensitivityEstimator (for analytical gradient computation).
 *
 * @tparam T Numeric type for flow/capacity computations.
 */
template <typename T>
class OptimalitySensitivity {
  using MatrixXd = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
  using MatrixHP = Eigen::Matrix<high_prec_sensitivity, Eigen::Dynamic, Eigen::Dynamic>;

public:
  /**
   * @brief Cached per-OD-pair data for efficient optimality function evaluation.
   *
   * Precomputes the route cost Jacobian inverse, the all-ones column eT*J^{-1}*e denominator,
   * and route-link membership sets. This avoids redundant matrix inversions when evaluating
   * phi_a across all design links.
   */
  struct OdCache {
    std::vector<std::vector<int>> routes;
    MatrixHP          jacobi_inverse;
    MatrixHP          e_col;
    MatrixHP          inv_transpose_e;
    high_prec_sensitivity denominator;
    T                 demand;
    std::vector<std::unordered_set<int>> route_link_sets;
    bool              valid = false;
  };

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
      c.e_col   = MatrixHP::Constant(rc, 1, high_prec_sensitivity(1));

      c.route_link_sets.resize(rc);
      for (int r = 0; r < rc; ++r) {
        c.route_link_sets[r].insert(c.routes[r].begin(), c.routes[r].end());
      }

      auto jacobi_hp = ConvertEigenMatrix(
          RoutesJacobiMatrix(c.routes, ctx.network.mutable_links()));
      c.jacobi_inverse = jacobi_hp.inverse();
      c.denominator    = (c.e_col.transpose() * c.jacobi_inverse * c.e_col)(0, 0);
      if (std::abs(c.denominator) <
          static_cast<high_prec_sensitivity>(CndOptimizationContext<T>::kDenominatorZeroGuard)) {
        continue;
      }
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
      const high_prec_sensitivity numerator =
          (c.inv_transpose_e.transpose() * cap_der)(0, 0);
      const T local_value = static_cast<T>(-numerator / c.denominator);
      if (!ctx.IsFiniteScalar(local_value)) return ctx.NaNValue();
      result += c.demand * local_value;
    }
    result /= ctx.budget_function_multiplier;
    return ctx.IsFiniteScalar(result) ? result : ctx.NaNValue();
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

  MatrixHP ConvertEigenMatrix(const MatrixXd& source) {
    MatrixHP target(source.rows(), source.cols());
    for (int i = 0; i < source.rows(); ++i) {
      for (int j = 0; j < source.cols(); ++j) {
        target(i, j) = static_cast<high_prec_sensitivity>(source(i, j));
      }
    }
    return target;
  }

  /// @brief Checks whether any OD pair has routes (needed for sensitivity-based gradient).
  bool HasRoutes(CndOptimizationContext<T>& ctx) {
    for (int od = 0; od < ctx.network.number_of_od_pairs(); ++od) {
      if (ctx.network.od_pairs()[od].GetRoutesCount() > 0) return true;
    }
    return false;
  }
};

}  // namespace TrafficAssignment

#endif  // OPTIMALITY_SENSITIVITY_H

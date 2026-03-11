#ifndef LINK_H
#define LINK_H

#include <vector>
#include <iostream>
#include <cmath>
#include <unordered_set>
#include <Eigen/Dense>

namespace TrafficAssignment {

  /**
   * Represents a directed link in a transportation network using the Bureau of Public Roads
   * (BPR) delay function: t_a(f) = t0 * (1 + b * (f/c)^p).
   *
   * The BPR function models congestion-dependent travel time on a link. As flow f approaches
   * or exceeds capacity c, travel time increases polynomially. This is the standard link
   * performance function used in traffic assignment. See Perederieieva et al. (2015) Section 2.
   *
   * @tparam T Data type for flow/capacity (e.g., double, long double).
   */
  template <class T>
  struct Link {
  static constexpr T MAX_ROUTE_DELAY = 1e9;

    // Network topology properties
    const int init;           ///< Source node ID (origin of the link).
    const int term;           ///< Target node ID (destination of the link).
    const int type;           ///< Link type/category (e.g., highway, arterial).

    // BPR function parameters
    T capacity;               ///< Link capacity c_a (vehicles/hour). Denominator in BPR ratio f/c.
    const T length;           ///< Physical length of the link (km/miles).
    const T free_flow_time;   ///< Travel time under free-flow conditions (minutes).
    const T b;                ///< BPR scaling factor (typical value: 0.15).
    const T power;            ///< BPR exponent (typical value: 4).
    const T speed;            ///< Speed limit (km/h or mph).
    const T toll;             ///< Toll cost for using the link (monetary units).
    const int power_int_;            ///< Integer cast of power for fast exponentiation via IntPow().
    const bool is_integer_power_;    ///< True if power is exactly integer; enables IntPow() optimization.

    T flow = 0; ///< Current traffic flow on the link (vehicles/hour).

    /**
     * @brief Constructs a Link object with given parameters.
     * @param init Source node ID.
     * @param term Target node ID.
     * @param capacity Maximum flow capacity.
     * @param length Physical length.
     * @param free_flow_time Free-flow travel time.
     * @param b BPR scaling factor.
     * @param power BPR exponent.
     * @param speed Speed limit.
     * @param toll Toll cost.
     * @param type Link type.
     */
    Link(int init, int term, T capacity, T length,
      T free_flow_time, T b, T power, T speed, T toll, int type) :
      init(init), term(term), type(type), capacity(capacity), length(length),
      free_flow_time(free_flow_time), b(b), power(power), speed(speed), toll(toll),
      power_int_(static_cast<int>(power)),
      is_integer_power_(power == static_cast<T>(static_cast<int>(power))), flow(0) { };

    /// @brief Default destructor (no dynamic memory to manage).
    ~Link() = default;

    /// @brief Fast integer exponentiation via binary exponentiation (squaring). O(log exp).
    static T IntPow(T base, int exp) {
      if (exp < 0) return T(1) / IntPow(base, -exp);
      T result = 1;
      while (exp > 0) {
        if (exp & 1) result *= base;
        base *= base;
        exp >>= 1;
      }
      return result;
    }

    /**
     * @brief Computes link travel time c_a(f_a) using the BPR delay function.
     *
     * Formula: t_a(f) = t0 * (1 + b * (f/c)^p)
     *
     * This is the link cost function used throughout traffic assignment.
     * See Perederieieva et al. (2015) Section 2.
     *
     * @param temp_flow Flow value to compute delay for. Defaults to current link flow.
     * @return Travel time on the link (minutes).
     */
    T Delay(T temp_flow = -1) const {
      if (temp_flow == -1) {
        temp_flow = flow;
      }
      if (capacity == 0 || b == 0) {
        return free_flow_time;
      }
      T ratio = temp_flow / capacity;
      if (ratio < 0) ratio = 0;
      if (is_integer_power_) {
        return free_flow_time + b * free_flow_time * IntPow(ratio, power_int_);
      }
      return free_flow_time + b * free_flow_time * std::pow(ratio, power);
    }

    /**
     * @brief Computes the Beckmann objective integrand: integral_0^f t_a(x) dx.
     *
     * Formula: t0 * (f + b * c * (f/c)^(p+1) / (p+1))
     *
     * This is the per-link contribution to the Beckmann objective function T(f).
     * The user equilibrium is the solution to min T(f) = sum_a integral_0^{f_a} t_a(x) dx.
     * See Bar-Gera (2010) Eq. 1, Perederieieva et al. (2015) Eq. 2.
     *
     * @param temp_flow Upper limit of integration. Defaults to current link flow.
     * @return Integral value (vehicle-minutes).
     */
    T DelayInteg(T temp_flow = -1) const {
      if (temp_flow == -1) {
        temp_flow = flow;
      }
      if (capacity == 0 || b == 0) {
        return free_flow_time * temp_flow;
      }
      T ratio = temp_flow / capacity;
      if (ratio < 0) ratio = 0;
      if (is_integer_power_) {
        return free_flow_time * (temp_flow + b * capacity * IntPow(ratio, power_int_ + 1) / (power + 1));
      }
      return free_flow_time * (temp_flow + b * capacity * std::pow(ratio, power + 1) / (power + 1));
    }

    /**
     * @brief Computes dc_a/df_a, the first derivative of the BPR delay function.
     *
     * Formula: t0 * b * p * (f/c)^(p-1) / c
     *
     * Used as the denominator in the Newton step for flow shifting:
     * delta_F = (C_max - C_min) / sum(dc_a/df_a). See Perederieieva et al. (2015) Eq. 14.
     *
     * @param temp_flow Flow value to evaluate at. Defaults to current link flow.
     * @return Rate of change of delay with respect to flow (minutes/vehicle).
     */
    T DelayDer(T temp_flow = -1) const {
      if (temp_flow == -1) {
        temp_flow = flow;
      }
      if (capacity == 0 || b == 0) {
        return 0;
      }
      T ratio = temp_flow / capacity;
      if (ratio < 0) ratio = 0;
      if (is_integer_power_) {
        return free_flow_time * b * power * IntPow(ratio, power_int_ - 1) / capacity;
      }
      return free_flow_time * b * power * std::pow(ratio, power - 1) / capacity;
    }

    /**
     * @brief Computes d^2 c_a / df_a^2, the second derivative of the BPR delay function.
     *
     * Formula: t0 * b * p * (p-1) * (f/c)^(p-2) / c^2
     *
     * Used in second-order methods (e.g., TapasAdvancedGradientDescentShiftMethod).
     *
     * @param temp_flow Flow value to evaluate at. Defaults to current link flow.
     * @return Second derivative of delay (minutes/vehicle^2).
     */
    T DelaySecondDer(T temp_flow = -1) const {
      if (temp_flow == -1) {
        temp_flow = flow;
      }
      if (capacity == 0) {
        return 0;
      }
      T ratio = temp_flow / capacity;
      if (ratio < 0) ratio = 0;
      if (is_integer_power_) {
        return free_flow_time * b * power * (power - 1) * IntPow(ratio, power_int_ - 2) / (capacity * capacity);
      }
      return free_flow_time * b * power * (power - 1) * std::pow(ratio, power - 2) / (capacity * capacity);
    }

    /**
     * @brief Computes total delay for a sequence of links.
     * @param links Vector of all links in the network.
     * @param links_list Indices of links to include in the calculation.
     * @return Sum of delays for the specified links (minutes).
     */
    static T GetLinksDelay(const std::vector <Link <T>>& links, const std::vector <int>& links_list) {
      T ans = 0;
      for (auto now : links_list) {
        ans += links[now].Delay();
      }
      return ans;
    }

    /**
     * @brief Computes total first derivative of delay for a sequence of links.
     * @param links Vector of all links in the network.
     * @param links_list Indices of links to include in the calculation.
     * @return Sum of delay derivatives for the specified links (minutes/vehicle).
     */
    static T GetLinksDelayDer(const std::vector <Link <T>>& links, const std::vector <int>& links_list) {
      T ans = 0;
      for (auto now : links_list) {
        ans += links[now].DelayDer();
      }
      return ans;
    }

    /**
     * @brief Computes total second derivative of delay for a sequence of links.
     * @param links Vector of all links in the network.
     * @param links_list Indices of links to include in the calculation.
     * @return Sum of second derivatives for the specified links (minutes/vehicle²).
     */
    static T GetLinksDelaySecondDer(const std::vector <Link <T>>& links, const std::vector <int>& links_list) {
      T ans = 0;
      for (auto now : links_list) {
        ans += links[now].DelaySecondDer();
      }
      return ans;
    }
    
    /**
     * @brief Computes dt_a/dc_a, the partial derivative of delay with respect to capacity.
     *
     * Formula: -p * b * t0 * (f/c)^p / c
     *
     * Used in bilevel CNDP optimization to evaluate sensitivity of travel time to
     * capacity changes (gradient of the lower-level objective w.r.t. design variables).
     *
     * @param temp_flow Flow value to evaluate at. Defaults to current link flow.
     * @return Partial derivative of delay w.r.t. capacity (negative: more capacity = less delay).
     */
    T DelayCapacityDer(T temp_flow = -1) const {
      if (temp_flow == -1) {
        temp_flow = flow;
      }
      if (capacity == 0 || b == 0) {
        return 0;
      }
      T ratio = temp_flow / capacity;
      if (ratio < 0) ratio = 0;
      if (is_integer_power_) {
        return -power * b * free_flow_time * IntPow(ratio, power_int_) / capacity;
      }
      return -power * b * free_flow_time * std::pow(ratio, power) / capacity;
    }

    /**
     * @brief Checks if any link in a list has non-zero capacity.
     * @param links Vector of all links in the network.
     * @param links_list Indices of links to check.
     * @return True if at least one link has capacity > 0, false otherwise.
     */
    static bool CheckNonZeroLinksCapacity(const std::vector <Link <T>>& links, const std::vector <int>& links_list) {
      bool ans = false;
      for (auto now : links_list) {
        if (links[now].capacity != 0) {
          ans = true;
        }
      }
      return ans;
    }
  };

  /**
   * @brief Computes element J(i,j) of the route cost Jacobian matrix.
   *
   * J(i,j) = sum of dc_a/df_a over links shared by routes i and j.
   * Used in the Krylatov (2023) shift method for Jacobian-based flow updates.
   */
  template <typename T>
  T CalculateJacobiElement(const std::vector<int>& route_i,
                           const std::vector<int>& route_j,
                           const std::vector <Link<T>>& links) {
    T sum = 0;
    std::unordered_set<int> links_i(route_i.begin(), route_i.end());
    
    for(int link_id : route_j) {
        if (links_i.count(link_id)) {
          sum += links[link_id].DelayDer();
        }
    }
    return sum;
  }

  /**
   * @brief Builds the full route cost Jacobian matrix J for a set of routes.
   *
   * J(i,j) = sum_{a in route_i intersect route_j} dc_a/df_a.
   * The matrix is symmetric since shared links contribute equally in both directions.
   * Used by the Krylatov (2023) method: delta_F = -J^{-1} * (C - e * C_target).
   */
  template <typename T>
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> RoutesJacobiMatrix(const std::vector<std::vector <int>>& routes, const std::vector <Link<T>>& links) {
    const size_t R = routes.size();
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> J(R, R);

    // Pre-build link sets for all routes to avoid repeated construction
    std::vector<std::unordered_set<int>> route_sets(R);
    for (size_t i = 0; i < R; ++i) {
      route_sets[i].insert(routes[i].begin(), routes[i].end());
    }

    // Exploit symmetry: J(i,j) = J(j,i) since both equal sum of DelayDer on shared links
    for (size_t i = 0; i < R; ++i) {
      // Diagonal
      T sum = 0;
      for (int link_id : routes[i]) {
        sum += links[link_id].DelayDer();
      }
      J(i, i) = sum;

      // Off-diagonal (upper triangle, mirror to lower)
      for (size_t j = i + 1; j < R; ++j) {
        T val = 0;
        for (int link_id : routes[j]) {
          if (route_sets[i].count(link_id)) {
            val += links[link_id].DelayDer();
          }
        }
        J(i, j) = val;
        J(j, i) = val;
      }
    }
    return J;
  }
}

#endif
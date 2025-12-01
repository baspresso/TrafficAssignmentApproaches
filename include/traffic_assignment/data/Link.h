#ifndef LINK_H
#define LINK_H

#include <vector>
#include <iostream>
#include <cmath>


namespace TrafficAssignment {

  /**
   * Represents a directed link in a transportation network using the BPR delay function.
   * @tparam T Data type for flow/capacity (e.g., double).
   */
  template <class T>
  struct Link {
  static constexpr T MAX_ROUTE_DELAY = 1e9;

    // Network topology properties
    const int init;           ///< Source node ID (origin of the link).
    const int term;           ///< Target node ID (destination of the link).
    const int type;           ///< Link type/category (e.g., highway, arterial).

    // BPR function parameters
    T capacity;        
    const T length;           ///< Physical length of the link (km/miles).
    const T free_flow_time;   ///< Travel time under free-flow conditions (minutes).
    const T b;                ///< BPR scaling factor (typical value: 0.15).
    const T power;            ///< BPR exponent (typical value: 4).
    const T speed;            ///< Speed limit (km/h or mph).
    const T toll;             ///< Toll cost for using the link (monetary units).

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
      free_flow_time(free_flow_time), b(b), power(power), speed(speed), toll(toll), flow(0) { };

    /// @brief Default destructor (no dynamic memory to manage).
    ~Link() = default;

    /**
     * @brief Computes travel time using the BPR delay function.
     * @param temp_flow Flow value to compute delay for. Defaults to current flow.
     * @return Travel time (delay).
     */
    T Delay(T temp_flow = -1) const {
      if (temp_flow == -1) {
        temp_flow = flow;
      }
      if (capacity == 0 || b == 0) {
        return free_flow_time;
      }
      return free_flow_time + b * free_flow_time * std::pow(temp_flow, power) / std::pow(capacity, power);
    }

    /**
     * @brief Computes the integral of the BPR delay function (used in objective functions).
     * @param temp_flow Upper limit of integration. Defaults to current flow.
     * @return Integral value (unit: vehicle-minutes).
     */
    T DelayInteg(T temp_flow = -1) const {
      if (temp_flow == -1) {
        temp_flow = flow;
      }
      if (capacity == 0 || b == 0) {
        return free_flow_time * temp_flow;
      }
      return free_flow_time * (temp_flow + b * capacity * std::pow(temp_flow / capacity, power + 1) / (power + 1));
    }

    /**
     * @brief Computes the first derivative of the BPR delay function (used in gradient calculations).
     * @param temp_flow Flow value to evaluate at. Defaults to current flow.
     * @return Rate of change of delay with respect to flow (minutes/vehicle).
     */
    T DelayDer(T temp_flow = -1) const {
      if (temp_flow == -1) {
        temp_flow = flow;
      }
      if (capacity == 0 || b == 0) {
        return 0;
      }
      return free_flow_time * b * power * std::pow(temp_flow / capacity, power - 1) / capacity;
    }

    /**
     * @brief Computes the second derivative of the BPR delay function (used in Hessian calculations).
     * @param temp_flow Flow value to evaluate at. Defaults to current flow.
     * @return Second derivative of delay (minutes/vehicle²).
     */
    T DelaySecondDer(T temp_flow = -1) const {
      if (temp_flow == -1) {
        temp_flow = flow;
      }
      if (capacity == 0) {
        return 0;
      }
      return free_flow_time * b * power * (power - 1) * std::pow(temp_flow, power - 2) / std::pow(capacity, power);
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
    static T GetLinksDelayDer(std::vector <Link <T>>& links, const std::vector <int>& links_list) {
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
    static T GetLinksDelaySecondDer(std::vector <Link <T>>& links, const std::vector <int>& links_list) {
      T ans = 0;
      for (auto now : links_list) {
        ans += links[now].DelaySecondDer();
      }
      return ans;
    }
    
    /**
     * @brief Checks if any link in a list has non-zero capacity.
     * @param links Vector of all links in the network.
     * @param links_list Indices of links to check.
     * @return True if at least one link has capacity > 0, false otherwise.
     */
    static bool CheckNonZeroLinksCapacity(std::vector <Link <T>>& links, const std::vector <int>& links_list) {
      bool ans = false;
      for (auto now : links_list) {
        if (links[now].capacity != 0) {
          ans = true;
        }
      }
      return ans;
    }
  };
}

#endif
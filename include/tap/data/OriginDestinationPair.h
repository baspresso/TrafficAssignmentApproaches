#ifndef ORIGIN_DESTINATION_PAIR_H
#define ORIGIN_DESTINATION_PAIR_H

#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include "Link.h"
#include <memory>

namespace TrafficAssignment {

template <typename T>
class Network;

template <class T>
class OriginDestinationPair;

} // namespace TrafficAssignment

// Include Network.h AFTER class declaration but before method implementations
#include "../core/Network.h"

namespace TrafficAssignment {

/**
* @brief Manages route assignments and flow distribution for a single origin-destination (O-D) pair.
* @tparam T Data type for flow/demand values (e.g., double).
*/
template <class T>
class OriginDestinationPair {
static constexpr T MAX_ROUTE_DELAY = 1e9; ///< Upper bound for route delays (prevents overflow).
public:
  /**
  * @brief Constructs an O-D pair with given parameters.
  * @param origin Origin node ID.
  * @param dest Destination node ID.
  * @param demand Total travel demand (vehicles/hour).
  * @param links Reference to the network's global link list (directly modified!).
  * @param adjacency_list Reference to the network's adjacency list for graph traversal.
  */
  OriginDestinationPair(int origin, int dest, T demand, 
                         Network<T>& network)
        : origin_(origin), dest_(dest), demand_(demand), 
          network_(network) { // Store shared_ptr directly
  }

  ~OriginDestinationPair() = default;

  // Basic getters
  // -------------

  /// @brief Returns total travel demand for this O-D pair.
  T GetDemand() const{
    return demand_;
  }

  /// @brief Returns (origin, destination) node IDs as a pair.
  std::pair <int, int> GetOriginDestination() const {
    return { origin_, dest_ };
  }

  /// @brief Returns the number of active routes.
  int GetRoutesCount() const {
    return routes_.size();
  }

  /// @brief Returns all routes (each represented as a list of link IDs).
  std::vector <std::vector <int>> GetRoutes() const {
    return routes_;
  }

  /// @brief Returns flow values for each route.
  std::vector <T> GetRoutesFlow() const {
    return routes_flow_;
  }

  /**
   * @brief Returns unique links used by all routes of this O-D pair.
   * @note O(n log n) due to sorting and deduplication.
   */
  std::vector <int> GetUsedLinks() {
    std::vector <int> used_links;
    for (int i = 0; i < routes_.size(); i++) {
      for (auto now : routes_[i]) {
        used_links.push_back(now);
      }
    }
    std::sort(used_links.begin(), used_links.end());
    used_links.resize(std::unique(used_links.begin(), used_links.end()) - used_links.begin());
    return used_links;
  }

  /**
   * @brief Finds the first route with all zero-capacity links and returns its delay.
   * @return Delay of the first zero-capacity route, or -1 if none exist.
   */
  T GetZeroCapacityRouteDelay() {
    for (int i = 0; i < this->routes_.size(); i++) {
      if (!Link<T>::CheckNonZeroLinksCapacity(this->network_.links(), this->routes_[i])) {
        return Link<T>::GetLinksDelay(this->network_.links(), this->routes_[i]);
      }
    }
    return -1;
  }

  std::vector<T> SetRoutesFlow(std::vector<T> routes_flow_update) {
      // 1. Remove previous contributions from network
      for (const auto& [link_id, flow] : links_flow_) {
          network_.SetLinkFlow(link_id, -flow); // Subtract old flow
      }
      links_flow_.clear();

      // 2. Filter invalid routes (non-positive flows)
      for (int i = routes_flow_update.size() - 1; i >= 0; --i) {
          if (routes_flow_update[i] <= 0) {
              routes_.erase(routes_.begin() + i);
              routes_flow_update.erase(routes_flow_update.begin() + i);
          }
      }

      if (routes_flow_update.empty()) {
          throw std::invalid_argument("All route flows are non-positive");
      }

      // 3. Normalize flows to meet total demand
      T total_flow = std::accumulate(routes_flow_update.begin(), 
                                    routes_flow_update.end(), T(0));
      routes_flow_.resize(routes_flow_update.size());
      
      for (size_t i = 0; i < routes_flow_update.size(); ++i) {
          routes_flow_[i] = (routes_flow_update[i] * demand_) / total_flow;
      }

      // 4. Update network with new flows
      for (size_t i = 0; i < routes_.size(); ++i) {
          const T flow = routes_flow_[i];
          for (int link_id : routes_[i]) {
              links_flow_[link_id] += flow; // Track local contribution
              network_.SetLinkFlow(link_id, flow); // Update network
           }
      }

      return routes_flow_;
  }  

  // Graph algorithms
  // ----------------

  bool AddNewRoute(std::vector<int> new_route) {
    // Check for duplicate routes
    for (const auto& existing_route : routes_) {
        if (new_route == existing_route) {
            return false;
        }
    }

    // Add new route with initial zero flow
    routes_.push_back(new_route);
    routes_flow_.push_back(0);

    // Handle first route initialization
    if (routes_.size() == 1) {
        // Assign full demand to the first route
        routes_flow_[0] = demand_;
        // Update network flows and track contributions
        for (int link_id : new_route) {
          network_.SetLinkFlow(link_id, demand_);
          links_flow_[link_id] = demand_;
        }
    }
    else {
        // Initialize new links in the route with zero contribution
        for (int link_id : new_route) {
          if (!links_flow_.count(link_id)) {
            links_flow_[link_id] = 0;
          }
        }
    }

    return true;
  }
  T RoutesDelta() {
    T delay_min = MAX_ROUTE_DELAY, delay_max = 0, delay_route = 0;
    if (routes_.size() > 0) {
      delay_route = Link<T>::GetLinksDelay(network_.links(), routes_[0]);
      delay_min = delay_route;
      delay_max = delay_route;
    }
    for (int i = 1; i < routes_.size(); i++) {
      delay_route = Link<T>::GetLinksDelay(network_.links(), routes_[i]);
      delay_min = std::min(delay_min, delay_route);
      delay_max = std::max(delay_max, delay_route);
    }
    return abs(delay_max - delay_min);
  }

  std::vector <T> RoutesDelays() {
    std::vector <T> routes_delays(routes_.size());
    for (int i = 0; i < routes_.size(); i++) {
      routes_delays[i] = Link<T>::GetLinksDelay(network_.links(), routes_[i]);
    }
    return routes_delays;
  }

  bool CheckNonZeroCapacityRoutes() {
    bool fl = true;
    for (int i = 0; i < this->routes_.size(); i++) {
      if (!Link<T>::CheckNonZeroLinksCapacity(network_.links(), this->routes_[i])) {
        fl = false;
      }
    }
    return fl;
  }

  void Reset() {
    routes_.clear();
    routes_flow_.clear();
    for (auto& [key, value] : links_flow_) {
      value = 0;
    }
  }
  private:
    // Core data members
    int origin_, dest_; ///< Node IDs for origin and destination
    T demand_;          ///< Total travel demand between O-D
    Network<T>& network_; 
    
    // Route-specific state
    std::vector<std::vector<int>> routes_; ///< Active routes (each is a list of link IDs)
    std::vector<T> routes_flow_;           ///< Flow assigned to each route
    std::unordered_map<int, T> links_flow_; ///< Flow on each link attributed to this O-D pair
  };
}

#endif // ORIGIN_DESTINATION_PAIR_H
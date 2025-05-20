#ifndef ORIGIN_DESTINATION_PAIR_H
#define ORIGIN_DESTINATION_PAIR_H

#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include "Link.h"

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
    OriginDestinationPair(int& origin, int& dest, T& demand, std::vector <Link <T>>& links, std::vector <std::vector <int>>& adjacency_list) :
      origin_(origin), dest_(dest), demand_(demand), links_(links), adjacency_list_(adjacency_list) { };

    ~OriginDestinationPair() = default;

    // Basic getters
    // -------------

    /// @brief Returns total travel demand for this O-D pair.
    T GetDemand() {
      return demand_;
    }

    /// @brief Returns (origin, destination) node IDs as a pair.
    std::pair <int, int> GetOriginDestination() {
      return { origin_, dest_ };
    }

    /// @brief Returns the number of active routes.
    int GetRoutesCount() {
      return routes_.size();
    }

    /// @brief Returns all routes (each represented as a list of link IDs).
    std::vector <std::vector <int>> GetRoutes() {
      return routes_;
    }

    /// @brief Returns flow values for each route.
    std::vector <T> GetRoutesFlow() {
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
        if (!Link<T>::CheckNonZeroLinksCapacity(this->links_, this->routes_[i])) {
          return Link<T>::GetLinksDelay(this->links_, this->routes_[i]);
        }
      }
      return -1;
    }

    /**
     * @brief Replaces all routes and their flows (resets link flows).
     * @param new_routes New list of routes (each as link IDs).
     * @param new_routes_flow Corresponding flows for each route.
     * @warning Overwrites global link flows (`links_`).
     */
    template <typename U>
    void SetNewRoutesInfo(const std::vector <std::vector <int>>& new_routes, const std::vector <U>& new_routes_flow) {
      ClearFlow();
      routes_ = new_routes;
      routes_flow_.resize(new_routes_flow.size());
      for (int i = 0; i < new_routes_flow.size(); i++) {
        routes_flow_[i] = new_routes_flow[i];
        for (auto now : new_routes[i]) {
          links_[now].flow += new_routes_flow[i];
          links_flow_[now] += new_routes_flow[i];
        }
      }
    }

    /**
     * @brief Redistributes flows across routes while maintaining demand constraints.
     * @param flow_update Proposed flow changes (may be modified during normalization).
     * @return Normalized flow values after removal of invalid routes.
     * @throws std::invalid_argument if flow normalization fails (e.g., all flows <= 0).
     */
     std::vector<T> SetRoutesFlow(std::vector<T> flow_update) {
        // Remove routes with non-positive flow
        for (int i = flow_update.size() - 1; i >= 0; --i) {
            if (flow_update[i] <= 0) {
                routes_.erase(routes_.begin() + i);
                routes_flow_.erase(routes_flow_.begin() + i);
                flow_update.erase(flow_update.begin() + i);
            }
        }

        if (flow_update.empty()) {
            throw std::invalid_argument("All route flows are non-positive after cleanup");
        }

        // Normalize flows to meet total demand
        T total_flow = std::accumulate(flow_update.begin(), flow_update.end(), T(0));
        for (size_t i = 0; i < flow_update.size(); ++i) {
            routes_flow_[i] = flow_update[i] * demand_ / total_flow;
            for (int link_id : routes_[i]) {
                links_[link_id].flow += routes_flow_[i] - (demand_ / routes_.size());
                links_flow_[link_id] += routes_flow_[i] - (demand_ / routes_.size());
            }
        }
        return flow_update;
    }

    // Graph algorithms
    // ----------------

    /**
     * @brief Computes the shortest path using Dijkstra's algorithm.
     * @return Ordered list of link IDs from origin to destination.
     * @note Assumes non-negative link costs (BPR delay functions).
     */
    std::vector<int> BestRoute() {
      std::priority_queue<std::pair<T, int>, 
                          std::vector<std::pair<T, int>>, 
                          std::greater<std::pair<T, int>>> pq;
      std::unordered_map<int, T> dist;
      std::unordered_map<int, int> prev_link;

      pq.push({0, origin_});
      dist[origin_] = 0;

      while (!pq.empty()) {
        auto [current_dist, u] = pq.top();
        pq.pop();

        if (u == dest_) break;
        if (current_dist > dist[u]) continue;

        for (int link_id : adjacency_list_[u]) {
          const Link<T>& link = links_[link_id];
          T new_dist = current_dist + link.Delay();
          if (!dist.count(link.term) || new_dist < dist[link.term]) {
            dist[link.term] = new_dist;
            prev_link[link.term] = link_id;
            pq.push({new_dist, link.term});
          }
        }
      }

      // Reconstruct path
      std::vector<int> path;
      int current = dest_;
      while (prev_link.count(current)) {
        path.push_back(prev_link[current]);
        current = links_[prev_link[current]].init;
      }
      std::reverse(path.begin(), path.end());
      return path;
    }

    void ClearFlow() {
      std::unordered_set <int> used_links;
      for (const auto& route : routes_) {
        for (const auto& link : route) {
          used_links.insert(link);
        }
      }
      for (const auto& link : used_links) {
        links_[link].flow -= links_flow_[link];
        links_flow_[link] = 0;
      }
    }

    T BestRouteDelay() {
      return Link<T>::GetLinksDelay(links_, BestRoute());
    }

    bool AddNewRoute(std::vector <int> new_route) {
      for (auto route : routes_) {
        if (new_route == route) {
          return false;
        }
      }
      routes_.push_back(new_route);
      routes_flow_.push_back(0);
      if (routes_.size() == 1) {
        routes_flow_[0] += demand_;
      }
      for (auto now : routes_[routes_.size() - 1]) {
        if (routes_.size() == 1) {
          links_[now].flow += demand_;
          links_flow_[now] = demand_;
        }
        else if (!links_flow_.count(now)) {
          links_flow_[now] = 0;
        }
      }
      return true;
    }

    T RoutesDelta() {
      T delay_min = MAX_ROUTE_DELAY, delay_max = 0, delay_route = 0;
      if (routes_.size() > 0) {
        delay_route = Link<T>::GetLinksDelay(links_, routes_[0]);
        delay_min = delay_route;
        delay_max = delay_route;
      }
      for (int i = 1; i < routes_.size(); i++) {
        delay_route = Link<T>::GetLinksDelay(links_, routes_[i]);
        delay_min = std::min(delay_min, delay_route);
        delay_max = std::max(delay_max, delay_route);
      }
      return abs(delay_max - delay_min);
    }

    std::vector <T> RoutesDelays() {
      std::vector <T> routes_delays(routes_.size());
      for (int i = 0; i < routes_.size(); i++) {
        routes_delays[i] = Link<T>::GetLinksDelay(links_, routes_[i]);
      }
      return routes_delays;
    }

    void SetDefaultFlow() {
      ClearFlow();
      for (int i = 0; i < routes_.size(); i++) {
        routes_flow_[i] = demand_ / routes_.size();
        for (auto now : routes_[i]) {
          links_[now].flow += demand_ / routes_.size();
          links_flow_[now] += demand_ / routes_.size();
        }
      }
    }

    bool CheckNonZeroCapacityRoutes() {
      bool fl = true;
      for (int i = 0; i < this->routes_.size(); i++) {
        if (!Link<T>::CheckNonZeroLinksCapacity(this->links_, this->routes_[i])) {
          fl = false;
        }
      }
      return fl;
    }

  private:
    // Core data members
    int origin_, dest_; ///< Node IDs for origin and destination
    T demand_;          ///< Total travel demand between O-D
    std::vector<Link<T>>& links_; ///< Reference to network's global links (CAUTION: modified directly!)
    const std::vector<std::vector<int>>& adjacency_list_; ///< Reference to network's adjacency list

    // Route-specific state
    std::vector<std::vector<int>> routes_; ///< Active routes (each is a list of link IDs)
    std::vector<T> routes_flow_;           ///< Flow assigned to each route
    std::unordered_map<int, T> links_flow_; ///< Flow on each link attributed to this O-D pair
  };
}

#endif // ORIGIN_DESTINATION_PAIR_H
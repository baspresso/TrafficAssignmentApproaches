// core/Network.h
#ifndef NETWORK_H
#define NETWORK_H

#include "../data/Link.h"
#include "../data/OriginDestinationPair.h"
#include <vector>
#include <unordered_map>
#include <map>

namespace TrafficAssignment {

/**
   * Represents a complete transportation network with nodes, links, and origin-destination pairs.
   * Provides infrastructure for traffic assignment algorithms and network analysis.
   * @tparam T Data type for flow/capacity values (e.g., double).
   */
template <typename T>
class Network {
public:

    /**
     * @brief Constructs a Network object with given parameters.
     * @param nodes Total number of nodes in the network.
     * @param zones Total number of traffic zones.
     * @param links Vector of links comprising the network.
     * @param trips Matrix of travel demands between zones.
     * @param adjacency Forward adjacency list representation.
     * @param reverse_adjacency Reverse adjacency list representation.
     */
    Network(int nodes, int zones,
            std::vector<Link<T>> links,
            std::vector<std::vector<T>> trips,
            std::vector<std::vector<int>> adjacency,
            std::vector<std::vector<int>> reverse_adjacency)
        : number_of_nodes_(nodes),
          number_of_zones_(zones),
          links_(std::move(links)),
          adjacency_list_(std::move(adjacency)),
          reverse_adjacency_list_(std::move(reverse_adjacency)) 
    {
        InitializeODPairs(trips);
    }

    /// @brief Default destructor (uses default container destruction).
    ~Network() = default;

    // Network metadata accessors
    // ---------------------------
    
    /**
     * @brief Returns total number of nodes in the network.
     * @return Number of nodes as integer.
     */
    int number_of_links() const noexcept { return links_.size(); }

    /**
     * @brief Returns current number of links in the network.
     * @return Number of links as integer.
     */
    int number_of_od_pairs() const noexcept { return origin_destination_pairs_.size(); }

    // Core data accessors
    // -------------------
    
    /**
     * @brief Provides read-only access to network links.
     * @return Const reference to vector of links.
     */
    const std::vector<Link<T>>& links() const { return links_; }

    /**
     * @brief Provides read-only access to OD pairs.
     * @return Const reference to vector of origin-destination pairs.
     */
    const std::vector<OriginDestinationPair<T>>& od_pairs() const { return origin_destination_pairs_; }

    /**
     * @brief Provides read-only access to forward adjacency list.
     * @return Const reference to adjacency list structure.
     */
    const std::vector<std::vector<int>>& adjacency() const { return adjacency_list_; }

    /**
     * @brief Provides read-only access to reverse adjacency list.
     * @return Const reference to reverse adjacency structure.
     */
    const std::vector<std::vector<int>>& reverse_adjacency() const { return reverse_adjacency_list_; }

    // Mutable accessors for algorithm manipulation
    // --------------------------------------------
    
    /**
     * @brief Provides mutable access to network links.
     * @warning Direct modification affects all dependent components.
     * @return Reference to mutable vector of links.
     */
    std::vector<Link<T>>& mutable_links() { return links_; }

    /**
     * @brief Provides mutable access to OD pairs.
     * @warning Changes affect traffic assignment results.
     * @return Reference to mutable vector of OD pairs.
     */
    std::vector<OriginDestinationPair<T>>& mutable_od_pairs() { return origin_destination_pairs_; }

    // Network analysis operations
    // ---------------------------
    
    /**
     * @brief Computes shortest paths from a single origin to all destinations.
     * @param origin Node ID of the starting zone.
     * @return Vector of pairs containing OD pair index and corresponding shortest path.
     * @note Uses Dijkstra's algorithm with priority queue implementation.
     */
    std::vector<std::pair<int, std::vector<int>>> SingleOriginBestRoutes(int origin) const {
        std::priority_queue <std::pair <T, int>, std::vector <std::pair <T, int>>, std::greater <std::pair <T, int>>> q;
        std::vector <std::pair <int, std::vector <int>>> best_routes;
        int u;
        q.push({ 0, -1 });
        std::unordered_set <int> processed;
        std::unordered_map <int, int> used_link;
        T cur_delay;
        int cnt_not_processed = origin_info_[origin].size();
        std::vector <int> destinations;
        while (cnt_not_processed > 0) {
          if (q.empty()) {
            std::cout << "???\n";
            break;
          }
          if (q.top().second != -1) {
            u = links_[q.top().second].term;
          }
          else {
            u = origin;
          }
          if (processed.count(u)) {
            q.pop();
            continue;
          }
          cur_delay = q.top().first;
          processed.insert(u);
          used_link[u] = q.top().second;
          q.pop();
          for (auto now : adjacency_list_[u]) {
            if (!processed.count(links_[now].term)) {
              q.push({ cur_delay + links_[now].Delay(), now });
            }
          }
          if (origin_info_[origin].count(u)) {
            cnt_not_processed--;
            destinations.push_back(u);
          }
        }
        for (auto dest : destinations) {
          best_routes.push_back({ origin_info_[origin].at(dest), RestoreRoute(origin, dest, used_link) });
        }
        return best_routes;
    }

private:
    // Network topology properties
    const int number_of_nodes_;      ///< Total number of nodes in the network.
    const int number_of_zones_;      ///< Total number of traffic zones (origins/destinations).

    // Graph structure components
    std::vector<Link<T>> links_;                ///< All links in the transportation network.
    std::vector<std::vector<int>> adjacency_list_;       ///< Adjacency list for forward graph traversal.
    std::vector<std::vector<int>> reverse_adjacency_list_; ///< Reverse adjacency list for backward traversal.

    // Travel demand information
    std::vector<OriginDestinationPair<T>> origin_destination_pairs_; ///< All active origin-destination pairs.
    std::vector<std::map<int, int>> origin_info_; ///< Mapping from origin-destination indices to pair IDs.

    /**
     * @brief Initializes origin-destination pairs from trip matrix.
     * @param trips Matrix of travel demands between zones.
     * @note Automatically filters zero-demand OD pairs.
     */
    void InitializeODPairs(const std::vector<std::vector<T>>& trips) {
      origin_destination_pairs_.clear();
      origin_destination_pairs_.reserve(number_of_zones_ * number_of_zones_);

      for(int origin = 0; origin < number_of_zones_; ++origin) {
          for(int dest = 0; dest < number_of_zones_; ++dest) {
              if(trips[origin][dest] > 0) {
                  origin_destination_pairs_.emplace_back(
                      origin, 
                      dest, 
                      trips[origin][dest],
                      links_,        // Now passed as const reference
                      adjacency_list_
                  );
              }
          }
      }
    }

    /**
     * @brief Reconstructs route from destination back to origin.
     * @param origin Starting node ID.
     * @param dest Target node ID.
     * @param used_link Map of node-to-link connections from pathfinding.
     * @return Ordered list of link IDs from origin to destination.
     */
    std::vector <int> RestoreRoute(int origin, int dest, const std::unordered_map <int, int>& used_link) {
        int now = dest;
        std::vector <int> new_route;
        while (now != origin) {
          new_route.push_back(used_link.at(now));
          now = links_[used_link.at(now)].init;
        }
        reverse(new_route.begin(), new_route.end());
        return new_route;
      }
};

} // namespace TrafficAssignment

#endif // NETWORK_H
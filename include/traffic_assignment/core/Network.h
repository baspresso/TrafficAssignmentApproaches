// core/Network.h
#ifndef NETWORK_H
#define NETWORK_H

#include "../data/Link.h"
#include "../data/OriginDestinationPair.h"
#include <vector>
#include <unordered_map>
#include <map>

namespace TrafficAssignment {

template <typename T>
class Network {
public:
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

    // Metadata access
    int number_of_nodes() const noexcept { return number_of_nodes_; }
    int number_of_zones() const noexcept { return number_of_zones_; }
    int number_of_links() const noexcept { return links_.size(); }
    int number_of_od_pairs() const noexcept { return origin_destination_pairs_.size(); }

    // Core accessors
    const std::vector<Link<T>>& links() const { return links_; }
    const std::vector<OriginDestinationPair<T>>& od_pairs() const { return origin_destination_pairs_; }
    const std::vector<std::vector<int>>& adjacency() const { return adjacency_list_; }
    const std::vector<std::vector<int>>& reverse_adjacency() const { return reverse_adjacency_list_; }

    // Mutable accessors for algorithms
    std::vector<Link<T>>& mutable_links() { return links_; }
    std::vector<OriginDestinationPair<T>>& mutable_od_pairs() { return origin_destination_pairs_; }

    // Network analysis
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
    int number_of_nodes_;
    int number_of_zones_;
    std::vector<Link<T>> links_;
    std::vector<std::vector<int>> adjacency_list_;
    std::vector<std::vector<int>> reverse_adjacency_list_;
    std::vector<OriginDestinationPair<T>> origin_destination_pairs_;
    std::vector<std::map<int, int>> origin_info_;

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

#endif
#ifndef TRAFFIC_ASSIGNMENT_APPROACH_H
#define TRAFFIC_ASSIGNMENT_APPROACH_H

#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include "../data/Link.h"
#include "../data/OriginDestinationPair.h"
#include "../utils/DataProcessor.h"
#include "../utils/StatisticsRecorder.h"

namespace TrafficAssignment {

  template <typename T>
  class TrafficAssignmentApproach {
  public:
    TrafficAssignmentApproach(std::string dataset_name, T alpha = 1e-6) :
      number_of_nodes_(0), number_of_zones_(0), number_of_links_(0), number_of_origin_destination_pairs_(0), alpha_(alpha), dataset_name_(dataset_name) {
      data_processor_.LoadData(dataset_name_);

      number_of_nodes_ = data_processor_.GetNumberOfNodes();
      number_of_zones_ = data_processor_.GetNumberOfZones();
      number_of_links_ = data_processor_.GetNumberOfLinks();

      std::vector <Link <long double>> links = data_processor_.GetLinks();
      for (int i = 0; i < number_of_links_; i++) {
        links_.emplace_back(links[i]);
      }

      adjacency_list_.resize(number_of_nodes_);
      for (int i = 0; i < number_of_links_; i++) {
        adjacency_list_[links_[i].init].push_back(i);
      }
      reverse_adjacency_list_.resize(number_of_nodes_);
      for (int i = 0; i < number_of_links_; i++) {
        reverse_adjacency_list_[links_[i].term].push_back(i);
      }

      std::vector <std::vector <T>> trips = data_processor_.GetTrips();
      origin_info_.resize(number_of_zones_);
      for (int origin = 0; origin < number_of_zones_; origin++) {
        for (int dest = 0; dest < number_of_zones_; dest++) {
          if (trips[origin][dest] > 0) {
            origin_destination_pairs_.emplace_back(origin, dest, trips[origin][dest], links_, adjacency_list_);
            number_of_origin_destination_pairs_++;
            origin_info_[origin][dest] = number_of_origin_destination_pairs_ - 1;
          }
        }
      }
    }

    virtual ~TrafficAssignmentApproach() {}

    std::vector <Link <T>> GetLinks() {
      return links_;
    }

    std::vector <std::vector <int>> GetaAdjacencyList() {
      return adjacency_list_;
    }

    std::vector <OriginDestinationPair <T>> GetOriginDestinationPairs() {
      return origin_destination_pairs_;
    }

    std::vector <std::vector <std::vector <int>>> GetOriginDestinationRoutes() {
      std::vector <std::vector <std::vector <int>>> result(this->number_of_origin_destination_pairs_);
      for (int i = 0; i < this->number_of_origin_destination_pairs_; i++) {
        result[i] = this->origin_destination_pairs_[i].GetRoutes();
      }
      return result;
    }

    std::vector <std::vector <T>> GetOriginDestinationRoutesFlows() {
      std::vector <std::vector <T>> result(this->number_of_origin_destination_pairs_);
      for (int i = 0; i < this->number_of_origin_destination_pairs_; i++) {
        result[i] = this->origin_destination_pairs_[i].GetRoutesFlow();
      }
      return result;
    }

    virtual std::string GetApproachName() = 0;

    virtual void ComputeTrafficFlows() = 0;

    T ObjectiveFunction() {
      T ans = 0;
      for (int i = 0; i < number_of_links_; i++) {
        ans += links_[i].DelayInteg();
      }
      return ans;
    }

    T RelativeGap() {
      T result = 1, numerator = 0, denominator = 0;
      for (int origin_index = 0; origin_index < number_of_zones_; origin_index++) {
        auto best_routes = SingleOriginBestRoutes(origin_index);
        for (auto now : best_routes) {
          numerator += origin_destination_pairs_[now.first].GetDemand() * Link<T>::GetLinksDelay(links_, now.second);
        }
      }
      for (auto& link : links_) {
        denominator += link.flow * link.Delay();
      }
      result -= numerator / denominator;
      return result;
    }

    T Delta() {
      T delta = 0;
      for (int t = 0; t < number_of_origin_destination_pairs_; t++) {
        delta = std::max(delta, origin_destination_pairs_[t].RoutesDelta());
      }
      return delta;
    }

    //template <typename U>
    //void SetPretunedSolution(TrafficAssignmentApproach <U>& solution) {
    //  this->start_ = solution.GetStart();
    //  this->time_on_statistics_ = solution.GetTimeOnStatistics();
    //  SetOriginDestinationPairsInfo(solution.GetOriginDestinationRoutes(), solution.GetOriginDestinationRoutesFlows());
    //}

    protected:
      T alpha_;
      std::vector <Link <T>> links_;
      // Stores outgoing links indeces for every node
      std::vector <std::vector <int>> adjacency_list_;
      // Stores incoming links indeces for every node
      std::vector <std::vector <int>> reverse_adjacency_list_;
      std::vector <OriginDestinationPair <T>> origin_destination_pairs_;
      // stores information about all destinations for the single origin
      std::vector <std::map <int, int>> origin_info_;
      int number_of_nodes_, number_of_zones_, number_of_links_, number_of_origin_destination_pairs_;
      std::string dataset_name_;
      DataProcessor data_processor_;
      StatisticsRecorder <T> statistics_recorder_;

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

      // Dijkstra's algorithm for finding shortest paths to destinatins for the single origin
      // every pair contains origin-destination pair index and the least cost route itself
      std::vector <std::pair <int, std::vector <int>>> SingleOriginBestRoutes(int origin) {
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
          best_routes.push_back({ origin_info_[origin][dest], RestoreRoute(origin, dest, used_link) });
        }
        return best_routes;
      }

      void AddNewOriginDestinationRoutes(const std::vector <std::pair <int, std::vector <int>>>& routes) {
        for (auto now : routes) {
          origin_destination_pairs_[now.first].AddNewRoute(now.second);
        }
      }

      template <typename U>
      void SetOriginDestinationPairsInfo(const std::vector <std::vector <std::vector <int>>>& od_routes, const std::vector <std::vector <U>>& od_routes_flows) {
        for (int i = 0; i < this->number_of_origin_destination_pairs_; i++) {
          this->origin_destination_pairs_[i].SetNewRoutesInfo(od_routes[i], od_routes_flows[i]);
        }
      }

      T RoutesLinksObjectiveFunction(const std::vector <std::vector <int>>& routes) {
        std::unordered_set <int> links;
        for (const auto& route : routes) {
          for (const auto& link_index : route) {
            links.insert(link_index);
          }
        }
        T objective_function_on_links = 0;
        for (auto link_index : links) {
          objective_function_on_links += links_[link_index].DelayInteg();
        }
        return objective_function_on_links;
      }
  };
}

#endif
#ifndef ORIGIN_DESTINATION_PAIR_H
#define ORIGIN_DESTINATION_PAIR_H

#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include "Link.h"

namespace TrafficAssignment {

  template <class T>
  class OriginDestinationPair {
#define MAX_ROUTE_DELAY 1e9
  public:
    OriginDestinationPair(int& origin, int& dest, T& demand, std::vector <Link <T>>& links, std::vector <std::vector <int>>& adjacency_list) :
      origin_(origin), dest_(dest), demand_(demand), links_(links), adjacency_list_(adjacency_list), current_delay_(0) { };

    ~OriginDestinationPair() { }

    T GetDemand() {
      return demand_;
    }

    std::pair <int, int> GetOriginDestination() {
      return { origin_, dest_ };
    }

    int GetRoutesCount() {
      return routes_.size();
    }

    std::vector <std::vector <int>> GetRoutes() {
      return routes_;
    }

    std::vector <T> GetRoutesFlow() {
      return routes_flow_;
    }

    std::vector <int> GetUsedLinks() {
      std::vector <int> used_links;
      for (int i = 0; i < routes_.size(); i++) {
        for (auto now : routes_[i]) {
          used_links.push_back(now);
        }
      }
      sort(used_links.begin(), used_links.end());
      used_links.resize(unique(used_links.begin(), used_links.end()) - used_links.begin());
      return used_links;
    }

    T GetZeroCapacityRouteDelay() {
      for (int i = 0; i < this->routes_.size(); i++) {
        if (!Link<T>::CheckNonZeroLinksCapacity(this->links_, this->routes_[i])) {
          return Link<T>::GetLinksDelay(this->links_, this->routes_[i]);
        }
      }
      return -1;
    }

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
    // Makes redistribution of routes flows, deletes routes that's are going to get a negative flow in flow_update
    // sets default flow for routes that are not in routes_indeces
    // sets given flow update for routes that are in routes indeces and perfomes normalisation of this flow
    // in order to satisfy demand requirement
    // this function assumes that routes_indeces are given in a sorted order
    // returns routes indeces
    std::vector <T> SetRoutesFlow(std::vector <T> flow_update) {
      int non_positive_cnt = 0;
      for (int i = 0; i < flow_update.size(); i++) {
        if (flow_update[i] <= 0) {
          non_positive_cnt++;
        }
      }
      if (non_positive_cnt > 0) {
        ClearFlow();
        for (int i = flow_update.size() - 1; i >= 0; i--) {
          if (flow_update[i] <= 0) {
            routes_.erase(routes_.begin() + i);
            routes_flow_.pop_back();
            flow_update.erase(flow_update.begin() + i);
          }
        }
      }
      SetDefaultFlow();
      T all_flow_update = 0;
      for (int i = 0; i < flow_update.size(); i++) {
        all_flow_update += flow_update[i];
      }
      for (int i = 0; i < flow_update.size(); i++) {
        flow_update[i] *= demand_ / all_flow_update;
        routes_flow_[i] = flow_update[i];
        for (auto now : routes_[i]) {
          links_[now].flow += flow_update[i] - demand_ / routes_.size();
          links_flow_[now] += flow_update[i] - demand_ / routes_.size();
        }
      }
      return flow_update;
    }

    std::vector <int> BestRoute() {
      std::priority_queue <std::pair <T, int>, std::vector <std::pair <T, int>>, std::greater <std::pair <T, int>>> q;
      int u;
      q.push({ 0, -1 });
      T cur_delay;
      std::unordered_set <int> processed;
      std::unordered_map <int, int> used_link;
      T temp;
      while (!q.empty()) {
        if (q.top().second != -1)
          u = links_[q.top().second].term;
        else
          u = origin_;
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
        if (u == dest_) {
          break;
        }
      }
      int now = dest_;
      std::vector <int> new_route;
      while (now != origin_) {
        new_route.push_back(used_link[now]);
        now = links_[used_link[now]].init;
      }
      reverse(new_route.begin(), new_route.end());
      return new_route;
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
        delay_min = min(delay_min, delay_route);
        delay_max = max(delay_max, delay_route);
      }
      //if (abs(delay_max - delay_min) > 1)
      //cout << abs(delay_max - delay_min) << '\n';
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

    void ShowRoutesFlow() {
      for (int i = 0; i < this->routes_.size(); i++) {
        std::cout << this->origin_ << ' ' << this->dest_ << ' ';
        std::cout << links_[routes_[i][0]].init << ',';
        for (auto now : routes_[i]) {
          std::cout << links_[now].term << ',';
        }
        std::cout << ' ' << this->routes_flow_[i] << '\n';
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
    int origin_, dest_;
    T demand_, current_delay_;
    std::vector <Link <T>>& links_;
    // Stores outgoing links indeces for every node
    const std::vector <std::vector <int>>& adjacency_list_;
    std::vector <std::vector <int>> routes_;
    std::vector <T> routes_flow_;
    std::unordered_map <int, T> links_flow_;
  };
}

#endif // ORIGIN_DESTINATION_PAIR_H
#ifndef ROUTE_BASED_KRYLATOV_2023_SHIFT_METHOD_H
#define ROUTE_BASED_KRYLATOV_2023_SHIFT_METHOD_H

#include <unordered_set>
#include "OriginDestinationPair.h"
#include "RouteBasedShiftMethod.h"
#include "Eigen\Dense"


namespace TrafficAssignment {
  template <typename T>
  class RouteBasedKrylatov2023ShiftMethod : public RouteBasedShiftMethod <T> {
    using MatrixXd = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
  public:
    RouteBasedKrylatov2023ShiftMethod(std::vector <Link <T>>& links, std::vector <OriginDestinationPair <T>>& origin_destination_pairs) :
      RouteBasedShiftMethod<T>(links, origin_destination_pairs) {}

    std::vector <T> FlowShift(int od_pair_index) {
      int routes_count = this->origin_destination_pairs_[od_pair_index].GetRoutesCount();
      std::vector <T> routes_flow = this->origin_destination_pairs_[od_pair_index].GetRoutesFlow();
      MatrixXd flow_delta = (-1) * RoutesJacobiMatrix(od_pair_index).inverse() * (RoutesDelayColumn(od_pair_index) - EColumn(od_pair_index) * TargetDelay(od_pair_index));
      std::vector <T> flow_shift(routes_count);

      for (int i = 0; i < routes_count; i++) {
        flow_shift[i] = routes_flow[i] + flow_delta(i, 0);
      }
      return flow_shift;
    }

    std::string GetApproachName() override {
      return "RouteBasedKrylatov2023";
    }

  private:
    MatrixXd RoutesJacobiMatrix(int od_pair_index) {
      std::vector <std::vector <int>> routes = this->origin_destination_pairs_[od_pair_index].GetRoutes();
      int routes_count = this->origin_destination_pairs_[od_pair_index].GetRoutesCount();
      MatrixXd jacobi_matrix(routes_count, routes_count);
      for (int i = 0; i < routes_count; i++) {
        std::unordered_set <int> links_cur_route;
        for (auto now : routes[i]) {
          links_cur_route.insert(now);
        }
        for (int j = 0; j < routes_count; j++) {
          jacobi_matrix(i, j) = 0;
          for (auto now : routes[j]) {
            if (links_cur_route.count(now)) {
              jacobi_matrix(i, j) += this->links_[now].DelayDer();
            }
          }
        }
      }
      return jacobi_matrix;
    }

    MatrixXd RoutesDelayColumn(int od_pair_index) {
      std::vector <std::vector <int>> routes = this->origin_destination_pairs_[od_pair_index].GetRoutes();
      int routes_count = this->origin_destination_pairs_[od_pair_index].GetRoutesCount();
      MatrixXd delay_column(routes_count, 1);
      for (int i = 0; i < routes_count; i++) {
        delay_column(i, 0) = Link<T>::GetLinksDelay(this->links_, routes[i]);
      }
      return delay_column;
    }

    MatrixXd EColumn(int od_pair_index) {
      int routes_count = this->origin_destination_pairs_[od_pair_index].GetRoutesCount();
      MatrixXd e(routes_count, 1);
      for (int i = 0; i < routes_count; i++) {
        e(i, 0) = 1;
      }
      return e;
    }    

    T CheckConstRouteDelay(int od_pair_index) {
      std::vector <std::vector <int>> routes = this->origin_destination_pairs_[od_pair_index].GetRoutes();
      for (const auto& route : routes) {
        bool current_route_result = true;
        for (const auto& now : route) {
          if (this->links_[now].power != 0) {
            current_route_result = false;
          }
        }
        if (current_route_result) {
          return true;
        }
      }
      return false;
    }

    T ConstRouteDelay(int od_pair_index) {
      std::vector <std::vector <int>> routes = this->origin_destination_pairs_[od_pair_index].GetRoutes();
      T result = 0;
      for (const auto& route : routes) {
        bool current_route_result = true;
        for (const auto& now : route) {
          if (this->links_[now].power != 0) {
            current_route_result = false;
          }
        }
        if (current_route_result) {
          result = Link<T>::GetLinksDelay(this->links_, route);;
        }
      }
      return result;
    }

    T TargetDelay(int od_pair_index) {
      T delay_const_route = ConstRouteDelay(od_pair_index);
      if (CheckConstRouteDelay(od_pair_index)) {
        return ConstRouteDelay(od_pair_index);
      }
      else {
        MatrixXd e = EColumn(od_pair_index);
        MatrixXd jacobi_inv = RoutesJacobiMatrix(od_pair_index).inverse();
        return (e.transpose() * jacobi_inv * RoutesDelayColumn(od_pair_index))(0, 0) / (e.transpose() * jacobi_inv * e)(0, 0);
      }
    }

  };
}

#endif ROUTE_BASED_KRYLATOV_2023_SHIFT_METHOD_H
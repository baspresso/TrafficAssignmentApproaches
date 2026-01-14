#ifndef ROUTE_BASED_KRYLATOV_2023_SHIFT_METHOD_H
#define ROUTE_BASED_KRYLATOV_2023_SHIFT_METHOD_H

#include <unordered_set>
#include "./RouteBasedShiftMethod.h"
#include "C:/tools/vcpkg/packages/eigen3_x64-windows/include/eigen3/Eigen/Dense"


namespace TrafficAssignment {
  template <typename T>
  class RouteBasedKrylatov2023ShiftMethod : public RouteBasedShiftMethod <T> {
    using MatrixXd = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    
  public:
    explicit RouteBasedKrylatov2023ShiftMethod(Network<T>& network)
        : RouteBasedShiftMethod<T>(network) {}

    std::vector <T> FlowShift(int od_index) {
      auto& od_pair = this->od_pairs()[od_index];
      const int routes_count = od_pair.GetRoutesCount();
        
      MatrixXd flow_delta = CalculateFlowDelta(od_index);
      return ApplyFlowUpdate(od_pair, flow_delta);
    }

  private:

    MatrixXd CalculateFlowDelta(int od_index) {
        MatrixXd J = RoutesJacobiMatrix(od_index);
        MatrixXd delays = RoutesDelayColumn(od_index);
        MatrixXd e = EColumn(od_index);
        
        return (-1) * J.inverse() * (delays - e * TargetDelay(od_index));
    }
    
    std::vector<T> ApplyFlowUpdate(OriginDestinationPair<T>& od_pair, 
                                    const MatrixXd& delta) {
        std::vector<T> new_flows;
        const auto& current_flows = od_pair.GetRoutesFlow();
        
        for(int i = 0; i < delta.rows(); ++i) {
            new_flows.push_back(current_flows[i] + delta(i, 0));
        }
        
        return new_flows;
    }

    MatrixXd RoutesJacobiMatrix(int od_index) {

      const auto& routes = this->GetODPair(od_index).GetRoutes();
      MatrixXd J(routes.size(), routes.size());
        
      for(size_t i = 0; i < routes.size(); ++i) {
        for(size_t j = 0; j < routes.size(); ++j) {
          J(i,j) = CalculateJacobiElement(routes[i], routes[j]);
        }
      }
      return J;
    }

    T CalculateJacobiElement(const std::vector<int>& route_i,
                              const std::vector<int>& route_j) {
        T sum = 0;
        std::unordered_set<int> links_i(route_i.begin(), route_i.end());
        
        for(int link_id : route_j) {
            if(links_i.count(link_id)) {
                sum += this->links()[link_id].DelayDer();
            }
        }
        return sum;
    }

    MatrixXd RoutesDelayColumn(int od_index) {
      const auto& routes = this->GetODPair(od_index).GetRoutes();
      MatrixXd delays(routes.size(), 1);
        
      for(size_t i = 0; i < routes.size(); ++i) {
        delays(i, 0) = this->CalculatePathCost(routes[i]);
      }
      return delays;
    }

    MatrixXd EColumn(int od_index) {
      const int size = this->GetODPair(od_index).GetRoutesCount();
      return MatrixXd::Constant(size, 1, 1.0);
    }    

    T TargetDelay(int od_index) {
        if(HasConstantRoute(od_index)) {
            return ConstantRouteDelay(od_index);
        }
        return WeightedAverageDelay(od_index);
    }

    bool HasConstantRoute(int od_index) const {
        const auto& routes = this->GetODPair(od_index).GetRoutes();
        return std::any_of(routes.begin(), routes.end(), [this](const auto& route) {
            return std::all_of(route.begin(), route.end(), [this](int link_id) {
                return this->network_.links()[link_id].power == 0;
            });
        });
    }

    T ConstantRouteDelay(int od_index) const {
        const auto& routes = this->GetODPair(od_index).GetRoutes();
        for(const auto& route : routes) {
            if(std::all_of(route.begin(), route.end(), [this](int link_id) {
                return this->network_.links()[link_id].power == 0;
            })) {
                return this->CalculatePathCost(route);
            }
        }
        return T(0);
    }

    T WeightedAverageDelay(int od_index) {
        MatrixXd J_inv = RoutesJacobiMatrix(od_index).inverse();
        MatrixXd e = EColumn(od_index);
        MatrixXd delays = RoutesDelayColumn(od_index);
        
        return (e.transpose() * J_inv * delays)(0, 0) / 
               (e.transpose() * J_inv * e)(0, 0);
    }
  };
}

#endif
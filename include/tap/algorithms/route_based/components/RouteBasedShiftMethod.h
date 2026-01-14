#ifndef ROUTE_BASED_SHIFT_METHOD_H
#define ROUTE_BASED_SHIFT_METHOD_H

#include <vector>
#include "../../../core/Network.h"

namespace TrafficAssignment {
  template <typename T>
  class RouteBasedShiftMethod {
  public:
    explicit RouteBasedShiftMethod(Network<T>& network) 
        : network_(network) {}

    virtual ~RouteBasedShiftMethod() = default;

    virtual std::vector<T> FlowShift(int od_index) = 0;

  protected:
    Network<T>& network_;

    // Helper methods for common operations
    std::vector<Link<T>>& links() { 
        return network_.mutable_links(); 
    }
    
    std::vector<OriginDestinationPair<T>>& od_pairs() { 
        return network_.mutable_od_pairs(); 
    }
    
    const OriginDestinationPair<T>& GetODPair(int index) const {
        return network_.od_pairs().at(index);
    }
    
    T CalculatePathCost(const std::vector<int>& path) const {
        T total = 0;
        for(int link_id : path) {
            total += network_.links()[link_id].Delay();
        }
        return total;
    }
  };
}

#endif
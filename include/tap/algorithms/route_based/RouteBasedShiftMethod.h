#ifndef ROUTE_BASED_SHIFT_METHOD_H
#define ROUTE_BASED_SHIFT_METHOD_H

#include <vector>
#include "../../core/Network.h"

namespace TrafficAssignment {
  /**
   * @brief Abstract interface for path-based flow shift methods.
   *
   * Given an OD pair with multiple routes, computes new flow allocations that
   * move toward user equilibrium by shifting flow from costlier to cheaper routes.
   * See Perederieieva et al. (2015) Section 4.4 for the shift method framework.
   *
   * @tparam T Numeric type for flow computations.
   */
  template <typename T>
  class RouteBasedShiftMethod {
  public:
    explicit RouteBasedShiftMethod(Network<T>& network)
        : network_(network) {}

    virtual ~RouteBasedShiftMethod() = default;

    /// @brief Computes new route flow allocations for the given OD pair.
    /// @param od_index Index of the OD pair in the network.
    /// @return New flow values for each route of the OD pair.
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
    
    /// @brief Computes total travel time along a path (sum of link delays).
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
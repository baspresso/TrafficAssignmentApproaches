#ifndef ROUTE_BASED_SHIFT_METHOD_H
#define ROUTE_BASED_SHIFT_METHOD_H

#include <vector>

namespace TrafficAssignment {
  template <typename T>
  class RouteBasedShiftMethod {
  public:
    RouteBasedShiftMethod(std::vector <Link <T>>& links, std::vector <OriginDestinationPair <T>>& origin_destination_pairs) :  
      links_(links), origin_destination_pairs_(origin_destination_pairs) {}

    virtual ~RouteBasedShiftMethod() = default;

    virtual std::vector <T> FlowShift(int od_pair_index) = 0;

  protected:
    std::vector <Link <T>>& links_;
    std::vector <OriginDestinationPair <T>>& origin_destination_pairs_;
  };
}

#endif
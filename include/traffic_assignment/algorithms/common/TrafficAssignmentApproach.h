#ifndef TRAFFIC_ASSIGNMENT_APPROACH_H
#define TRAFFIC_ASSIGNMENT_APPROACH_H

#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include "../../data/Link.h"
#include "../../data/OriginDestinationPair.h"
#include "../../utils/DataProcessor.h"
#include "../../utils/StatisticsRecorder.h"
#include "../../core/Network.h"

namespace TrafficAssignment {

  template <typename T>
  class TrafficAssignmentApproach {
  public:
    TrafficAssignmentApproach(Network<T>& network, T alpha = 1e-6)
      : network_(network),
        alpha_(alpha),
        statistics_recorder_(network) {}

    virtual ~TrafficAssignmentApproach() = default;
    
    virtual std::string GetApproachName() = 0;

    virtual void ComputeTrafficFlows() = 0;

    Network<T>& network() { return network_; }

    protected:
      Network<T>& network_;
      T alpha_;
      StatisticsRecorder <T> statistics_recorder_;
  };
}

#endif
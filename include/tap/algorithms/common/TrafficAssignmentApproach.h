#ifndef TRAFFIC_ASSIGNMENT_APPROACH_H
#define TRAFFIC_ASSIGNMENT_APPROACH_H

#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <cstddef>
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

    virtual void ComputeTrafficFlows(bool statistics_recording = false) = 0;

    // Optional runtime tuning hook for approaches that support route-search parallelism.
    virtual void SetRouteSearchThreadCount(std::size_t thread_count) {
      (void)thread_count;
    }

    virtual std::size_t GetRouteSearchThreadCount() const {
      return 1;
    }

    void SetOutputRoot(const std::string& root) {
      statistics_recorder_.SetOutputRoot(root);
    }

    virtual void Reset() {
      network_.Reset();
    }

    Network<T>& network() { return network_; }

    protected:
      Network<T>& network_;
      T alpha_;
      StatisticsRecorder <T> statistics_recorder_;
  };
}

#endif

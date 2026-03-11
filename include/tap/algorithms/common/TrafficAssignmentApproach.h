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

  /**
   * @brief Abstract base class for Traffic Assignment Problem (TAP) solvers.
   *
   * Solves the user equilibrium problem: min T(f) = sum_a integral_0^{f_a} t_a(x) dx
   * subject to flow conservation constraints (sum of route flows = OD demand).
   * See Perederieieva et al. (2015) Eq. 2, Bar-Gera (2010) Eq. 1.
   *
   * Subclasses implement specific algorithms:
   * - RouteBasedApproach: path-based with explicit route enumeration
   * - TapasApproach: bush-based TAPAS with paired alternative segments
   *
   * @tparam T Numeric type for flow computations.
   */
  template <typename T>
  class TrafficAssignmentApproach {
  public:
    TrafficAssignmentApproach(Network<T>& network, T alpha = 1e-6)
      : network_(network),
        alpha_(alpha),
        statistics_recorder_(network) {}

    virtual ~TrafficAssignmentApproach() = default;
    
    virtual std::string GetApproachName() = 0;

    /// @brief Runs the TAP algorithm to compute user equilibrium flows.
    /// @param statistics_recording If true, records convergence trace to CSV.
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

    /// @brief Resets the network to initial state (zero flows, no routes).
    virtual void Reset() {
      network_.Reset();
    }

    Network<T>& network() { return network_; }

    protected:
      Network<T>& network_;
      T alpha_;                                   ///< RGAP convergence threshold.
      StatisticsRecorder <T> statistics_recorder_; ///< Optional convergence trace recorder.
  };
}

#endif

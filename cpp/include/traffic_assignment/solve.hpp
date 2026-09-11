#pragma once

#include <traffic_assignment/network.hpp>
#include <traffic_assignment/results.hpp>

#include <memory>
#include <vector>

namespace traffic_assignment {

class TrafficAssignmentApproach {
 public:
  virtual ~TrafficAssignmentApproach();
  TrafficAssignmentApproach(const TrafficAssignmentApproach&) = delete;
  TrafficAssignmentApproach& operator=(const TrafficAssignmentApproach&) = delete;
  std::shared_ptr<Network> network() const;
  std::string GetApproachName() const;
  void ComputeTrafficFlows(bool statistics_recording = false);
  void Reset();
  void SetOutputRoot(const std::string& path);
  void SetRouteSearchThreadCount(std::size_t count);
  std::size_t GetRouteSearchThreadCount() const;

 protected:
  TrafficAssignmentApproach(std::shared_ptr<Network> network, const TapOptions& options);

 private:
  friend struct detail::Access;
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

class TapasApproach : public TrafficAssignmentApproach {
 public:
  TapasApproach(std::shared_ptr<Network> network, Real alpha = 1e-14L,
                int max_iterations = 200, Real mu = 0.5L, Real v = 0.25L);
};

class RouteBasedApproach : public TrafficAssignmentApproach {
 public:
  RouteBasedApproach(std::shared_ptr<Network> network, Real alpha = 1e-14L,
                     const std::string& shift_method = "NewtonStep",
                     std::size_t route_search_threads = 1, int max_iterations = 200,
                     int full_iteration_count = 3, int origin_iteration_count = 1,
                     Real ema_alpha = 0.7L);
};

std::shared_ptr<TrafficAssignmentApproach> MakeApproach(
    std::shared_ptr<Network> network, const TapOptions& options = {});

TapResult SolveTap(TrafficAssignmentApproach& approach, bool reset = false,
                   bool statistics_recording = false);
TapResult SolveTap(std::shared_ptr<Network> network, const TapOptions& options = {});

class BilevelCND {
 public:
  BilevelCND(std::shared_ptr<Network> network,
             std::shared_ptr<TrafficAssignmentApproach> approach,
             const std::vector<LinkConstraint>& constraints,
             const std::vector<StepConfig>& pipeline,
             const CndpOptions& options = {});
  ~BilevelCND();
  BilevelCND(const BilevelCND&) = delete;
  BilevelCND& operator=(const BilevelCND&) = delete;

  void SetVerbose(bool enabled);
  void SetStatisticsEnabled(bool enabled);
  void SetFinalDiagnosticsEnabled(bool enabled);
  void SetRouteSearchThreadCount(std::size_t count);
  std::size_t active_constraint_count() const;
  CndpResult ComputeNetworkDesign();

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace traffic_assignment

#ifndef OPTIMIZATION_STEP_H
#define OPTIMIZATION_STEP_H

#include <string>

namespace TrafficAssignment {

struct OptimizationStepConfig {
  std::string type;
  std::string name;
  int max_iterations = 100;
  double tolerance = 1e-4;

  // NLopt-specific
  std::string algorithm;
  std::string local_algorithm;
  int local_max_iterations = 0;
  double local_tolerance = 0.0;
};

struct StepResult {
  double objective = 0.0;
  double total_travel_time = 0.0;
  double budget = 0.0;
  bool success = true;
  int evaluations = 0;
};

template <typename T>
class CndOptimizationContext;

template <typename T>
class OptimizationStep {
public:
  virtual ~OptimizationStep() = default;
  virtual StepResult Execute(CndOptimizationContext<T>& ctx) = 0;
  virtual std::string GetName() const = 0;
};

}  // namespace TrafficAssignment

#endif  // OPTIMIZATION_STEP_H

#ifndef OPTIMIZATION_PIPELINE_H
#define OPTIMIZATION_PIPELINE_H

#include <vector>
#include <memory>
#include <string>
#include <iostream>
#include <stdexcept>
#include "OptimizationStep.h"
#include "CndOptimizationContext.h"
#include "OptimizationStepFactory.h"

namespace TrafficAssignment {

/**
 * @brief Sequential executor for CNDP optimization steps.
 *
 * Runs [step.0], [step.1], ... in order, with each step receiving the network state
 * (link capacities) left by the previous step as a warm start. After each step,
 * performs a consistency TA solve and logs a post-step quality point.
 *
 * @tparam T Numeric type for flow/capacity computations.
 */
template <typename T>
class OptimizationPipeline {
public:
  OptimizationPipeline() = default;

  /// @brief Appends an optimization step to the pipeline.
  void AddStep(std::unique_ptr<OptimizationStep<T>> step) {
    steps_.push_back(std::move(step));
  }

  /// @brief Constructs a pipeline from a list of step configurations using OptimizationStepFactory.
  static OptimizationPipeline<T> BuildFromConfigs(
      const std::vector<OptimizationStepConfig>& configs) {
    OptimizationPipeline<T> pipeline;
    for (const auto& config : configs) {
      pipeline.AddStep(OptimizationStepFactory<T>::Create(config));
    }
    return pipeline;
  }

  /// @brief Executes all steps sequentially, with post-step TA consistency checks and logging.
  void Execute(CndOptimizationContext<T>& ctx) {
    for (std::size_t i = 0; i < steps_.size(); ++i) {
      auto& step = steps_[i];
      std::cout << step->GetName() << "-start\n";
      StepResult result = step->Execute(ctx);

      // Post-step: TA computation for consistency + log quality time point
      double post_step_ta_seconds = 0.0;
      if (!ctx.ComputeTrafficFlowsSafely(post_step_ta_seconds, "post_step_" + std::to_string(i))) {
        throw std::runtime_error("TA computation failed after step '" + step->GetName() + "'.");
      }
      const double post_total_travel_time = ctx.network.TotalTravelTime();
      const double post_budget = static_cast<double>(ctx.Budget());
      const double post_objective = post_total_travel_time + post_budget;
      if (!CndOptimizationContext<T>::IsFiniteNumber(post_total_travel_time) ||
          !CndOptimizationContext<T>::IsFiniteNumber(post_budget) ||
          !CndOptimizationContext<T>::IsFiniteNumber(post_objective)) {
        ctx.PrintInvalidStateWarning("post_step_" + std::to_string(i), "objective_invalid");
        throw std::runtime_error("Post-step objective contains non-finite values after step '" +
                                 step->GetName() + "'.");
      }
      ctx.UpdateRuntimeAndBudgetMetrics(post_step_ta_seconds, post_budget);
      ctx.LogQualityTimePoint("post_" + step->GetName(),
                              0,
                              post_objective,
                              post_total_travel_time,
                              post_budget,
                              post_step_ta_seconds,
                              true);
    }
  }

  std::size_t StepCount() const { return steps_.size(); }

private:
  std::vector<std::unique_ptr<OptimizationStep<T>>> steps_;
};

}  // namespace TrafficAssignment

#endif  // OPTIMIZATION_PIPELINE_H

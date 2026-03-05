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

template <typename T>
class OptimizationPipeline {
public:
  OptimizationPipeline() = default;

  void AddStep(std::unique_ptr<OptimizationStep<T>> step) {
    steps_.push_back(std::move(step));
  }

  static OptimizationPipeline<T> BuildFromConfigs(
      const std::vector<OptimizationStepConfig>& configs) {
    OptimizationPipeline<T> pipeline;
    for (const auto& config : configs) {
      pipeline.AddStep(OptimizationStepFactory<T>::Create(config));
    }
    return pipeline;
  }

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

#ifndef OPTIMIZATION_STEP_FACTORY_H
#define OPTIMIZATION_STEP_FACTORY_H

#include <memory>
#include <string>
#include <algorithm>
#include <stdexcept>
#include "OptimizationStep.h"
#include "steps/NloptOptimizationStep.h"
#include "steps/OptimalityConditionStep.h"

namespace TrafficAssignment {

template <typename T>
class OptimizationStepFactory {
public:
  static std::unique_ptr<OptimizationStep<T>> Create(const OptimizationStepConfig& config) {
    std::string type = config.type;
    std::transform(type.begin(), type.end(), type.begin(),
      [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (type == "nlopt") {
      return std::make_unique<NloptOptimizationStep<T>>(config);
    }
    if (type == "optimality_condition") {
      return std::make_unique<OptimalityConditionStep<T>>(config);
    }
    throw std::runtime_error("Unknown optimization step type: '" + config.type + "'");
  }
};

}  // namespace TrafficAssignment

#endif  // OPTIMIZATION_STEP_FACTORY_H

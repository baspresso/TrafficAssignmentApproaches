#ifndef OPTIMIZATION_STEP_FACTORY_H
#define OPTIMIZATION_STEP_FACTORY_H

#include <memory>
#include <string>
#include <algorithm>
#include <stdexcept>
#include "OptimizationStep.h"
#include "steps/NloptOptimizationStep.h"
#include "steps/OptimalityConditionStep.h"
#include "steps/OptimlibOptimizationStep.h"
#include "steps/GradientDescentStep.h"

namespace TrafficAssignment {

/**
 * @brief Factory that creates OptimizationStep instances from configuration.
 *
 * Maps step type strings to concrete step classes:
 * - "nlopt" -> NloptOptimizationStep (derivative-free: COBYLA, BOBYQA, ISRES, ...)
 * - "optimality_condition" -> OptimalityConditionStep (sensitivity-based, Chiou 2005)
 * - "optimlib" -> OptimlibOptimizationStep (population-based: DE, PSO, NM, ...)
 * - "gradient_descent" -> GradientDescentStep (projected GD with finite-difference gradients)
 *
 * @tparam T Numeric type for flow/capacity computations.
 */
template <typename T>
class OptimizationStepFactory {
public:
  /// @brief Creates a step instance from config. Throws on unknown type.
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
    if (type == "optimlib") {
      return std::make_unique<OptimlibOptimizationStep<T>>(config);
    }
    if (type == "gradient_descent") {
      return std::make_unique<GradientDescentStep<T>>(config);
    }
    throw std::runtime_error("Unknown optimization step type: '" + config.type + "'");
  }
};

}  // namespace TrafficAssignment

#endif  // OPTIMIZATION_STEP_FACTORY_H

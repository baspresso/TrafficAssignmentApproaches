#pragma once

#include <traffic_assignment/options.hpp>

#include <cmath>
#include <stdexcept>
#include <string>

namespace traffic_assignment::detail {

template <typename T>
void ValidateNonnegative(T value, const std::string& field, bool positive = false) {
  if (!std::isfinite(value) || (positive ? value <= 0 : value < 0)) {
    throw std::invalid_argument(field + (positive
        ? " must be finite and strictly positive" : " must be finite and nonnegative"));
  }
}

inline void ValidateConstraint(const LinkConstraint& constraint, const std::string& field) {
  ValidateNonnegative(constraint.lower_bound, field + ".lower_bound", true);
  ValidateNonnegative(constraint.upper_bound, field + ".upper_bound", true);
  ValidateNonnegative(constraint.investment_cost_param, field + ".investment_cost_param");
  if (constraint.lower_bound > constraint.upper_bound) {
    throw std::invalid_argument(field + ".lower_bound must be less than or equal to upper_bound");
  }
}

}  // namespace traffic_assignment::detail

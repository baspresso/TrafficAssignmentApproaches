#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace traffic_assignment::detail::csv {

inline std::string Trim(const std::string& value) {
  const auto first = value.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return {};
  return value.substr(first, value.find_last_not_of(" \t\r\n") - first + 1);
}

inline std::vector<std::string> Split(const std::string& line, char delimiter = ',') {
  std::vector<std::string> fields;
  std::size_t start = 0;
  while (true) {
    const auto end = line.find(delimiter, start);
    fields.push_back(Trim(line.substr(start, end == std::string::npos ? end : end - start)));
    if (end == std::string::npos) break;
    start = end + 1;
  }
  return fields;
}

inline double Number(const std::string& token, const std::string& field) {
  std::size_t consumed = 0;
  double value;
  try {
    value = std::stod(token, &consumed);
  } catch (const std::exception&) {
    throw std::invalid_argument(field + " must be a finite number");
  }
  if (consumed != token.size() || !std::isfinite(value)) {
    throw std::invalid_argument(field + " must be a finite number");
  }
  return value;
}

inline int Integer(const std::string& token, const std::string& field, int minimum = 0) {
  const auto value = Number(token, field);
  if (value != std::floor(value) || value < minimum || value > std::numeric_limits<int>::max()) {
    throw std::invalid_argument(field + " must be an integer in [" +
        std::to_string(minimum) + ", " + std::to_string(std::numeric_limits<int>::max()) + "]");
  }
  return static_cast<int>(value);
}

template <typename T>
std::vector<T> OrderByLinkIndex(std::vector<T> values, const std::vector<int>& indices) {
  if (indices.size() != values.size()) throw std::invalid_argument("one link_index required per row");
  auto ordered = values;
  std::vector<bool> seen(values.size(), false);
  for (std::size_t i = 0; i < indices.size(); ++i) {
    const auto position = static_cast<std::size_t>(indices[i]);
    if (position >= values.size() || seen[position]) {
      throw std::invalid_argument("link_index must contain every index in [0, " +
          std::to_string(values.size()) + ") exactly once");
    }
    seen[position] = true;
    ordered[position] = std::move(values[i]);
  }
  return ordered;
}

}  // namespace traffic_assignment::detail::csv

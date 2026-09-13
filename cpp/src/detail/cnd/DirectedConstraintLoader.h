#pragma once

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "LinkConstraint.h"
#include "../csv.hpp"
#include "../input_validation.hpp"

namespace TrafficAssignment {

/// Files without link_index retain row order. Native endpoints are zero-based.
class DirectedConstraintLoader {
 public:
  void SetVerbose(bool verbose) { verbose_ = verbose; }

  std::vector<DirectedLinkCapacityConstraint> LoadFromFile(
      const std::string& filepath, int node_index_base = 1) const {
    namespace csv = traffic_assignment::detail::csv;
    if (node_index_base != 0 && node_index_base != 1) {
      throw std::invalid_argument("node_index_base must be 0 or 1");
    }
    std::ifstream file(filepath);
    if (!file) throw std::runtime_error("Cannot open constraint file: " + filepath);

    std::vector<DirectedLinkCapacityConstraint> constraints;
    std::vector<int> indices;
    std::map<std::string, int> columns;
    int init = -1, term = -1, lower = -1, upper = -1, cost = -1, index = -1;
    char delimiter = ',';
    std::size_t line_number = 0, column_count = 0;
    std::string line;
    try {
      while (std::getline(file, line)) {
        ++line_number;
        line = csv::Trim(line);
        if (line.empty() || line.front() == '#') continue;
        if (columns.empty()) {
          for (char candidate : {'\t', ';', '|', ','}) {
            if (line.find(candidate) != std::string::npos) {
              delimiter = candidate;
              break;
            }
          }
          const auto header = csv::Split(line, delimiter);
          column_count = header.size();
          for (std::size_t i = 0; i < header.size(); ++i) {
            auto name = header[i];
            std::transform(name.begin(), name.end(), name.begin(),
                           [](unsigned char ch) { return std::tolower(ch); });
            if (name.empty() || !columns.emplace(name, static_cast<int>(i)).second) {
              throw std::invalid_argument("constraint columns must have unique, nonempty names");
            }
          }
          auto column = [&](std::initializer_list<const char*> aliases, bool required = true) {
            for (const auto* alias : aliases) {
              if (const auto found = columns.find(alias); found != columns.end()) return found->second;
            }
            if (required) throw std::invalid_argument(std::string("missing required column: ") + *aliases.begin());
            return -1;
          };
          init = column({"init_node", "from_node", "source", "origin", "i", "init"});
          term = column({"term_node", "to_node", "destination", "sink", "j", "term"});
          lower = column({"lower_bound", "lb", "lower", "min_capacity", "min"});
          upper = column({"upper_bound", "ub", "upper", "max_capacity", "max"});
          cost = column({"investment_cost_param", "cost_param", "cost", "investment_cost", "unit_cost"}, false);
          index = column({"link_index"}, false);
          continue;
        }
        const auto fields = csv::Split(line, delimiter);
        if (fields.size() != column_count) {
          throw std::invalid_argument("constraint row must have one field per header column");
        }
        DirectedLinkCapacityConstraint constraint(
            csv::Integer(fields[init], "init_node", node_index_base) - node_index_base,
            csv::Integer(fields[term], "term_node", node_index_base) - node_index_base,
            csv::Number(fields[lower], "lower_bound"),
            csv::Number(fields[upper], "upper_bound"),
            cost < 0 ? 1.0 : csv::Number(fields[cost], "investment_cost_param"));
        traffic_assignment::detail::ValidateConstraint(constraint, "constraint");
        constraints.push_back(constraint);
        if (index >= 0) indices.push_back(csv::Integer(fields[index], "link_index"));
      }
      if (columns.empty()) throw std::invalid_argument("constraints file must contain a header");
      if (constraints.empty()) throw std::invalid_argument("constraints must not be empty");
      if (index >= 0) {
        constraints = csv::OrderByLinkIndex(std::move(constraints), indices);
      }
    } catch (const std::exception& error) {
      throw std::runtime_error(filepath + ": line " + std::to_string(line_number) + ": " + error.what());
    }
    if (verbose_) std::cout << "Loaded " << constraints.size() << " directed constraints from " << filepath << '\n';
    return constraints;
  }

 private:
  bool verbose_ = false;
};

}  // namespace TrafficAssignment

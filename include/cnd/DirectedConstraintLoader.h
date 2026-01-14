#ifndef DIRECTED_CONSTRAINT_LOADER_H
#define DIRECTED_CONSTRAINT_LOADER_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <map>
#include <algorithm>
#include <optional>
#include <iomanip>
#include "BilevelCND.h"

namespace TrafficAssignment {

/**
 * @class DirectedLinkCapacityConstraint
 * @brief Represents a constraint for a directed link (arc) in the network
 *
 * This structure represents capacity constraints for directed arcs,
 * identified by their origin (init_node) and destination (term_node) nodes.
 */
struct DirectedLinkCapacityConstraint {
  std::size_t init_node;               // Origin node
  std::size_t term_node;               // Destination node
  double lower_bound;                  // Minimum capacity
  double upper_bound;                  // Maximum capacity
  double investment_cost_param;        // Cost per unit capacity

  DirectedLinkCapacityConstraint() = default;

  DirectedLinkCapacityConstraint(
    std::size_t origin,
    std::size_t destination,
    double lb,
    double ub,
    double cost)
    : init_node(origin),
      term_node(destination),
      lower_bound(lb),
      upper_bound(ub),
      investment_cost_param(cost) {}

  /**
   * @brief Get string representation for logging
   */
  std::string ToString() const {
    return std::string("(") + std::to_string(init_node) + "->" +
           std::to_string(term_node) + ")";
  }

  /**
   * @brief Get unique arc ID (for compatibility with undirected networks)
   * Format: "init_node_term_node"
   */
  std::string GetArcID() const {
    return std::to_string(init_node) + "_" + std::to_string(term_node);
  }
};

/**
 * @class DirectedConstraintCSVFormat
 * @brief Defines the expected CSV format for directed constraint files
 */
struct DirectedConstraintCSVFormat {
  enum class Delimiter {
    COMMA = ',',
    SEMICOLON = ';',
    TAB = '\t',
    PIPE = '|'
  };

  Delimiter delimiter = Delimiter::COMMA;
  bool has_header = true;
  bool skip_empty_lines = true;

  // Column indices for directed network
  int col_init_node = 0;
  int col_term_node = 1;
  int col_lower_bound = 2;
  int col_upper_bound = 3;
  int col_investment_cost = 4;

  // Optional comment character
  char comment_char = '#';
  bool skip_comment_lines = true;
};

/**
 * @class DirectedConstraintLoader
 * @brief Loads directed link capacity constraints from CSV files
 *
 * Adapted for networks with directed arcs (init_node → term_node).
 * Features:
 * - Multiple CSV format support
 * - Flexible column ordering via header parsing
 * - Comprehensive error handling and validation
 * - Logging and statistics reporting
 * - Conversion to/from LinkCapacityConstraint for compatibility
 *
 * Example CSV format:
 *
 *   init_node,term_node,lower_bound,upper_bound,investment_cost_param
 *   1,2,23310.180576,28490.220704000003,0
 *   1,3,21063.125871,25743.820509,0
 *   2,1,23310.180576,28490.220704000003,0
 *   2,6,4462.3628352,5453.9990208,0
 *   3,1,21063.125871,25743.820509,0
 */
class DirectedConstraintLoader {
public:
  /**
   * @brief Loading statistics for reporting
   */
  struct LoadingStatistics {
    std::size_t total_lines_read = 0;
    std::size_t header_lines = 0;
    std::size_t comment_lines = 0;
    std::size_t empty_lines = 0;
    std::size_t valid_constraints = 0;
    std::size_t validation_errors = 0;
    std::size_t parsing_errors = 0;

    void Print(std::ostream& os = std::cout) const {
      os << "=== Directed Constraint Loading Statistics ===" << std::endl;
      os << "Total lines read:       " << total_lines_read << std::endl;
      os << "Header lines:           " << header_lines << std::endl;
      os << "Comment lines:          " << comment_lines << std::endl;
      os << "Empty lines:            " << empty_lines << std::endl;
      os << "Valid constraints:      " << valid_constraints << std::endl;
      os << "Validation errors:      " << validation_errors << std::endl;
      os << "Parsing errors:         " << parsing_errors << std::endl;
    }
  };

  /**
   * @brief Constructor with default format
   */
  DirectedConstraintLoader() : format_(), statistics_() {}

  /**
   * @brief Constructor with custom format
   */
  explicit DirectedConstraintLoader(const DirectedConstraintCSVFormat& format)
    : format_(format), statistics_() {}

  /**
   * @brief Set the CSV format specification
   */
  void SetFormat(const DirectedConstraintCSVFormat& format) {
    format_ = format;
  }

  /**
   * @brief Get the current CSV format specification
   */
  const DirectedConstraintCSVFormat& GetFormat() const {
    return format_;
  }

  /**
   * @brief Get loading statistics from last load operation
   */
  const LoadingStatistics& GetStatistics() const {
    return statistics_;
  }

  /**
   * @brief Enable/disable verbose output during loading
   */
  void SetVerbose(bool verbose) {
    verbose_ = verbose;
  }

  /**
   * @brief Load directed constraints from CSV file
   *
   * Automatically detects CSV delimiter and column order from header.
   *
   * @param filepath Path to CSV file
   * @return Vector of DirectedLinkCapacityConstraint objects
   * @throws std::runtime_error if file cannot be opened or parsing fails
   */
  std::vector<DirectedLinkCapacityConstraint> LoadFromFile(
    const std::string& filepath) {
    
    std::ifstream file(filepath);
    if (!file.is_open()) {
      throw std::runtime_error("Cannot open constraint file: " + filepath);
    }

    if (verbose_) {
      std::cout << "Loading directed constraints from: " << filepath << std::endl;
    }

    statistics_ = LoadingStatistics();

    std::vector<DirectedLinkCapacityConstraint> constraints;
    std::string line;
    int line_number = 0;

    try {
      // Read first line to detect format
      if (!std::getline(file, line)) {
        throw std::runtime_error("File is empty");
      }

      line_number = 1;
      statistics_.total_lines_read++;

      // Auto-detect delimiter
      auto actual_delimiter = format_.delimiter;
      if (!IsValidDelimiter(line)) {
        actual_delimiter = DetectDelimiter(line);
        if (verbose_) {
          std::cout << "Auto-detected delimiter: '" 
                    << static_cast<char>(actual_delimiter) << "'" << std::endl;
        }
      }

      // Parse header if present
      if (format_.has_header) {
        statistics_.header_lines++;
        ParseHeaderLine(line, static_cast<char>(actual_delimiter));
        if (verbose_) {
          std::cout << "Header parsed: init_node[" << format_.col_init_node
                    << "], term_node[" << format_.col_term_node
                    << "], lb[" << format_.col_lower_bound
                    << "], ub[" << format_.col_upper_bound
                    << "], cost[" << format_.col_investment_cost << "]" << std::endl;
        }
      } else {
        // Parse first line as data if no header
        auto constraint = ParseConstraintLine(
          line, static_cast<char>(actual_delimiter), line_number);
        if (constraint.has_value()) {
          constraints.push_back(constraint.value());
          statistics_.valid_constraints++;
        }
      }

      // Parse remaining lines
      while (std::getline(file, line)) {
        line_number++;
        statistics_.total_lines_read++;

        if (format_.skip_empty_lines && line.empty()) {
          statistics_.empty_lines++;
          continue;
        }

        line = TrimWhitespace(line);

        if (format_.skip_comment_lines && !line.empty() &&
            line[0] == format_.comment_char) {
          statistics_.comment_lines++;
          continue;
        }

        if (line.empty()) {
          statistics_.empty_lines++;
          continue;
        }

        auto constraint = ParseConstraintLine(
          line, static_cast<char>(actual_delimiter), line_number);
        if (constraint.has_value()) {
          constraints.push_back(constraint.value());
          statistics_.valid_constraints++;
        } else {
          statistics_.parsing_errors++;
        }
      }

      file.close();

      if (verbose_) {
        std::cout << "Successfully loaded " << statistics_.valid_constraints
                  << " directed constraints" << std::endl;
      }

      if (statistics_.parsing_errors > 0) {
        std::cout << "WARNING: " << statistics_.parsing_errors
                  << " lines had parsing errors" << std::endl;
      }

      if (statistics_.validation_errors > 0) {
        std::cout << "WARNING: " << statistics_.validation_errors
                  << " constraints failed validation" << std::endl;
      }

      return constraints;

    } catch (const std::exception& e) {
      file.close();
      throw std::runtime_error(
        "Error loading constraints at line " + std::to_string(line_number) +
        ": " + std::string(e.what()));
    }
  }

  /**
   * @brief Directly add constraints to BilevelCND from CSV file
   *
   * Convenience method that loads constraints and adds them to solver.
   *
   * @param filepath Path to CSV constraint file
   * @param arc_link_mapping Optional mapping from arc (init,term) to link_id
   * @throws std::runtime_error if loading or validation fails
   */
  void LoadAndAddConstraints(
    const std::string& filepath,
    //BilevelCND<double>& solver,
    const std::map<std::pair<std::size_t, std::size_t>, std::size_t>*
      arc_link_mapping = nullptr) {
    
    auto constraints = LoadFromFile(filepath);
  }

  static std::size_t GenerateArcID(std::size_t init_node, std::size_t term_node) {
     // Simple approach for typical networks (nodes < 1000)
     return init_node * 10000 + term_node;
  }

private:
  DirectedConstraintCSVFormat format_;
  LoadingStatistics statistics_;
  bool verbose_ = false;

  /**
   * @brief Trim leading and trailing whitespace from string
   */
  static std::string TrimWhitespace(const std::string& str) {
    const auto start = str.find_first_not_of(" \t\n\r");
    if (start == std::string::npos) return "";
    const auto end = str.find_last_not_of(" \t\n\r");
    return str.substr(start, end - start + 1);
  }

  /**
   * @brief Split string by delimiter
   */
  static std::vector<std::string> Split(
    const std::string& line,
    char delimiter) {
    
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string field;

    while (std::getline(ss, field, delimiter)) {
      fields.push_back(TrimWhitespace(field));
    }

    return fields;
  }

  /**
   * @brief Auto-detect CSV delimiter from line
   */
  static DirectedConstraintCSVFormat::Delimiter DetectDelimiter(
    const std::string& line) {
    
    int comma_count = 0, semicolon_count = 0, tab_count = 0, pipe_count = 0;

    for (char c : line) {
      if (c == ',') comma_count++;
      else if (c == ';') semicolon_count++;
      else if (c == '\t') tab_count++;
      else if (c == '|') pipe_count++;
    }

    if (tab_count > 0) return DirectedConstraintCSVFormat::Delimiter::TAB;
    if (semicolon_count > 0)
      return DirectedConstraintCSVFormat::Delimiter::SEMICOLON;
    if (pipe_count > 0) return DirectedConstraintCSVFormat::Delimiter::PIPE;
    return DirectedConstraintCSVFormat::Delimiter::COMMA;
  }

  /**
   * @brief Check if delimiter is valid for the line
   */
  bool IsValidDelimiter(const std::string& line) const {
    char delim = static_cast<char>(format_.delimiter);
    return line.find(delim) != std::string::npos;
  }

  /**
   * @brief Parse header line and identify column positions
   */
  void ParseHeaderLine(const std::string& line, char delimiter) {
    auto fields = Split(line, delimiter);

    std::map<std::string, int> column_map;
    for (std::size_t i = 0; i < fields.size(); ++i) {
      std::string field = fields[i];
      std::transform(field.begin(), field.end(), field.begin(), ::tolower);
      column_map[field] = static_cast<int>(i);
    }

    // Map header names to column indices
    std::vector<std::string> init_node_aliases = {"init_node", "from_node", "source",
                                                   "origin", "i", "init"};
    std::vector<std::string> term_node_aliases = {"term_node", "to_node", "destination",
                                                   "sink", "j", "term"};
    std::vector<std::string> lb_aliases = {"lower_bound", "lb", "lower",
                                            "min_capacity", "min"};
    std::vector<std::string> ub_aliases = {"upper_bound", "ub", "upper",
                                            "max_capacity", "max"};
    std::vector<std::string> cost_aliases = {"investment_cost_param", "cost_param",
                                              "cost", "investment_cost", "unit_cost"};

    for (const auto& alias : init_node_aliases) {
      if (column_map.find(alias) != column_map.end()) {
        format_.col_init_node = column_map[alias];
        break;
      }
    }

    for (const auto& alias : term_node_aliases) {
      if (column_map.find(alias) != column_map.end()) {
        format_.col_term_node = column_map[alias];
        break;
      }
    }

    for (const auto& alias : lb_aliases) {
      if (column_map.find(alias) != column_map.end()) {
        format_.col_lower_bound = column_map[alias];
        break;
      }
    }

    for (const auto& alias : ub_aliases) {
      if (column_map.find(alias) != column_map.end()) {
        format_.col_upper_bound = column_map[alias];
        break;
      }
    }

    for (const auto& alias : cost_aliases) {
      if (column_map.find(alias) != column_map.end()) {
        format_.col_investment_cost = column_map[alias];
        break;
      }
    }
  }

  /**
   * @brief Parse a single constraint data line
   */
  std::optional<DirectedLinkCapacityConstraint> ParseConstraintLine(
    const std::string& line,
    char delimiter,
    int line_number) {
    
    auto fields = Split(line, delimiter);

    int required_cols = std::max({format_.col_init_node, format_.col_term_node,
                                  format_.col_lower_bound, format_.col_upper_bound,
                                  format_.col_investment_cost}) +
                       1;

    if (static_cast<int>(fields.size()) < required_cols) {
      if (verbose_) {
        std::cerr << "ERROR Line " << line_number << ": Expected at least "
                  << required_cols << " columns, got " << fields.size()
                  << std::endl;
      }
      return std::nullopt;
    }

    try {
      std::size_t init_node = std::stoul(fields[format_.col_init_node]);
      std::size_t term_node = std::stoul(fields[format_.col_term_node]);
      double lower_bound = std::stod(fields[format_.col_lower_bound]);
      double upper_bound = std::stod(fields[format_.col_upper_bound]);
      double investment_cost = std::stod(fields[format_.col_investment_cost]);

      DirectedLinkCapacityConstraint constraint(init_node, term_node, lower_bound,
                                                upper_bound, investment_cost);

      return constraint;

    } catch (const std::exception& e) {
      if (verbose_) {
        std::cerr << "ERROR Line " << line_number
                  << ": Failed to parse fields - " << e.what() << std::endl;
        std::cerr << "  Fields: ";
        for (const auto& f : fields) {
          std::cerr << "[" << f << "] ";
        }
        std::cerr << std::endl;
      }
      return std::nullopt;
    }
  }
};

} // namespace TrafficAssignment

#endif // DIRECTED_CONSTRAINT_LOADER_H
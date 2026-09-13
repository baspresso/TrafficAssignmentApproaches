#pragma once

#include <traffic_assignment/network.hpp>

#include <filesystem>
#include <fstream>
#include <map>
#include <tuple>
#include <vector>

#include "../data/Link.h"
#include "../../csv.hpp"

namespace TrafficAssignment {
namespace fs = std::filesystem;

/// Parses preprocessed CSV files into values for the public Network constructor.
class NetworkBuilder {
 public:
  static constexpr const char* kDefaultDataRoot = "data/TransportationNetworks";

  std::tuple<int, int, std::vector<traffic_assignment::LinkData>> LoadNetworkData(
      const std::string& dataset, const fs::path& root = kDefaultDataRoot) const {
    namespace csv = traffic_assignment::detail::csv;
    const auto path = root / dataset / (dataset + "_net.csv");
    std::ifstream file(path);
    if (!file) throw std::runtime_error("Cannot open network file: " + path.string());
    std::size_t line_number = 1;
    try {
      std::string line;
      if (!std::getline(file, line)) throw std::invalid_argument("missing network metadata");
      const auto metadata = ParseMetadata(line);
      const int zones = RequiredMetadata(metadata, "NUMBER OF ZONES");
      const int nodes = RequiredMetadata(metadata, "NUMBER OF NODES");
      const int count = RequiredMetadata(metadata, "NUMBER OF LINKS");
      if (zones <= 0 || nodes < zones) throw std::invalid_argument("need 0 < n_zones <= n_nodes");

      ++line_number;
      if (!std::getline(file, line)) throw std::invalid_argument("missing network header");
      const auto header = csv::Split(line);
      std::map<std::string, std::size_t> columns;
      for (std::size_t i = 0; i < header.size(); ++i) {
        if (header[i].empty() || !columns.emplace(header[i], i).second) {
          throw std::invalid_argument("network columns must have unique, nonempty names");
        }
      }
      for (const auto* field : {"init_node", "term_node", "capacity", "free_flow_time", "b", "power"}) {
        if (!columns.contains(field)) throw std::invalid_argument(std::string("missing required column: ") + field);
      }
      std::vector<traffic_assignment::LinkData> links;
      std::vector<int> indices;
      while (std::getline(file, line)) {
        ++line_number;
        line = csv::Trim(line);
        if (line.empty() || line.front() == '#') continue;
        const auto fields = csv::Split(line);
        if (fields.size() != header.size()) throw std::invalid_argument("network row must have one field per header column");
        auto number = [&](const char* field, double fallback = 0) {
          const auto found = columns.find(field);
          return found == columns.end() ? fallback : csv::Number(fields[found->second], field);
        };
        auto integer = [&](const char* field, int fallback, int minimum = 0) {
          const auto found = columns.find(field);
          return found == columns.end() ? fallback : csv::Integer(fields[found->second], field, minimum);
        };
        links.push_back({integer("init_node", 0, 1) - 1, integer("term_node", 0, 1) - 1,
                         number("capacity"), number("length"), number("free_flow_time"),
                         number("b"), number("power"), number("speed"), number("toll"),
                         integer("link_type", 1)});
        if (columns.contains("link_index")) indices.push_back(integer("link_index", 0));
      }
      if (links.size() != static_cast<std::size_t>(count)) {
        throw std::invalid_argument("network row count must match NUMBER OF LINKS");
      }
      if (columns.contains("link_index")) links = csv::OrderByLinkIndex(std::move(links), indices);
      return {nodes, zones, std::move(links)};
    } catch (const std::exception& error) {
      throw std::runtime_error(path.string() + ": line " + std::to_string(line_number) + ": " + error.what());
    }
  }

  template <typename T>
  std::vector<std::vector<T>> LoadTripData(const std::string& dataset, int zones,
                                         const fs::path& root = kDefaultDataRoot) const {
    namespace csv = traffic_assignment::detail::csv;
    const auto path = root / dataset / (dataset + "_trips.csv");
    std::ifstream file(path);
    if (!file) throw std::runtime_error("Cannot open trip file: " + path.string());
    std::size_t line_number = 1;
    try {
      std::string line;
      if (!std::getline(file, line)) throw std::invalid_argument("missing demand metadata");
      const auto metadata = ParseMetadata(line);
      if (zones <= 0 || RequiredMetadata(metadata, "NUMBER OF ZONES") != zones) {
        throw std::invalid_argument("demand NUMBER OF ZONES must match the network");
      }
      std::vector<std::vector<T>> trips;
      while (std::getline(file, line)) {
        ++line_number;
        line = csv::Trim(line);
        if (line.empty() || line.front() == '#') continue;
        const auto fields = csv::Split(line);
        if (fields.size() != static_cast<std::size_t>(zones)) {
          throw std::invalid_argument("demand must have one column per zone");
        }
        std::vector<T> row;
        row.reserve(fields.size());
        for (const auto& field : fields) row.push_back(static_cast<T>(csv::Number(field, "demand")));
        trips.push_back(std::move(row));
      }
      if (trips.size() != static_cast<std::size_t>(zones)) {
        throw std::invalid_argument("demand must have one row per zone");
      }
      return trips;
    } catch (const std::exception& error) {
      throw std::runtime_error(path.string() + ": line " + std::to_string(line_number) + ": " + error.what());
    }
  }

  template <typename T>
  std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>> BuildAdjacencyLists(
      const std::vector<Link<T>>& links, int node_count) const {
    std::vector<std::vector<int>> adjacency(node_count), reverse(node_count);
    for (std::size_t i = 0; i < links.size(); ++i) {
      adjacency[links[i].init].push_back(static_cast<int>(i));
      reverse[links[i].term].push_back(static_cast<int>(i));
    }
    return {std::move(adjacency), std::move(reverse)};
  }

 private:
  static std::map<std::string, int> ParseMetadata(const std::string& line) {
    namespace csv = traffic_assignment::detail::csv;
    const auto text = csv::Trim(line);
    if (text.empty() || text.front() != '#') throw std::invalid_argument("metadata must start with #");
    std::map<std::string, int> result;
    for (const auto& item : csv::Split(text.substr(1))) {
      const auto fields = csv::Split(item, ':');
      if (fields.size() != 2 || !result.emplace(fields[0], csv::Integer(fields[1], fields[0])).second) {
        throw std::invalid_argument("invalid or duplicate metadata field");
      }
    }
    return result;
  }

  static int RequiredMetadata(const std::map<std::string, int>& metadata, const std::string& key) {
    const auto found = metadata.find(key);
    if (found == metadata.end()) throw std::invalid_argument("missing metadata: " + key);
    return found->second;
  }
};

}  // namespace TrafficAssignment

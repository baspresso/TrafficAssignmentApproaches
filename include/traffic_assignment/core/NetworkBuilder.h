// include/traffic_assignment/core/NetworkBuilder.h
#ifndef NETWORK_BUILDER_H
#define NETWORK_BUILDER_H

#include "Network.h"
#include <filesystem>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

namespace fs = std::filesystem;

namespace TrafficAssignment {

class NetworkBuilder {
public:
    NetworkBuilder() = default;

    template <typename T>
    Network<T> BuildFromDataset(const std::string& dataset_name) {
        std::string name = dataset_name;
        auto [nodes, zones, links] = LoadNetworkData<T>(dataset_name);
        auto trips = LoadTripData<T>(dataset_name, zones);
        auto [adjacency, reverse_adjacency] = BuildAdjacencyLists<T>(links, nodes);
        
        return Network<T>(name, nodes, zones, 
            links, trips, 
            adjacency, reverse_adjacency
        );
    }

private:
    template <typename T>
    std::tuple<int, int, std::vector<Link<T>>> 
    LoadNetworkData(const std::string& dataset_name) {
        namespace fs = std::filesystem;
        
        // Path construction
        fs::path net_path = GetDatasetPath(dataset_name) / (dataset_name + "_net.csv");
        std::ifstream net_file(net_path);
        
        if(!net_file.is_open()) {
            throw std::runtime_error("Cannot open network file: " + net_path.string());
        }

        // Read metadata
        int nodes = 0, zones = 0, link_count = 0;
        std::string line;
        std::getline(net_file, line);
        ParseNetworkMetadata(line, nodes, zones, link_count);
        // Skip line
        std::getline(net_file, line);
        // Read links
        std::vector<Link<T>> links;
        links.reserve(link_count);
        
        while(std::getline(net_file, line)) {
            links.push_back(ParseLinkLine<T>(line));
        }

        return {nodes, zones, links};
    }

    template <typename T>
    std::vector<std::vector<T>> LoadTripData(const std::string& dataset_name, int zones) {
        namespace fs = std::filesystem;
        
        fs::path trip_path = GetDatasetPath(dataset_name) / (dataset_name + "_trips.csv");
        std::ifstream trip_file(trip_path);
        
        if(!trip_file.is_open()) {
            throw std::runtime_error("Cannot open trip file: " + trip_path.string());
        }

        std::vector<std::vector<T>> trips(zones, std::vector<T>(zones, 0));
        std::string line;
        
        // Skip metadata line
        std::getline(trip_file, line);
        
        for(int origin = 0; origin < zones; ++origin) {
            std::getline(trip_file, line);
            std::stringstream ss(line);
            std::string value;
            
            for(int dest = 0; dest < zones; ++dest) {
                std::getline(ss, value, ',');
                trips[origin][dest] = ConvertValue<T>(value);
            }
        }

        return trips;
    }

    template <typename T>
    std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>> 
    BuildAdjacencyLists(const std::vector<Link<T>>& links, int node_count) {
        std::vector<std::vector<int>> adjacency(node_count);
        std::vector<std::vector<int>> reverse_adjacency(node_count);

        for(int i = 0; i < links.size(); ++i) {
            const auto& link = links[i];
            adjacency[link.init].push_back(i);
            reverse_adjacency[link.term].push_back(i);
        }

        return {adjacency, reverse_adjacency};
    }

    // Helper methods
    fs::path GetDatasetPath(const std::string& dataset_name) {
        static const fs::path base_path = "C:/Projects/TrafficAssignmentApproaches/data/TransportationNetworks";
        return base_path / dataset_name;
    }

    void ParseNetworkMetadata(const std::string& line, int& nodes, int& zones, int& links) {
        std::stringstream ss(line.substr(1)); // Skip '#'
        std::string token;
        
        std::getline(ss, token, ':'); // "NUMBER OF ZONES"
        ss >> zones;
        ss.ignore(2); // Skip ", "
        
        std::getline(ss, token, ':'); // "NUMBER OF NODES"
        ss >> nodes;
        ss.ignore(2);
        
        std::getline(ss, token, ':'); // "NUMBER OF LINKS"
        ss >> links;
    }

    template <typename T>
    Link<T> ParseLinkLine(const std::string& line) {
        std::stringstream ss(line);
        std::string token;
        int init, term, type;
        T capacity, length, free_flow_time, b, power, speed, toll;
        // Input files use 1-based indexing, convert to 0-based
        std::getline(ss, token, ','); init = std::stoi(token) - 1;
        std::getline(ss, token, ','); term = std::stoi(token) - 1;
        std::getline(ss, token, ','); capacity = std::stod(token);
        std::getline(ss, token, ','); length = std::stod(token);
        std::getline(ss, token, ','); free_flow_time = std::stod(token);
        std::getline(ss, token, ','); b = std::stod(token);
        std::getline(ss, token, ','); power = std::stod(token);
        std::getline(ss, token, ','); speed = std::stod(token);
        std::getline(ss, token, ','); toll = std::stod(token);
        std::getline(ss, token); type = std::stoi(token);


        return Link<T>(
            init, term, capacity, length, free_flow_time,
            b, power, speed, toll, type
        );
    }

    template <typename T>
    T ConvertValue(const std::string& str) {
        if constexpr (std::is_floating_point_v<T>) {
            return std::stod(str);
        } else {
            return static_cast<T>(std::stoll(str));
        }
    }
};

} // namespace TrafficAssignment

#endif // NETWORK_BUILDER_H
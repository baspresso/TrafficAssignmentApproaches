#ifndef DATA_PROCESSOR_H
#define DATA_PROCESSOR_H

#include "../data/Link.h"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <filesystem>

namespace TrafficAssignment {
  class DataProcessor {
  public:
    DataProcessor() = default;

    bool LoadData(const std::string& dataset_name) {
      //std::filesystem::path current_path = std::filesystem::current_path();
      std::filesystem::path current_path = "C:/Projects/TrafficAssignmentApproaches";
      //while (current_path.filename() != "out") {
      //  current_path = current_path.parent_path();
      //}
      //current_path = current_path.parent_path();
      auto net_file_path = current_path / "data" / "TransportationNetworks" / dataset_name / (dataset_name + "_net.csv");
      auto trips_file_path = current_path / "data/" / "TransportationNetworks" / dataset_name / (dataset_name + "_trips.csv");

      LoadNet(net_file_path);
      LoadTrips(trips_file_path);
      return true;
    }

    const int GetNumberOfZones() const {
      return number_of_zones_;
    }

    const int GetNumberOfNodes() const {
      return number_of_nodes_;
    }

    const int GetNumberOfLinks() const {
      return number_of_links_;
    }

    const std::vector <Link <long double>> GetLinks() const {
      return links_;
    }

    const std::vector <std::vector <long double>> GetTrips() const {
      return trips_;
    }

  private:
    int number_of_zones_, number_of_nodes_, number_of_links_;
    std::vector <Link <long double>> links_;
    std::vector <std::vector <long double>> trips_;

    bool LoadNet(const std::filesystem::path& net_file_path) {
      std::ifstream file(net_file_path);
      std::string line;


      if (file.is_open()) {
        // Reading metadata

        std::getline(file, line);
        std::stringstream metadata_ss(line);

        int vertices, links;
        std::string temp;

        if (line[0] == '#') {
          std::getline(metadata_ss, temp, ':'); // Skipping "# NUMBER OF ZONES:"    
          std::getline(metadata_ss, temp, ','); number_of_zones_ = std::stoi(temp);
          std::getline(metadata_ss, temp, ':'); // Skipping "NUMBER OF NODES:"
          std::getline(metadata_ss, temp, ','); number_of_nodes_ = std::stoi(temp);
          std::getline(metadata_ss, temp, ':'); // Skipping "NUMBER OF LINKS:"
          std::getline(metadata_ss, temp); number_of_links_ = std::stoi(temp);
        }

        links_.clear();
        std::getline(file, line);
        // Reading links information
        while (std::getline(file, line)) {
          std::stringstream ss(line);
          int init, term, type;
          long double capacity, length, free_flow_time, b, power, speed, toll;
          std::string temp;

          std::getline(ss, temp, ','); init = std::stoi(temp) - 1;
          std::getline(ss, temp, ','); term = std::stoi(temp) - 1;
          std::getline(ss, temp, ','); capacity = std::stod(temp);
          std::getline(ss, temp, ','); length = std::stod(temp);
          std::getline(ss, temp, ','); free_flow_time = std::stod(temp);
          std::getline(ss, temp, ','); b = std::stod(temp);
          std::getline(ss, temp, ','); power = std::stod(temp);
          std::getline(ss, temp, ','); speed = std::stod(temp);
          std::getline(ss, temp, ','); toll = std::stod(temp);
          std::getline(ss, temp); type = std::stoi(temp);

          links_.emplace_back(init, term, capacity, length, free_flow_time, b, power, speed, toll, type);
        }

        file.close();
        return true;
      }
      else {
        std::cerr << "Unable to open file: " << net_file_path << std::endl;
        return false;
      }
    }

    bool LoadTrips(const std::filesystem::path& trips_file_path) {
      std::ifstream file(trips_file_path);
      std::string line;

      if (file.is_open()) {
        // Reading metadata
        std::getline(file, line);
        std::stringstream metadata_ss(line);
        std::string temp;

        if (line[0] == '#') {
          std::getline(metadata_ss, temp, ':'); // Skipping "# NUMBER OF ZONES:"
          std::getline(metadata_ss, temp); number_of_zones_ = std::stoi(temp);
        }

        trips_.clear();

        trips_.resize(number_of_zones_);
        for (int i = 0; i < number_of_zones_; i++) {
          trips_[i].resize(number_of_zones_, 0);

          std::getline(file, line);
          std::stringstream ss(line);
          std::string temp;
          long double demand;
          for (int j = 0; j < number_of_zones_ - 1; j++) {
            std::getline(ss, temp, ','); demand = std::stod(temp);
            trips_[i][j] = demand;
          }
          std::getline(ss, temp); demand = std::stod(temp);
          trips_[i][number_of_zones_ - 1] = demand;
        }

        file.close();
        return true;
      }
      else {
        std::cerr << "Unable to open file: " << trips_file_path << std::endl;
        return false;
      }

      trips_.clear();
      return true;
    }

  };
}
#endif
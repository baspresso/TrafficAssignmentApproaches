#include "DataProcessor.h"
#include <fstream>
#include <sstream>
#include <iostream>

namespace TrafficAssignment {

  bool DataProcessor::LoadNet(const std::string& filename) {
    std::ifstream file(filename);
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

      // Reading links information
      while (std::getline(file, line)) {
        std::stringstream ss(line);
        int init, term, type;
        long double capacity, length, free_flow_time, b, power, speed, toll;
        std::string temp;

        std::getline(ss, temp, ','); init = std::stoi(temp);
        std::getline(ss, temp, ','); term = std::stoi(temp);
        std::getline(ss, temp, ','); capacity = std::stoi(temp);
        std::getline(ss, temp, ','); free_flow_time = std::stoi(temp);
        std::getline(ss, temp, ','); b = std::stoi(temp);
        std::getline(ss, temp, ','); speed = std::stoi(temp);
        std::getline(ss, temp, ','); toll = std::stoi(temp);
        std::getline(ss, temp); type = std::stoi(temp);

        links_.emplace_back(init, term, capacity, length, free_flow_time, b, power, speed, toll, type);
      }

      file.close();
      return true;
    }
    else {
      std::cerr << "Unable to open file: " << filename << std::endl;
      return false;
    }
  }

  bool DataProcessor::LoadTrips(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;

    if (file.is_open()) {
      // Reading metadata
      std::getline(file, line);
      std::stringstream metadata_ss(line);

      int vertices, links;
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
          std::getline(ss, temp, ','); demand = std::stoi(temp);
          trips_[i][j] = demand;
        }
        std::getline(ss, temp); demand = std::stoi(temp);
        trips_[i][number_of_zones_ - 1] = demand;
      }

      file.close();
      return true;
    }
    else {
      std::cerr << "Unable to open file: " << filename << std::endl;
      return false;
    }

    trips_.clear();
  }

  bool DataProcessor::LoadData(const std::string& dataset_name) {
    LoadNet(dataset_name + "_net.csv");
    LoadTrips(dataset_name + "_trips.csv");
    return true;
  }
}
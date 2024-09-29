#include "DataProcessor.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <filesystem>

namespace TrafficAssignment {

  bool DataProcessor::LoadNet(const std::filesystem::path& net_file_path) {
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
      //std::cout << number_of_zones_ << ' ' << number_of_nodes_ << ' ' << number_of_links_ << '\n';

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
        std::getline(ss, temp, ','); capacity = std::stoi(temp);
        std::getline(ss, temp, ','); free_flow_time = std::stoi(temp);
        std::getline(ss, temp, ','); b = std::stoi(temp);
        std::getline(ss, temp, ','); speed = std::stoi(temp);
        std::getline(ss, temp, ','); toll = std::stoi(temp);
        std::getline(ss, temp); type = std::stoi(temp);

        links_.emplace_back(init, term, capacity, length, free_flow_time, b, power, speed, toll, type);

        //std::cout << init << ' ' << term << ' ' << capacity << ' ' << length << ' ' << free_flow_time << ' ' << b << ' ' << power << ' ' << speed << ' ' << toll << ' ' << type << '\n';
      }

      file.close();
      return true;
    }
    else {
      std::cerr << "Unable to open file: " << net_file_path << std::endl;
      return false;
    }
  }

  bool DataProcessor::LoadTrips(const std::filesystem::path& trips_file_path) {
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
      std::cerr << "Unable to open file: " << trips_file_path << std::endl;
      return false;
    }

    trips_.clear();
  }

  bool DataProcessor::LoadData(const std::string& dataset_name) {
    std::filesystem::path current_path = std::filesystem::current_path();

    while (current_path.filename() != "out") {
      current_path = current_path.parent_path();
    }
    current_path = current_path.parent_path() / "TrafficAssignmentApproaches";
    auto net_file_path = current_path / "data/" / dataset_name / (dataset_name + "_net.csv");
    auto trips_file_path = current_path / "data/" / dataset_name / (dataset_name + "_trips.csv");

    LoadNet(net_file_path);
    LoadTrips(trips_file_path);
    return true;
  }
}
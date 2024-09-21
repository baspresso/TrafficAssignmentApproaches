#ifndef DATA_PROCESSOR_H
#define DATA_PROCESSOR_H

#include <vector>
#include <string>
#include "Link.h"

namespace TrafficAssignment {
  class DataProcessor {
  public:
    DataProcessor() = default;

    bool LoadData(const std::string& dataset_name);

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

    bool LoadNet(const std::string& filename);

    bool LoadTrips(const std::string& filename);

  };
}
#endif DATA_PROCESSOR_H
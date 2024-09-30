#ifndef STATISTICS_RECORDER_H
#define STATISTICS_RECORDER_H

#include <chrono>
#include <filesystem>
#include <string>
#include <fstream>

namespace TrafficAssignment {

  template <typename T>
  class TrafficAssignmentApproach;

  template <typename T>
  class StatisticsRecorder {
  public:
    StatisticsRecorder(): 
      traffic_assignment_approach_(nullptr), time_on_statistics_(0) {
      recording_start_ = std::chrono::high_resolution_clock::now();
    }

    ~StatisticsRecorder() = default;

    void StartRecording(TrafficAssignmentApproach<T>* traffic_assignment_approach, std::string dataset_name) {
      file_path = std::filesystem::current_path();
      while (file_path.filename() != "out") {
        file_path = file_path.parent_path();
      }
      file_path = file_path.parent_path() / "TrafficAssignmentApproaches" / "performance_results" / dataset_name;
      std::filesystem::create_directories(file_path);
      file_path /= (traffic_assignment_approach->GetApproachName() + ".csv");
      file.open(file_path, std::ios::out);

      file << "RGAP,Objective Function,Time\n";

      traffic_assignment_approach_ = traffic_assignment_approach;
      time_on_statistics_ = 0;

      recording_start_ = std::chrono::high_resolution_clock::now();
      file.close();
    }

    void RecordStatistics() {
      auto statistics_start = std::chrono::high_resolution_clock::now();
      file.open(file_path, std::ios::app);
      T rgap = traffic_assignment_approach_->RelativeGap();
      T objective_function = traffic_assignment_approach_->ObjectiveFunction();
      file << rgap << "," << objective_function << ",";
      time_on_statistics_ += 0.001 * (std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - statistics_start)).count();
      file << std::setprecision(10) << 0.001 * ((std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - recording_start_)).count()) - time_on_statistics_ << '\n';
      file.close();
    }

  private:
    std::chrono::time_point<std::chrono::high_resolution_clock> recording_start_;
    long double time_on_statistics_;
    TrafficAssignmentApproach<T>* traffic_assignment_approach_;
    std::ofstream file;
    std::filesystem::path file_path;
  };
}

#endif STATISTICS_RECORDER_H
#ifndef STATISTICS_RECORDER_H
#define STATISTICS_RECORDER_H

#include <chrono>
#include <string>
#include <fstream>
#include <iomanip>
#include <memory>
#include "../core/Network.h"

namespace TrafficAssignment {

/**
 * @brief Records TAP convergence metrics (RGAP, objective, elapsed time) to CSV.
 *
 * Tracks wall-clock time (excluding time spent computing statistics) and writes
 * periodic snapshots of convergence state. Used by TrafficAssignmentApproach
 * implementations to produce performance traces for analysis.
 *
 * Output CSV columns: RelativeGap, ObjectiveFunction, ElapsedTime(s).
 */
template <typename T>
class StatisticsRecorder {
public:
    StatisticsRecorder(Network<T>& network) :
        network_(network),
        time_on_statistics_(0.0),
        is_recording_(false)
    {
        recording_start_ = std::chrono::high_resolution_clock::now();
    }

    void SetOutputRoot(const std::string& root) { output_root_ = root; }

    /// @brief Initializes output CSV file and starts the wall-clock timer.
    void StartRecording(std::string approach_name) {
        if(is_recording_) return;

        dataset_name_ = network_.name();
        approach_name_ = approach_name;
        InitializeOutputFile();
        ResetTimers();
        is_recording_ = true;
    }

    /// @brief Appends one row (RGAP, objective, elapsed time) to the trace CSV.
    void RecordStatistics() {
        if(!is_recording_) return;
        
        PauseRecording();

        auto statistics_start = std::chrono::high_resolution_clock::now();
        
        file.open(file_path_, std::ios::app);
        const T rgap = network_.RelativeGap();
        const T objective = network_.ObjectiveFunction();
        const T elapsed = GetElapsedTime();
        
        file << std::setprecision(10)
             << rgap << ","
             << objective << ","
             << elapsed << "\n";
        
        file.close();

        ResumeRecording();
    }

    void StopRecording() {
        if(!is_recording_) return;
        
        file.close();
        is_recording_ = false;
        // Note: Removing 'network_ = nullptr;' as it's a reference and cannot be assigned nullptr.
    }

    /// @brief Pauses the wall-clock timer (to exclude statistics computation time).
    void PauseRecording() {
        if(!is_recording_) return;
        pause_start_ = std::chrono::high_resolution_clock::now();
    }

    /// @brief Resumes the wall-clock timer, accumulating the paused duration.
    void ResumeRecording() {
        if(!is_recording_) return;
        auto duration = std::chrono::high_resolution_clock::now() - pause_start_;
        time_on_statistics_ += std::chrono::duration<double>(duration).count();
    }

private:
    Network<T>& network_;
    std::chrono::time_point<std::chrono::high_resolution_clock> recording_start_;
    std::chrono::time_point<std::chrono::high_resolution_clock> pause_start_;
    double time_on_statistics_;
    bool is_recording_;
    std::ofstream file;
    std::string file_path_;
    std::string dataset_name_;
    std::string approach_name_;
    std::string output_root_;

    void InitializeOutputFile() {
        // Create results directory
        const std::string base = output_root_.empty() ? "../performance_results" : output_root_;
        const std::string dir_path = base + "/" + dataset_name_;
        std::filesystem::create_directories(dir_path);
        
        // Create unique filename
        const auto now = std::chrono::system_clock::now();
        const auto timestamp = std::chrono::duration_cast<std::chrono::milliseconds>(
            now.time_since_epoch()).count();
        
        file_path_ = dir_path + "/" + approach_name_ + ".csv";
        
        // Write CSV header
        file.open(file_path_, std::ios::out);
        file << "RelativeGap,ObjectiveFunction,ElapsedTime(s)\n";
        file.close();
    }

    void ResetTimers() {
        time_on_statistics_ = 0.0;
        recording_start_ = std::chrono::high_resolution_clock::now();
    }

    T GetElapsedTime() const {
        auto elapsed = pause_start_ - recording_start_;
        double elapsed_seconds = std::chrono::duration<double>(elapsed).count();
        return static_cast<T>(elapsed_seconds - time_on_statistics_);
    }
};

} // namespace TrafficAssignment

#endif
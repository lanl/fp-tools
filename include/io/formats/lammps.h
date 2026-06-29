#pragma once

#include "base.h"
#include "traj_prov.h"

#include <array>
#include <fstream>
#include <string>

class LAMMPSTrajectory : public TrajectoryProvider {
public:
    explicit LAMMPSTrajectory(const std::string& file_path);
    bool parse_snapshot(Frame& frame) override;

    // Non-copyable: owns an open file handle.
    LAMMPSTrajectory(const LAMMPSTrajectory&)            = delete;
    LAMMPSTrajectory& operator=(const LAMMPSTrajectory&) = delete;

    // Movable: allows for transfer of ownership if needed.
    LAMMPSTrajectory(LAMMPSTrajectory&&)            = default;
    LAMMPSTrajectory& operator=(LAMMPSTrajectory&&) = default;

    ~LAMMPSTrajectory() = default;

private:
    std::ifstream file_;

    // Read buffer for the file stream. Increasing this from the default
    // (~8KB) to 4MB significantly reduces syscall overhead on large files.
    static constexpr size_t BUFFER_SIZE = 4 * 1024 * 1024;
    std::array<char, BUFFER_SIZE> read_buffer_;
};

#pragma once

#include "base.h"
#include "traj_prov.h"

#include <array>
#include <fstream>
#include <string>

class BinaryTrajectory : public TrajectoryProvider {
public:
    explicit BinaryTrajectory(const std::string& file_path);
    bool parse_snapshot(Frame& frame) override;

    // Non-copyable: owns an open file handle.
    BinaryTrajectory(const BinaryTrajectory&)            = delete;
    BinaryTrajectory& operator=(const BinaryTrajectory&) = delete;

    // Movable: allows for transfer of ownership if needed.
    BinaryTrajectory(BinaryTrajectory&&)            = default;
    BinaryTrajectory& operator=(BinaryTrajectory&&) = default;

    ~BinaryTrajectory() = default;

private:
    std::ifstream file_;

    int natoms_;

    // Read buffer for the file stream. Increasing this from the default
    // (~8KB) to 4MB significantly reduces syscall overhead on large files.
    static constexpr size_t BUFFER_SIZE = 4 * 1024 * 1024;
    std::array<char, BUFFER_SIZE> read_buffer_;
};

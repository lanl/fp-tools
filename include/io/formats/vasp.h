#pragma once

#include "base.h"
#include "traj_prov.h"

#include <array>
#include <fstream>
#include <string>

class VASPTrajectory : public TrajectoryProvider {
public:
    explicit VASPTrajectory(const std::string& file_path);

    bool parse_snapshot(Frame& frame) override;

    // Non-copyable: owns an open file handle
    VASPTrajectory(const VASPTrajectory&)            = delete;
    VASPTrajectory& operator=(const VASPTrajectory&) = delete;

    // Movable: allows for transfer of ownership if needed
    VASPTrajectory(VASPTrajectory&&)            = default;
    VASPTrajectory& operator=(VASPTrajectory&&) = default;

    ~VASPTrajectory() = default;

private:
    std::ifstream file_;
    int snap_counter_{};
    int total_atoms_{};

    // Box geometry
    double lx_{}, ly_{}, lz_{};
    double box_bounds_[3][2]{};

    // TODO: find a better way to do this than creating a list of atom_num
    // length that contains strings of all of the atom types
    std::vector<std::string> atom_types_;

    void parse_header();
};

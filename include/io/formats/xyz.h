#pragma once

#include "base.h"
#include "traj_prov.h"

#include <fstream>
#include <string>

class XYZTrajectory : public TrajectoryProvider {
public:
    explicit XYZTrajectory(const std::string& filename);

    // Override the pure virtual method from base class
    bool parse_snapshot(Frame& frame) override;

private:
    std::ifstream file{};
    int timestep_counter{0};
};


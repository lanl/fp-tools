#pragma once

#include "base.h"

class TrajectoryProvider {
public:
    virtual ~TrajectoryProvider() = default;

    // Pure virtual method to be implemented by derived classes
    virtual bool parse_snapshot(Frame& frame) = 0;
};

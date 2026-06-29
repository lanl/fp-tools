#pragma once

#include "base.h"

#include <array>
#include <cstddef>
#include <vector>

// Vector Operations
inline double squared_displacement(const Atom& a, const Atom& b) {
  constexpr double scale = 1;
  double dx = (a.x - b.x) / scale;
  double dy = (a.y - b.y) / scale;
  double dz = (a.z - b.z) / scale;
  return dx*dx + dy*dy + dz*dz;
}

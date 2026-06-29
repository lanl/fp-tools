#pragma once

#include <string>
#include <vector>

// Define a data structure for each atom
struct Atom {
    int id{};
    std::string type{};
    std::string species{};
    
    // Positions; these should ALWAYS be populated.
    double x{}, y{}, z{};

    // Velocities; these should be zero if not present in the trajectory file.
    double vx{}, vy{}, vz{};

    // Forces; these should be zero if not present in the trajectory file.
    double fx{}, fy{}, fz{};

    // Molecule ID; -1 if not present in the trajectory file.
    int molecule_id = -1;
};

// Define a data structure for each trajectory frame (variable names are based on the format of the dump files)
struct Frame {
    int timestep{};
    int atomnum{};

    double box_bounds[3][2]{};
    double box_tilt[3]{};
    
    bool is_triclinic{false};

    std::vector<Atom> atoms{};
};

struct Vector3D {
    double x, y, z;
};


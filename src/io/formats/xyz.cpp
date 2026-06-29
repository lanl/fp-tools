//  TODO: Finish documentation
//  TODO: Make sure the parsing of types is correct
//  TODO: Make sure all of the parsed values match base.h types

/* 
* -----------------------------------------------------------------------------
* XYZ Trajectory Parser
* ----------------------
* This file implements the XYZTrajectory class for reading VASP XDACAR-style
* trajectory files. It supports only rectangular (orthorhombic) boxes and
* assumes atomic positions are provided in fractional ("Direct") format.
*
* The parser reads the box geometry, element types, atom counts, and then
* returns one Frame at a time, scaled into real coordinates.
*
* Authors: leahghartman
* -----------------------------------------------------------------------------
*/

#include "xyz.h"

#include <sstream>

XYZTrajectory::XYZTrajectory(const std::string& filename) {
    file.open(filename);
    if (!file.is_open()) {
      throw std::runtime_error("Could not open XYZ file: " + filename);
    }
}

bool XYZTrajectory::parse_snapshot(Frame& frame) {
    std::string line{};
    int natoms{};

    // Read number of atoms
    if (!std::getline(file, line)) return false;  // End of file
    std::istringstream(line) >> natoms;

    // Skip comment line
    if (!std::getline(file, line)) return false;

    frame.atoms.clear();
    frame.atoms.reserve(natoms);
    frame.timestep = timestep_counter++;
    frame.atomnum = natoms;

    // Parse atoms
    for (int i{0}; i < natoms; ++i) {
        if (!std::getline(file, line)) {
            throw std::runtime_error("Unexpected end of file while reading atom coordinates.");
        }

        std::istringstream iss(line);
        std::string element;
        double x, y, z;

        if (!(iss >> element >> x >> y >> z)) {
            throw std::runtime_error("Invalid atom line in XYZ file");
        }

        Atom atom;
        atom.type = "";  // Optional: convert `element` to int if needed
        atom.x = x;
        atom.y = y;
        atom.z = z;

        // Extended XYZ check: Attempt to read vx, vy, vz from the same stream line
        double vx{0.0}, vy{0.0}, vz{0.0};
        if (iss >> vx >> vy >> vz) {
            atom.vx = vx;
            atom.vy = vy;
            atom.vz = vz;
        } else {
            // Default to 0.0 if the file doesn't include velocities
            atom.vx = 0.0;
            atom.vy = 0.0;
            atom.vz = 0.0;
        }

        frame.atoms.push_back(atom);
    }
    return true;
}

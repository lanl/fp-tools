#include "base.h"
#include "binary.h"

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <cstdint>
#include <cmath>
#include <cstring>

BinaryTrajectory::BinaryTrajectory(const std::string& file_path) {
    file_.open(file_path, std::ios::binary);
    if (!file_.is_open()) {
        throw std::runtime_error("Could not open binary trajectory file: " + file_path);
    }
    file_.rdbuf()->pubsetbuf(read_buffer_.data(), read_buffer_.size());

    // READ NATOMS EXACTLY ONCE HERE
    int32_t m1 = 0, m2 = 0, natoms = 0;
    file_.read(reinterpret_cast<char*>(&m1), 4);
    file_.read(reinterpret_cast<char*>(&natoms), 4);
    file_.read(reinterpret_cast<char*>(&m2), 4);

    if (!file_ || m1 != m2 || m1 != 4) {
        throw std::runtime_error("Failed to read valid Fortran natoms header from file start.");
    }

    // Save this to a member variable in your class header (e.g., int natoms_;)
    this->natoms_ = natoms; 
}

bool BinaryTrajectory::parse_snapshot(Frame& frame)
{
    if (!file_)
        return false;

    // Read the opening marker of the first atom to check for EOF
    int32_t r1 = 0;
    file_.read(reinterpret_cast<char*>(&r1), sizeof(int32_t));
    
    if (file_.gcount() == 0 || file_.eof()) {
        return false; 
    }

    int32_t natoms = this->natoms_; 
    frame.atomnum = natoms;
    
    frame.atoms.clear(); 
    frame.atoms.resize(natoms); 

    struct AtomRecord {
        double x, y, z;
        double type;
        double id;
    };

    constexpr int32_t EXPECTED_REC_SIZE = 40;

    // Process the first atom (r1 is already read)
    if (r1 != EXPECTED_REC_SIZE) {
        throw std::runtime_error("Bad atom record size at frame start. Expected 40 but got " + std::to_string(r1));
    }
    
    AtomRecord data;
    int32_t r2 = 0;
    
    file_.read(reinterpret_cast<char*>(&data), sizeof(AtomRecord));
    file_.read(reinterpret_cast<char*>(&r2), sizeof(int32_t));
    if (!file_ || r1 != r2) throw std::runtime_error("Atom marker mismatch on first atom.");

    // Directly assign to the pre-allocated vector slot to avoid push_back copies
    frame.atoms[0].x = data.x;
    frame.atoms[0].y = data.y;
    frame.atoms[0].z = data.z;
    // OPTIMIZATION: If possible, avoid std::to_string here. 
    // If Atom::type must be a string, this line is your bottleneck:
    frame.atoms[0].type = std::to_string(static_cast<int>(data.type));
    frame.atoms[0].id = static_cast<int>(data.id);

    // Loop through the remaining atoms
    for (int i = 1; i < natoms; ++i) {
        file_.read(reinterpret_cast<char*>(&r1), sizeof(int32_t));
        if (r1 != EXPECTED_REC_SIZE) {
            throw std::runtime_error("Bad atom record size at index " + std::to_string(i));
        }

        file_.read(reinterpret_cast<char*>(&data), sizeof(AtomRecord));
        file_.read(reinterpret_cast<char*>(&r2), sizeof(int32_t));

        if (r1 != r2) {
            throw std::runtime_error("Atom record marker mismatch at index " + std::to_string(i));
        }

        frame.atoms[i].x = data.x;
        frame.atoms[i].y = data.y;
        frame.atoms[i].z = data.z;
        frame.atoms[i].type = std::to_string(static_cast<int>(data.type)); 
        frame.atoms[i].id = static_cast<int>(data.id);
    }

    return true;
}

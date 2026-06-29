#include "vasp.h"
#include "base.h"

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <unordered_map>

VASPTrajectory::VASPTrajectory(const std::string& file_path) {
    
    // Open the file and throw immediately if this doesn't succeed
    file_.open(file_path);
    if (!file_.is_open()) {
        throw std::runtime_error(
            "Could not open VASP trajectory file: " + file_path
        );
    }

    parse_header();
}

void VASPTrajectory::parse_header() {
    std::string line{};

    // Helper lambda to get the next line that contains actual data
    auto get_next_valid_line = [this](std::string& out_line) -> bool {
        while (std::getline(this->file_, out_line)) {
            // Trim any trailing carriage returns, newlines, or spaces
            size_t last = out_line.find_last_not_of(" \t\r\n");
            if (last != std::string::npos) {
                out_line = out_line.substr(0, last + 1);
                
                // If it also has leading spaces, check if it's completely whitespace
                if (out_line.find_first_not_of(" \t\r\n") != std::string::npos) {
                    return true; // Found a line with actual content
                }
            }
        }
        return false; // Reached End-Of-File without finding data
    };

    // -- SYSTEM NAME ---------------------------
    // Safely skips any initial blank lines until it finds the system name (e.g., "Li")
    if (!get_next_valid_line(line)) {
        throw std::runtime_error("XDATCAR: unexpected end of file in header");
    }

    // -- SCALING FACTOR ------------------------
    // Grabs the next valid line containing data
    if (!get_next_valid_line(line)) {
        throw std::runtime_error("XDATCAR: missing scaling factor");
    }

    char* end;
    const double scale = std::strtod(line.c_str(), &end);
    if (end == line.c_str() || scale <= 0) {
        throw std::runtime_error("XDATCAR: invalid scaling factor: '" + line + "'");
    }

    // -- LATTICE VECTORS -----------------------
    double lattice[3][3]{};
    for (int i{0}; i < 3; ++i) {
        if (!get_next_valid_line(line)) {
            throw std::runtime_error(
                "XDATCAR: missing lattice vector " + std::to_string(i + 1));
        }
        std::istringstream iss(line);
        if (!(iss >> lattice[i][0] >> lattice[i][1] >> lattice[i][2])) {
            throw std::runtime_error(
                "XDATCAR: malformed lattice vector: " + line);
        }
    }

    // Check for orthogonality
    const double tol = 1e-6;
    const bool orthogonal = 
        std::abs(lattice[0][1]) < tol && std::abs(lattice[0][2]) < tol &&
        std::abs(lattice[1][0]) < tol && std::abs(lattice[1][2]) < tol &&
        std::abs(lattice[2][0]) < tol && std::abs(lattice[2][1]) < tol;

    if (!orthogonal) {
        throw std::runtime_error(
            "XDATCAR: non-orthogonal simulation cells are not yet supported "
            "in FP-Tools. Lattice vectors must be diagonal."
        );
    }

    lx_ = lattice[0][0] * scale;
    ly_ = lattice[1][1] * scale;
    lz_ = lattice[2][2] * scale;

    box_bounds_[0][0] = 0.0;  box_bounds_[0][1] = lx_;
    box_bounds_[1][0] = 0.0;  box_bounds_[1][1] = ly_;
    box_bounds_[2][0] = 0.0;  box_bounds_[2][1] = lz_;

    // -- SPECIES NAME -----------------------------
    if (!get_next_valid_line(line)) {
        throw std::runtime_error("XDATCAR: missing species names");
    }
    std::vector<std::string> species;
    {
        std::istringstream iss(line);
        std::string s;
        while (iss >> s) species.push_back(s);
    }
    if (species.empty()) {
        throw std::runtime_error("XDATCAR: missing atom counts");
    }

    // -- ATOMS PER SPECIES ------------------------
    if (!get_next_valid_line(line)) {
        throw std::runtime_error("XDATCAR: missing atom counts");
    }
    std::vector<int> counts;
    {
        std::istringstream iss(line);
        int n;
        while (iss >> n) {
            if (n <= 0) {
                throw std::runtime_error("XDATCAR: atom count must be > 0, got " + std::to_string(n));
            }
            counts.push_back(n);
        }
    }
    if (counts.size() != species.size()) {
        throw std::runtime_error(
            "XDATCAR: species name count (" + std::to_string(species.size())
            + ") does not match atom count entries (" 
            + std::to_string(counts.size()) + ")");
    }

    total_atoms_ = 0;
    for (size_t s = 0; s < species.size(); ++s) {
        for (int n{0}; n < counts[s]; ++n) {
            atom_types_.push_back(species[s]);
        }
        total_atoms_ += counts[s];
    }
}

bool VASPTrajectory::parse_snapshot(Frame& frame) {
    std::string line{};

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 2; ++j) {
            frame.box_bounds[i][j] = box_bounds_[i][j];
        }
    }
    frame.atomnum    = total_atoms_; 

    // -- CONFIGURATION HEADER --------------------
    bool found = false;
    bool is_direct = true;

    while (std::getline(file_, line)) {
        if (line.find("configuration=") != std::string::npos) {
            if (line.find("Direct") != std::string::npos 
                || line.find("direct") != std::string::npos) {
                is_direct = true;
            } else if (line.find("Cartesian") != std::string::npos
                       || line.find("cartesian") != std::string::npos) {
                is_direct = false;
            } else {
                throw std::runtime_error(
                    "XDATCAR: unrecognized coordinate type in line: " + line);
            }
            found = true;
            break;
        }
    }
    if (!found) return false;

    // -- SNAPSHOT METADATA ----------------------
    frame.timestep = snap_counter_++;

    // -- ATOM DATA ------------------------------
    frame.atoms.clear();
    frame.atoms.reserve(static_cast<size_t>(total_atoms_));
    
    for (int i{0}; i < total_atoms_; ++i) {
        if (!std::getline(file_, line)) {
            throw std::runtime_error(
                "XDATCAR: unexpected end of file while reading atoms in "
                "snapshot " + std::to_string(frame.timestep)
                + "(expected " + std::to_string(total_atoms_)
                + "atoms, got " + std::to_string(i));
        }

        // Parse three coordinate values from the char buffer
        const char* ptr = line.c_str();
        char* num_end;
        
        const double c1 = std::strtod(ptr, &num_end);
        if (num_end == ptr) {
            throw std::runtime_error(
                "XDATCAR: malformed coordinate line in snapshot "
                + std::to_string(frame.timestep)
                + ", atom " + std::to_string(i));
        }
        ptr = num_end;

        const double c2 = std::strtod(ptr, &num_end);
        if (num_end == ptr) {
            throw std::runtime_error(
                "XDATCAR: malformed coordinate line in snapshot "
                + std::to_string(frame.timestep)
                + ", atom " + std::to_string(i));
        }
        ptr = num_end;

        const double c3 = std::strtod(ptr, &num_end);
        if (num_end == ptr) {
            throw std::runtime_error(
                "XDATCAR: malformed coordinate line in snapshot "
                + std::to_string(frame.timestep)
                + ", atom " + std::to_string(i));
        }

        Atom atom;
        atom.id = i + 1;
        atom.type = atom_types_[i];
        
        if (is_direct) {
            // Convert fractional to Cartesian for an orthogonal box:
            // x_real = frac_x * lx (off-diagonal terms are 0)
            atom.x = c1 * lx_;
            atom.y = c2 * ly_;
            atom.z = c3 * lz_;
        } else {
            // Cartesian coordinates already in real units
            atom.x = c1;
            atom.y = c2;
            atom.z = c3;
        }
        frame.atoms.push_back(std::move(atom));
    }
    return true;
}

#include "base.h"
#include "lammps.h"

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <unordered_map>


LAMMPSTrajectory::LAMMPSTrajectory(const std::string& file_path) {

    // Open the file and throw immediately if not successful.
    file_.open(file_path);
    if (!file_.is_open()) {
        throw std::runtime_error(
            "Could not open LAMMPS trajectory file: " + file_path
        );
    }

    // Increase the read buffer from the default (~8KB) to 4MB. Fewer syscalls
    // should mean much faster file reads.
    file_.rdbuf()->pubsetbuf(read_buffer_.data(), read_buffer_.size());
}

bool LAMMPSTrajectory::parse_snapshot(Frame& frame) {
    std::string line{};

    // Find the start of the next frame. LAMMPS dump frames always begin with
    // "ITEM: TIMESTEP".
    while (std::getline(file_, line)) {
        if (line.find("ITEM: TIMESTEP") != std::string::npos) break;
    }
    if (file_.eof()) return false;

    // -- TIMESTEP ------------------------------
    if (!std::getline(file_, line)) return false;
    frame.timestep = static_cast<int>(std::strtol(line.c_str(), nullptr, 10));

    // -- NUMBER OF ATOMS -----------------------
    if (!std::getline(file_, line)) return false;  // Skip over "ITEM: NUMBER OF ATOMS"
    if (!std::getline(file_, line)) return false;
    frame.atomnum = static_cast<int>(std::strtol(line.c_str(), nullptr, 10));

    // -- BOX BOUNDS ----------------------------
    // Header line format (orthogonal): "ITEM: BOX BOUNDS pp pp pp"
    // Header line format (triclinic):  "ITEM: BOX BOUNDS xy xz yz pp pp pp"
    // ------------------------------------------
    if (!std::getline(file_, line)) return false;

    frame.is_triclinic = (line.find("xy") != std::string::npos);
    
    // Loops once per spatial dimension. First parses the lo values and stores
    // them in box_bounds[i][0], then parses the hi values and stores them in
    // box_bounds[i][1].
    for (int i{0}; i < 3; ++i) {
        if (!std::getline(file_, line)) return false;
        const char* ptr = line.c_str();
        char* end;
        frame.box_bounds[i][0] = std::strtod(ptr, &end); ptr = end;
        frame.box_bounds[i][1] = std::strtod(ptr, &end); ptr = end;

        // If the dump file is in a triclinic geometry, then store its tilt factors.
        if (frame.is_triclinic) {
            frame.box_tilt[i] = std::strtod(ptr, nullptr);
        } else {
            frame.box_tilt[i] = 0.0;
        }
    }
    
    // Parse the box lengths and origins for corrdinate unscaling and unwrapping.
    const double xlo = frame.box_bounds[0][0];
    const double ylo = frame.box_bounds[1][0];
    const double zlo = frame.box_bounds[2][0];
    const double lx  = frame.box_bounds[0][1] - xlo;
    const double ly  = frame.box_bounds[1][1] - ylo;
    const double lz  = frame.box_bounds[2][1] - zlo;
    
    // -- ATOM HEADER ---------------------------
    // Parses "ITEM: ATOMS id type x y z ..." into a name->index map.
    // ------------------------------------------
    if (!std::getline(file_, line)) return false; // "ITEM: ATOMS ..."
    std::unordered_map<std::string, int> col;
    {
        std::istringstream iss(line.substr(12));  // Skip the "ITEM: ATOMS " portion
        std::string name{};
        int idx{0};
        while (iss >> name) col[name] = idx++;
    }

    // -- RESOLVE COLUMN INDICES ----------------
    // An anonymous function defined inline to assist in finding the coordinate
    // type of the output data. [&] means it captures all local variables by
    // reference without them being passed as arguments. "-> int" is the
    // explicit return type.
    // ------------------------------------------
    auto find_coord = [&](const std::string& base) -> int {
        for (const std::string& suffix : {"u", "su", "s", ""}) {
            auto it = col.find(base + suffix);
            if (it != col.end()) return it->second;
        }
        throw std::runtime_error(
            "LAMMPS snapshot at timestep " + std::to_string(frame.timestep)
            + " is missing coordinate column " + base);
    };

    // Detect whether the found coordinates are scaled. Note that here it's
    // being assumed that LAMMPS is using the same coordinate type for all 
    // three dimensions.
    auto coord_is_scaled = [&](const std::string& base) -> bool {
        return col.count(base + "su") > 0 || col.count(base + "s") > 0;
    };

    // Resolve all of the column indices before the actual atom loop. For
    // example, if you had a file with a header "id type xu yu zu", x_idx
    // would be 2, y_idx would be 3, and z_idx would be 4.
    const int x_idx = find_coord("x");
    const int y_idx = find_coord("y");
    const int z_idx = find_coord("z");
    const bool coords_scaled = coord_is_scaled("x");

    // -- ATOM ID -------------------------------
    // col.count() counts how many entries in the map have the key "id". If 1,
    // then its column index will be the value stored for "id" in the map.
    // ------------------------------------------
    const bool has_id = col.count("id") > 0;
    const int  id_idx = has_id ? col["id"] : -1;

    // -- TYPE-----------------------------------
    // Check for two possible column names -- "type" or "element"
    // ------------------------------------------
    int type_idx{-1};
    if      (col.count("type"))    type_idx = col["type"];
    else if (col.count("element")) type_idx = col["element"];

    // -- VELOCITIES (optional) -----------------
    const int vx_idx = col.count("vx") ? col["vx"] : -1;
    const int vy_idx = col.count("vy") ? col["vy"] : -1;
    const int vz_idx = col.count("vz") ? col["vz"] : -1;
    const bool has_velocity = (vx_idx >= 0 && vy_idx >= 0 && vz_idx >= 0);

    // -- FORCES (optional) ---------------------
    const int fx_idx = col.count("fx") ? col["fx"] : -1;
    const int fy_idx = col.count("fy") ? col["fy"] : -1;
    const int fz_idx = col.count("fz") ? col["fz"] : -1;
    const bool has_forces = (fx_idx >= 0 && fy_idx >= 0 && fz_idx >= 0);

    // -- IMAGE FLAGS (optional) ----------------
    // Some users dump "x y z ix iy iz" instead of "xu yu zu". We can
    // reconstruct unwrapped coordinates using this information. This is
    // equivalent to dumping xu/yu/zu directly, but produces smaller files
    // since the integer image flags compress better than floats.
    // ------------------------------------------
    const int ix_idx = col.count("ix") ? col["ix"] : -1;
    const int iy_idx = col.count("iy") ? col["iy"] : -1;
    const int iz_idx = col.count("iz") ? col["iz"] : -1;
    const bool has_image_flags = (ix_idx >= 0 && iy_idx >= 0 && iz_idx >= 0);

    // -- MOLECULE ID (optional) ----------------
    const int mol_idx = col.count("mol") ? col["mol"] : -1;

    // -- ATOM DATA -----------------------------
    // Lines are tokenized by ...
    // ------------------------------------------

    // Removes atoms from the previous snapshot without releasing the vector's
    // allocated memory. Then, reserve() preallocates exactly the right space
    // so no reallocations happen during the atom loop.
    frame.atoms.clear();
    frame.atoms.reserve(static_cast<size_t>(frame.atomnum));

    // The tok vector is reused across lines to avoid repeated allocation.
    std::vector<char*> tok;
    tok.reserve(16);

    for (int i{0}; i < frame.atomnum; ++i) {
        // A file ending mid-snapshot is a malformed file, so throw an error.
        if (!std::getline(file_, line)) {
            throw std::runtime_error(
                "Unexpected end of file while reading atoms at timestep "
                + std::to_string(frame.timestep));
        }

        // The tokenizer walks through the char buffer. *ptr dereferences the
        // pointer and reads the character at the current position.
        tok.clear();
        char* ptr = const_cast<char*>(line.c_str());
        while (*ptr) {
            // Skips white space
            while (*ptr == ' ' || *ptr == '\t') ++ptr;
            if (!*ptr) break;

            // Records a pointer to the start of the current token
            tok.push_back(ptr);

            // Advances past the token to the next whitespace
            while (*ptr != ' ' && *ptr != '\t' && *ptr != '\0') ++ptr;
            if (*ptr) *ptr++ = '\0';
        }

        Atom atom;

        // Uses the already-found id_idx to index directly into the token 
        // pointer array. If no "id" column exists, fall back to the sequential
        // index (i + 1).
        atom.id = has_id
            ? static_cast<int>(std::strtol(tok[id_idx], nullptr, 10))
            : (i + 1);

        // Must store the type of the atom as a string as it could be anything
        // from "C" to "1" to "Al", or anything else.
        atom.type = (type_idx >= 0)
            ? std::string(tok[type_idx])
            : "unknown";

        // Parse the raw coordinate values.
        double raw_x = std::strtod(tok[x_idx], nullptr);
        double raw_y = std::strtod(tok[y_idx], nullptr);
        double raw_z = std::strtod(tok[z_idx], nullptr);

        if (coords_scaled) {
            // Converts fractional coordinates to real Cartesian coordinates.
            atom.x = raw_x * lx + xlo;
            atom.y = raw_y * ly + ylo;
            atom.z = raw_z * lz + zlo;

            // Handle possible floating point overshoot at periodic boundaries
            if (atom.x < xlo)      atom.x += lx;
            if (atom.y < ylo)      atom.y += ly;
            if (atom.z < zlo)      atom.z += lz;
            if (atom.x > xlo + lx) atom.x -= lx;
            if (atom.y > ylo + ly) atom.y -= ly;
            if (atom.z > zlo + lz) atom.z -= lz;

        } else {

            // If the coordinates are unwrapped, then the actual positions
            // are just the raw values we already read above.
            atom.x = raw_x;
            atom.y = raw_y;
            atom.z = raw_z;

            // If image flags are present, we can reconstruct unwrapped
            // positions using x_unwrapped = x_wrapped + ix * lx. 
            if (has_image_flags) {
                const int ix = static_cast<int>(
                        std::strtol(tok[ix_idx], nullptr, 10));
                const int iy = static_cast<int>(
                        std::strtol(tok[iy_idx], nullptr, 10));
                const int iz = static_cast<int>(
                        std::strtol(tok[iz_idx], nullptr, 10));

                atom.x += ix * lx;
                atom.y += iy * ly;
                atom.z += iz * lz;
            }
        }

        // If the dump file contains velocities, then parse and store them.
        if (has_velocity) {
            atom.vx = std::strtod(tok[vx_idx], nullptr);
            atom.vy = std::strtod(tok[vy_idx], nullptr);
            atom.vz = std::strtod(tok[vz_idx], nullptr);
        } else {
            atom.vx = atom.vy = atom.vz = 0.0;
        }

        // If the dump file contains forces, then parse and store them.
        if (has_forces) {
            atom.fx = std::strtod(tok[fx_idx], nullptr);
            atom.fy = std::strtod(tok[fy_idx], nullptr);
            atom.fz = std::strtod(tok[fz_idx], nullptr);
        } else {
            atom.fx = atom.fy = atom.fz = 0.0;
        }

        // If the dump file contains a molecule ID, then parse and store it.
        atom.molecule_id = (mol_idx >= 0)
            ? static_cast<int>(std::strtol(tok[mol_idx], nullptr, 10))
            : -1;

        frame.atoms.push_back(std::move(atom));
    }
    return true;
}

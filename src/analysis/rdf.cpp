#include "rdf.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <algorithm>

// =======================================================
// RDF Constructor
// =======================================================

RDF::RDF(double dr, double r_max)
    : dr_{dr}, r_max_{r_max}, 
      num_bins_{static_cast<size_t>(std::ceil(r_max / dr))} {

    // r_vals_ stores bin centers
    r_vals_.resize(num_bins_);
    for (size_t i = 0; i < num_bins_; ++i) {
        r_vals_[i] = (i + 0.5) * dr_;
    }
}

// =======================================================
// RDF::calculate
// =======================================================

void RDF::calculate(const std::vector<Frame>& frames,
                    const std::vector<std::vector<std::string>>& pairs_input,
                    Logger *logger) {

    // Saves the input pairs_ to the internal pairs_ vector so other methods can
    this->pairs_ = pairs_input;

    // Initializes the histogram; one row per pair, one column per bin
    histogram_.assign(pairs_.size(), std::vector<double>(num_bins_, 0.0));

    // Keeps track of total and completed frames to show progress with pacman
    size_t total = frames.size();
    size_t processed = 0;

    // Iterate through each frame
    for (const auto &fr : frames) {
        const auto &atoms = fr.atoms;

        // Extract box lengths from bounds in each direction. These are used for
        // PBCs
        double xlo = fr.box_bounds[0][0], xhi = fr.box_bounds[0][1];
        double ylo = fr.box_bounds[1][0], yhi = fr.box_bounds[1][1];
        double zlo = fr.box_bounds[2][0], zhi = fr.box_bounds[2][1];
        double lx = xhi - xlo, ly = yhi - ylo, lz = zhi - zlo;

        // Organize atoms by their type for faster lookup during pairwise
        // comparisons
        std::unordered_map<std::string, std::vector<size_t>> atoms_by_type;
        for (size_t i = 0; i < atoms.size(); ++i) {
            atoms_by_type[atoms[i].type].push_back(i);
        }

        // Loop over each pair (like H-H or O-H, etc.)
        for (size_t p = 0; p < pairs_.size(); ++p) {

            // Extract atom type strings for the current pair
            const auto &pair = pairs_[p];

            if (pair.size() == 1 && pair[0] == "all") {

                // Loop over every atom in the box without filtering types
                for (size_t i{0}; i < atoms.size(); ++i) {
                    const Atom& A = atoms[i];

                    // Loop over atoms AFTER atom A to avoid self- and double-counting
                    for (size_t j = i + 1; j < atoms.size(); ++j) {
                        const Atom& B = atoms[j];

                        // Compute the displacement between atoms using the MIC (for PBCs)
                        double dx = A.x - B.x;
                        dx -= lx * std::nearbyint(dx / lx);
                        double dy = A.y - B.y;
                        dy -= ly * std::nearbyint(dy / ly);
                        double dz = A.z - B.z;
                        dz -= lz * std::nearbyint(dz / lz);
                        double r2 = dx * dx + dy * dy + dz * dz;

                        // Ignore pairs_ beyond the maximum RDF radius
                        if (r2 >= r_max_ * r_max_)
                            continue;

                        // Convert r to a bin index
                        int bin = static_cast<int>(std::sqrt(r2) / dr_);

                        if (std::sqrt(r2) > std::min<double>({lx, ly, lz}) / 2) continue;

                        // If valid, increment the histogram by 2.0 to account for both i-j
                        // and j-i
                        if (bin >= 0 && bin < (int)num_bins_) {
                            histogram_[p][bin] += 2.0;
                        }
                    }
                }
                continue;
            }

            const std::string &type1 = pair[0];
            const std::string &type2 = pair[1];

            // Find the indices of the atoms with the specified types
            auto it1 = atoms_by_type.find(type1);
            auto it2 = atoms_by_type.find(type2);
            if (it1 == atoms_by_type.end() || it2 == atoms_by_type.end()) {
                continue; // no atoms of one of the types
            }

            // 'second' contains all the atoms of type one or two depending on the
            // statement
            const auto &group1 = it1->second;
            const auto &group2 = it2->second;

            // Loop over all unique pairs_ to avoid double counting
            if (type1 == type2) {
                // Fetch atoms of the same type
                for (size_t i = 0; i < group1.size(); ++i) {
                    const Atom &A = atoms[group1[i]];

                    // Loops over atoms AFTER atom A to avoid double and self-counting
                    for (size_t j = i + 1; j < group1.size(); ++j) {
                        const Atom &B = atoms[group1[j]];

                        // Compute the displacement between atoms using the MIC (for PBCs)
                        double dx = A.x - B.x;
                        dx -= lx * std::nearbyint(dx / lx);
                        double dy = A.y - B.y;
                        dy -= ly * std::nearbyint(dy / ly);
                        double dz = A.z - B.z;
                        dz -= lz * std::nearbyint(dz / lz);
                        double r2 = dx * dx + dy * dy + dz * dz;

                        // Ignore pairs_ beyond the maximum RDF radius
                        if (r2 >= r_max_ * r_max_)
                            continue;

                        // Convert r to a bin index
                        int bin = static_cast<int>(std::sqrt(r2) / dr_);

                        // If valid, increment the histogram by 2.0 to account for both i-j
                        // and j-i
                        if (bin >= 0 && bin < (int)num_bins_) {
                            histogram_[p][bin] += 2.0;
                        }
                    }
                }
            } else {
                // Loop over all combinations of atoms from type1 and type2
                for (size_t i = 0; i < group1.size(); ++i) {
                    const Atom &A = atoms[group1[i]];

                    // No need to skip i==j since the types differ
                    for (size_t j = 0; j < group2.size(); ++j) {
                        const Atom &B = atoms[group2[j]];

                        // Same distance and binning logic as the above case
                        double dx = A.x - B.x;
                        dx -= lx * std::nearbyint(dx / lx);
                        double dy = A.y - B.y;
                        dy -= ly * std::nearbyint(dy / ly);
                        double dz = A.z - B.z;
                        dz -= lz * std::nearbyint(dz / lz);
                        double r2 = dx * dx + dy * dy + dz * dz;

                        if (r2 >= r_max_ * r_max_)
                            continue;
                
                        int bin = static_cast<int>(std::sqrt(r2) / dr_);

                        // Increment by 1.0 since we're not avoiding duplicates
                        if (bin >= 0 && bin < (int)num_bins_) {
                            histogram_[p][bin] += 1.0;
                        }
                    }
                }
            }
        }
        // Report progress using pacman
        if (logger && ((processed + 1) % 10 == 0 || processed + 1 == total)) {
            print_progress_bar(processed + 1, total, 40, "RDF accumulation");
        }
        ++processed;
    }
}

void RDF::normalize(
    double box_volume, size_t num_frames,
    const std::unordered_map<std::string, size_t> &atoms_per_type_counts) {
    
    if (rdf_.size() != pairs_.size()) {
        rdf_.assign(pairs_.size(), std::vector<double>(num_bins_, 0.0));
    }

    // Total number of particles in the system (regardless of type)
    size_t total_num_atoms = 0;
    for (const auto &kv : atoms_per_type_counts) {
        total_num_atoms += kv.second;
    }
    double N_all = static_cast<double>(total_num_atoms);
    double rho_all = N_all / box_volume;

    double number_density = static_cast<double>(total_num_atoms) / box_volume;

    for (size_t p_idx = 0; p_idx < pairs_.size(); ++p_idx) {
        const auto &pair = pairs_[p_idx];

        double N_i{0.0};
        double rho_j{0.0};

        if (pair.size() == 1 && pair[0] == "all") {
            N_i = N_all;
            rho_j = rho_all;
        } else {
            // Standard individual pair normalization
            const std::string &type_i = pair[0];
            const std::string &type_j = pair[1];
            N_i = static_cast<double>(atoms_per_type_counts.at(type_i));
            double N_j = static_cast<double>(atoms_per_type_counts.at(type_j));
            rho_j = N_j / box_volume;
        }

        for (size_t bin = 0; bin < num_bins_; ++bin) {
            double r_lo = bin * dr_;
            double r_hi = r_lo + dr_;
            double shell_vol = (4.0 / 3.0) * M_PI * (std::pow(r_hi, 3) - std::pow(r_lo, 3));

            double raw = static_cast<double>(histogram_[p_idx][bin]);
            rdf_[p_idx][bin] = raw / (num_frames * N_i * shell_vol * rho_j);
        }
    }
}

void RDF::write(const std::string &output_dir, double& min_box_length) const {
    namespace fs = std::filesystem;
    const auto &pair_list = this->pairs_; // filtered pairs_

    for (std::size_t p = 0; p < pair_list.size(); ++p) {
        std::string filename{};

        if (pair_list[p].size() == 1 && pair_list[p][0] == "all") {
            filename = output_dir + "/rdf-all.dat";
        } else {
            filename = output_dir + "/rdf-" + pair_list[p][0] + "_" + pair_list[p][1] + ".dat";
        }
        std::ofstream out(filename);
        if (!out) {
            std::cerr << "Error: cannot open " << filename << '\n';
            continue;
        }

        if (pair_list[p].size() == 1 && pair_list[p][0] == "all") {
             out << "# r  g(r) for pair " << pair_list[p][0] << '\n';
        } else {
            out << "# r  g(r) for pair " << pair_list[p][0] << "-" << pair_list[p][1] << '\n';
        }

        for (std::size_t i = 0; i < num_bins_; ++i) {
            if (r_vals_[i] >= min_box_length / 2) {
                continue;
            } else {
                out << r_vals_[i] << ' ' << rdf_[p][i] << '\n';
            }
        }
    }
}

// ----------- Plotting -------------

void RDF::plot(
    const std::vector<std::vector<std::string>>& input_file_groups,
    const std::vector<std::vector<std::string>>& label_groups,
    const std::string& output_dir, const std::string& format) const {

    const std::string script_path = std::string(SOURCE_DIR) + "/src/fp_plot.py";

    for (size_t i = 0; i < input_file_groups.size(); ++i) {
        std::string command = "python3 " + script_path + " rdf " + std::to_string(i);

        const auto& files = input_file_groups[i];
        const auto& labels = label_groups[i];

        for (size_t j = 0; j < files.size(); ++j) {
            command += " " + files[j] + " \"" + labels[j] + "\"";
        }

        command += " " + output_dir + " " + format;

        std::cout << "RDF plot command: " << command << std::endl;

        int result = std::system(command.c_str());
        if (result != 0) {
            std::cerr << "Warning: RDF plot command failed for group " << i << "\n";
        }
    }
}

// ----------- Factory functions -------------

RDF build_rdf(const std::vector<Frame> &frames, double dr,
                           double r_max, double box_volume,
                           std::vector<std::vector<std::string>> pairs_,
                           Logger *logger) {
    RDF rdf(dr, r_max);

    std::unordered_map<std::string, size_t> atom_type_counts;
    for (const auto &atom : frames[0].atoms) {
        atom_type_counts[atom.type]++;
    }

    rdf.calculate(frames, pairs_, logger);
    rdf.normalize(box_volume, frames.size(), atom_type_counts);
    return rdf;
}

/**
 * -----------------------------------------------------------------------------
 * @file mfpt.cpp
 * @brief Implementation of Mean First Passage Time (MFPT) analysis tools.
 *
 * This file contains the implementation of:
 *   - FPTAccumulator methods for collecting and finalizing first-passage time
 *     distributions.
 *   - Functions to compute FPT distributions over molecular dynamics
 * trajectories.
 *   - Functions to write MFPT summary data and F_r(t) distributions to files.
 *   - Utilities to compute and plot the local log-derivative D(r) of MFPT data.
 *
 * Authors: Leah Hartman
 * -----------------------------------------------------------------------------
 */

#include "mfpt.h"
#include "logger.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// Include OpenMP only if hte user turns it on
#ifdef USE_OMP
#include <atomic>
#include <omp.h>
#endif

// =============================================================================
// ===                     FPTAccumulator Class Methods                      ===
// =============================================================================

// Constructor: initialize sums and counters for first-passage times
FPTAccumulator::FPTAccumulator(int num_bins)
    : fpt_sum(num_bins, 0.0), counts(num_bins, 0), finalized(false) {}

// Accumulate a first-passage time (FPT) into the specified lag bin
void FPTAccumulator::accumulate(int bin_index, double fpt) {
    if (bin_index >= 0 && static_cast<size_t>(bin_index) < fpt_sum.size()) {
        fpt_sum[bin_index] += fpt;
        counts[bin_index] += 1;
    }
}

std::string get_species_suffix(const std::string& species) {
    if (species.empty() || species =="all") return "_all";
    return "_" + species;
}

// Finalize probability distribution F_r(t) by converting counts to
// probabilities using the total number of first-passage events
// recorded(total_counts).
void FPTAccumulator::finalize(double dt) {
    if (finalized)
        return; // Prevent double-finalization

    fpt_dist.resize(counts.size(), 0.0);
    total_counts = 0;
    for (int c : counts)
        total_counts += c; // Sum over all lag bins

    if (total_counts > 0) {
        for (size_t i{0}; i < counts.size(); ++i) {
            fpt_dist[i] = counts[i] / (total_counts * dt); // Normalize to PDF
        }
    }
    finalized = true;
}

// Accessor methods
const std::vector<double> &FPTAccumulator::get_fpt_dist() const {
    return fpt_dist;
}
int FPTAccumulator::get_total_counts() const { return total_counts; }

const std::vector<double> FPTAccumulator::get_fpt_sum() const {
    return fpt_sum;
}

const std::vector<int> FPTAccumulator::get_counts() const {
    return counts;
}


// =============================================================================
// ===                     Writing FPT and MFPT Outputs                      ===
// =============================================================================

/**
 * @brief Write FPT distributions and summarized MFPT data to files.
 *
 * For each radius bin:
 *   - Writes F_r(t) distributions to individual files (ftp_N.out).
 *   - Computes MFPT statistics t1/t0 and writes to summary file mfpt.dat.
 *
 * @param output_file Path for the MFPT summary file (mfpt.dat).
 * @param mfpt_accums Vector of accumulators containing finalized FPT data.
 * @param dr Radius bin width.
 * @param dt Time step between trajectory frames.
 */

void write_fpt(const std::string &mfpt_file, const std::string &fpt_dir,
               const std::vector<FPTAccumulator> &fpt_accums, double dr,
               double dt, double box_volume, size_t total_atoms,
               const std::string& species) {

    // Generate a unique suffix for this species combination
    std::string suffix = get_species_suffix(species);

    std::ofstream mfpt_out(mfpt_file);
    mfpt_out << "# mean first passage time T(r) vs r for species: " 
             << (species.empty() ? "all" : species) << "\n";
    mfpt_out << "# 1) r/a , 2) MFPT , 3) t1 , 4) t2 , 5) counts\n";

    double density = static_cast<double>(total_atoms) / box_volume;

    for (size_t radius_idx{0}; radius_idx < fpt_accums.size(); ++radius_idx) {
        const auto &accum = fpt_accums[radius_idx];
        const auto &fpt_dist = accum.get_fpt_dist();
        int total_counts = accum.get_total_counts();

        double r = dr * (radius_idx + 1);

        std::string fpt_filename =
            fpt_dir + "/fpt" + suffix + "_" + std::to_string(radius_idx + 1) + ".out";
        std::ofstream fpt_out(fpt_filename);
        fpt_out << "# Probability distribution F_r(t) of first passages at "
                   "distance r = " << r << "\n";
        fpt_out << "# 1) t , 2) F_r(t)\n";
        fpt_out << std::scientific << std::setprecision(15);

        // Accumulate integrals for MFPT calculation
        double t0 = 0.0, t1 = 0.0, t2 = 0.0;
        for (size_t ibin = 0; ibin < fpt_dist.size(); ++ibin) {
            double time = (ibin + 1) * dt;
            double frt = fpt_dist[ibin];
            fpt_out << time << "\t" << frt << "\n";

            t0 += frt * dt;               // Zeroth moment
            t1 += time * frt * dt;        // First moment (mean)
            t2 += time * time * frt * dt; // Second moment
        }
        fpt_out.close();
	
        const auto &fpt_sum = accum.get_fpt_sum();
        const auto &counts = accum.get_counts();

        // Compute MFPT as whatever it actually is
        double mfpt = (total_counts > 0 && t0 > 0) ? t1/t0 : 0.0;

        // Write summary entry to mfpt.dat
        mfpt_out << std::scientific << std::setprecision(15) 
                 << r << "\t" << mfpt << "\t" << t1 << "\t" << t2 << "\t" 
                 << total_counts << "\n";
    }
    mfpt_out.close();
}

// =============================================================================
// ===               Compute FPT Distributions over Trajectory               ===
// =============================================================================

/**
 * @brief Accumulate FPT distributions across the trajectory frames.
 *
 * For each radius bin, accumulates first passage times for all particles,
 * starting from every possible initial frame, and records lag times at which
 * particles first cross the radius threshold.
 *
 * @param frames Trajectory frames to process.
 * @param dr Radius bin width.
 * @param r_max Maximum radius threshold.
 * @param dt Time step between frames.
 * @return Vector of MFPTAccumulator objects with accumulated FPT data.
 */

std::vector<FPTAccumulator> compute_fpt(const std::vector<Frame> &frames,
                                        double dr, double r_max, double dt,
                                        const std::string& species) {
    size_t n_frames = frames.size();

    // -- SPECIES FILTERING -----------------------------------
    std::vector<size_t> atom_indices;
    if (species.empty() || species == "all") {
        // No filter then; include every atom
        atom_indices.resize(frames[0].atoms.size());
        for (size_t i = 0; i < atom_indices.size(); ++i) atom_indices[i] = i;
    } else {
        for (size_t i = 0; i < frames[0].atoms.size(); ++i) {
            if (frames[0].atoms[i].type == species) {
                atom_indices.push_back(i);
            }
        }
    }

    if (atom_indices.empty()) {
        throw std::runtime_error(
            "compute_fpt: no atoms found matching requested species");
    }

    const size_t n_atoms       = atom_indices.size();
    const int    n_bins_radius = static_cast<int>(std::ceil(r_max / dr));
    const int    max_lag_steps = static_cast<int>(n_frames - 1);

    std::vector<double> fpt (n_bins_radius, 0.0);
    std::vector<int> counts (n_bins_radius, 0);

    // Prepare accumulators: one per radius threshold
    std::vector<FPTAccumulator> accumulators (
        n_bins_radius, FPTAccumulator(max_lag_steps + 1));

    // TODO: make this configurable via origin_stride
    // Every 5th snapshot is used as a time origin. I'm trying to trade statistical 
    // accuracy for speed.
    const size_t origin_stride = 5;
    size_t total_origins =
        n_bins_radius * (n_frames / origin_stride + (n_frames % origin_stride != 0));
    size_t done_origins = 0;    // cumulative progress counter
  
    // Loop over radius thresholds (distance epsilon)
    for (int radius_idx = 0; radius_idx < n_bins_radius; ++radius_idx) {
        const double epsilon = dr * (radius_idx + 1);
        const double epsilon2 = epsilon * epsilon;

        // Loop over all possible time origins t0
        for (size_t t0 = 0; t0 < n_frames; t0+=origin_stride) {
            const Frame& frame_start = frames[t0];

            // Record initial positions of all atoms at t0
            std::vector<double> x0(n_atoms), y0(n_atoms), z0(n_atoms);
            for (size_t a = 0; a < n_atoms; ++a) {
                const size_t atom_idx = atom_indices[a];
                x0[a] = frame_start.atoms[atom_idx].x;
                y0[a] = frame_start.atoms[atom_idx].y;
                z0[a] = frame_start.atoms[atom_idx].z;
            }
       
            // For each atom, find first passage time crossing epsilon
            for (size_t a = 0; a < n_atoms; ++a) {
                const size_t atom_idx = atom_indices[a];
	            int first_crossing_lag = -1;

                for (size_t t = t0 + 1; t < n_frames; ++t) {
                    const auto& atom = frames[t].atoms[atom_idx];
                    double dx = atom.x - x0[a];
                    double dy = atom.y - y0[a];
                    double dz = atom.z - z0[a];
                    double dist2 = dx*dx + dy*dy + dz*dz;

                    if (dist2 >= epsilon2) {
                        first_crossing_lag = static_cast<int>(t - t0 - 1);
	                    fpt[radius_idx]    += first_crossing_lag * dt;
	                    counts[radius_idx] += 1;
                        break;
                    }
                }
                if (first_crossing_lag >= 0) {
                    accumulators[radius_idx].accumulate(first_crossing_lag, dt * first_crossing_lag);
                }
            }

            // Update cumulative progress and print progress bar once every 10 or at
            // end
            ++done_origins;
            if (done_origins % 10 == 0 || done_origins == total_origins) {
                print_progress_bar(done_origins, total_origins, 40, "MFPT origins");
            }
        }
    }

    std::cout << std::endl; // End the progress bar line

    // Finalize the probability distribution for each radius
    for (auto& accum : accumulators) accum.finalize(dt);
    return accumulators;
}

// =============================================================================
// ===                         Plotting Utilities                            ===
// =============================================================================

/**
 * @brief Execute Python script to plot MFPT data.
 */
void plot_mfpt(const std::string &summary_file, const std::string &output_dir,
               const std::string &format) {
    namespace fs = std::filesystem;
    fs::path abs_summary_file = fs::absolute(summary_file);
    fs::path abs_output_dir = fs::absolute(output_dir);

    const std::string script_path = std::string(SOURCE_DIR) + "/src/fp_plot.py";
    std::string command = "python3 " + script_path + " mfpt " +
                        abs_summary_file.string() + " " +
                        abs_output_dir.string() + " " + format;

    int result = std::system(command.c_str());
}

// =============================================================================
// ===               Log-Derivative Computation D(r) and Plot                ===
// =============================================================================

/**
 * @brief Compute local log-derivative D(r) from MFPT summary file.
 */
void compute_dofr(const std::string &mfpt_file,
                  const std::string &output_file) {

    std::ifstream infile(mfpt_file);
    std::vector<double> r_vals, t_vals;
    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        double r, t;
        if (!(iss >> r >> t)) continue;
        r_vals.push_back(r);
        t_vals.push_back(t);
    }
    infile.close();

    std::ofstream outfile(output_file);
    outfile << "# Local log-derivative D(r)\n";
    outfile << "# 1) r (midpoint), 2) D(r)\n";

    for (size_t i = 1; i < r_vals.size(); ++i) {
        double r1 = r_vals[i - 1], t1 = t_vals[i - 1];
        double r2 = r_vals[i], t2 = t_vals[i];

        if (r1 <= 0 || r2 <= 0 || t1 <= 0 || t2 <= 0) continue;

        double log_r1 = std::log(r1);
        double log_r2 = std::log(r2);
        if (std::abs(log_r2 - log_r1) < 1e-6) continue;

        double derivative = (std::log(t2) - std::log(t1)) / (log_r2 - log_r1);
        double r_mid = 0.5 * (r1 + r2);
        outfile << r_mid << " " << derivative << "\n";
    }
    outfile.close();
}

/**
 * @brief Execute Python script to plot D(r) data.
 */
void plot_dofr(const std::string &log_derivative_file,
               const std::string &output_dir, const std::string &format) {
    namespace fs = std::filesystem;
    fs::path abs_log_file = fs::absolute(log_derivative_file);
    fs::path abs_output_dir = fs::absolute(output_dir);

    const std::string script_path = std::string(SOURCE_DIR) + "/src/fp_plot.py";
    std::string command = "python3 " + script_path + " dlog " +
                        abs_log_file.string() + " " + abs_output_dir.string() +
                        " " + format;

    int result = std::system(command.c_str());
}

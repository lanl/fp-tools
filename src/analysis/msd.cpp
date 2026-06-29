/*
 * -----------------------------------------------------------------------------
 * Mean Squared Displacement (MSD) Accumulator
 * -------------------------------------------
 * This file implements the MSDAccumulator class for computing Mean Squared
 * Displacement (MSD) from molecular dynamics trajectories. It accumulates
 * squared displacements for different time lags, averages them, writes the
 * results to disk, and can call a Python script to plot the MSD curve.
 *
 * Authors: leahghartman
 * -----------------------------------------------------------------------------
 */

#include "logger.h"
#include "basic_math.h"
#include "msd.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <system_error>
#include <unordered_map>
#include <unordered_set>

#include <numeric>

// Include OpenMP only if the user turns it on
#ifdef USE_OMP
#include <atomic>
#include <omp.h>
#endif

/*
 * Constructor for MSDAccumulator
 *
 * Initializes vectors to store summed squared displacements (msd_sum) and
 * the number of samples contributing to each lag (counts).
 *
 * Parameters:
 *   - max_lag: maximum lag time (number of frames) to calculate MSD over
 */
MSDAccumulator::MSDAccumulator(int max_lag, int num_groups) : finalized(false) {
  msd_sum.resize(num_groups, std::vector<double>(max_lag + 1, 0.0));
  counts.resize(num_groups, std::vector<int>(max_lag + 1, 0));
}

/*
 * Accumulate a squared displacement (dr2) for a given lag index.
 *
 * Adds dr2 to msd_sum[lag] and increments counts[lag].
 * Throws an exception if lag exceeds allocated size.
 *
 * Parameters:
 *   - lag: index corresponding to the time lag
 *   - dr2: squared displacement for this lag
 */
void MSDAccumulator::accumulate(int group, int lag, double dr2) {
  if (group >= msd_sum.size() || lag >= msd_sum[group].size()) {
    throw std::runtime_error(
        "MSDAccumulator::accumulate - Invalid group or lag index");
  }
  msd_sum[group][lag] += dr2;
  counts[group][lag] += 1;
}

/*
 * Finalize the MSD calculation.
 *
 * After all accumulations, this averages the summed squared displacements
 * by dividing each msd_sum[lag] by counts[lag], resulting in the mean
 * squared displacement for each lag.
 *
 * This is called once after accumulation is complete.
 */
void MSDAccumulator::finalize() {
  if (finalized)
    return;
  for (size_t g{0}; g < msd_sum.size(); ++g) {
    for (size_t i{0}; i < msd_sum[g].size(); ++i) {
      if (counts[g][i] > 0) {
        msd_sum[g][i] /= counts[g][i];
      }
    }
  }
  finalized = true;
}

/*
 * Write MSD results to a text file.
 *
 * The output file will contain three columns per line:
 *   - Lag index
 *   - Mean squared displacement at that lag
 *   - Number of samples (counts) contributing to that lag
 *
 * Parameters:
 *   - output_file: path to the output file (e.g., "output/msd.dat")
 */
void MSDAccumulator::write(const std::string &base_filename,
                           const std::vector<std::string> &group_labels,
                           double dt) const {

    for (size_t g = 0; g < msd_sum.size(); ++g) {
        // Generate file name: e.g., base + "_1_2_3.dat"
        std::string output_file = base_filename + "_" + group_labels[g] + ".dat";

        std::ofstream out(output_file);
        if (!out) {
            throw std::runtime_error("Failed to open file for writing MSD");
        }

        out << "# Fp-Tools MSD data file\n";
        out << "# Group: " << group_labels[g] << "\n";
        out << "# Time\t|\tMSD\n";
        out << std::scientific << std::setprecision(15);

        for (size_t lag = 0; lag < msd_sum[g].size(); ++lag) {
            double time = lag * dt;
            double avg_msd = msd_sum[g][lag];
            out << time << "\t" << avg_msd << "\n";
        }
        out.close();
    }
}

/*
 * Plot the MSD curve by calling an external Python script.
 *
 * Invokes 'plot_msd.py' with two arguments:
 *   - input_file: the MSD data file (e.g., "msd.dat")
 *   - format: desired output image format (e.g., "png", "pdf")
 *
 * Emits a warning if the Python script fails to execute.
 *
 * Parameters:
 *   - input_file: path to the MSD data file
 *   - format:     image format to generate
 */

void MSDAccumulator::plot(
    const std::vector<std::vector<std::string>> &input_file_groups,
    const std::vector<std::vector<std::string>> &label_groups,
    const std::string &output_dir, const std::string &format) const {

  const std::string script_path = std::string(SOURCE_DIR) + "/scripts/fp_plot.py";

  // If no plotting groups are provided, fall back to plotting "__all__"
  if (input_file_groups.empty() && label_groups.empty()) {
    std::string command =
        "python3 " + script_path + " msd 0 " + output_dir + "/msd_all.dat \"All\" " +
        output_dir + " " + format;

    int result = std::system(command.c_str());
    if (result != 0) {
      std::cerr << "Warning: Failed to execute MSD plotting script for __all__.\n";
    }
    return;
  }

  for (size_t i = 0; i < input_file_groups.size(); ++i) {
    std::string command =
        "python3 " + script_path + " msd " + std::to_string(i);

    const auto &files = input_file_groups[i];
    const auto &labels = label_groups[i];

    for (size_t j = 0; j < files.size(); ++j) {
      command += " " + files[j] + " \"" + labels[j] + "\"";
    }

    command += " " + output_dir + "/fig" + " " + format;

    int result = std::system(command.c_str());
    if (result != 0) {
      std::cerr << "Warning: Failed to execute MSD plotting script for group "
                << i << ".\n";
    }
  }
}


void MSDAccumulator::clear() {
  for (auto &group_vec : msd_sum) {
    std::fill(group_vec.begin(), group_vec.end(), 0.0);
  }
  for (auto &group_vec : counts) {
    std::fill(group_vec.begin(), group_vec.end(), 0);
  }
}

void MSDAccumulator::merge_group(const MSDAccumulator &other, int group_index) {
  int num_lags = static_cast<int>(other.msd_sum[0].size());
  for (int lag = 0; lag < num_lags; ++lag) {
    this->msd_sum[group_index][lag] += other.msd_sum[0][lag];
    this->counts[group_index][lag] += other.counts[0][lag];
  }
}

MSDAccumulator compute_msd(const std::vector<Frame> &frames, int max_lag,
                           std::vector<std::vector<std::string>> groups,
                           Logger *logger) {
  if (frames.empty() || groups.empty())
    return MSDAccumulator(0);

  int max_allowed_lag = static_cast<int>(frames.size()) - 1;
  if (max_lag <= 0 || max_lag > max_allowed_lag) {
    max_lag = max_allowed_lag;
  }

  MSDAccumulator accumulator(max_lag, static_cast<int>(groups.size()));

  // Build all atom trajectories once
  std::unordered_map<int, std::vector<Atom>> atom_trajectories;
  for (const Frame &frame : frames) {
    for (const Atom &atom : frame.atoms) {
      atom_trajectories[atom.id].push_back(atom);
    }
  }

  std::size_t group_index = 0;
  for (const auto &group : groups) {
    // Select matched trajectories for this group
    std::vector<const std::vector<Atom> *> matched_trajectories;
    matched_trajectories.reserve(atom_trajectories.size());

    if (group.size() == 1 && group[0] == "all") {
      for (const auto &[atom_id, traj] : atom_trajectories) {
        if (!traj.empty()) {
          matched_trajectories.push_back(&traj);
        }
      }
    } else if (group.size() == 1) {
      const std::string &target_type = group[0];
      for (const auto &[atom_id, traj] : atom_trajectories) {
        if (!traj.empty() && traj[0].type == target_type) {
          matched_trajectories.push_back(&traj);
        }
      }
    } else {
      std::unordered_set<std::string> type_set(group.begin(), group.end());
      for (const auto &[atom_id, traj] : atom_trajectories) {
        if (!traj.empty() && type_set.count(traj[0].type)) {
          matched_trajectories.push_back(&traj);
        }
      }
    }

    std::size_t total = matched_trajectories.size();

#ifdef USE_OMP
    std::atomic<int> processed(0);

#pragma omp parallel
    {
      MSDAccumulator local_accumulator(max_lag, 1); // Per-thread accumulator

#pragma omp for schedule(dynamic)
      for (int idx = 0; idx < static_cast<int>(total); ++idx) {
        local_accumulator.clear();

        const auto &traj = *matched_trajectories[idx];
        int N = static_cast<int>(traj.size());

        for (int lag = 0; lag <= max_lag; ++lag) {
          int n_valid = N - lag;
          if (n_valid <= 0)
            continue;

          for (int i = 0; i < n_valid; ++i) {
            double dr2 = squared_displacement(traj[i + lag], traj[i]);
            local_accumulator.accumulate(0, lag, dr2);
          }
        }

        int current_done = ++processed;

        if (omp_get_thread_num() == 0) {
          print_progress_bar(current_done, total, 40,
                             "MSD (Group " + std::to_string(group_index) + ")",
                             "particles");
        }

#pragma omp critical
        {
          accumulator.merge_group(local_accumulator, group_index);
        }
      }

#pragma omp barrier
      if (omp_get_thread_num() == 0 && processed == total) {
        print_progress_bar(total, total, 40,
                           "MSD (Group " + std::to_string(group_index) + ")",
                           "particles");
      }
    }
#else
    // Serial fallback
    for (int idx = 0; idx < static_cast<int>(total); ++idx) {
      const auto &traj = *matched_trajectories[idx];
      int N = static_cast<int>(traj.size());

      // So instead of going all the way to the last frame as 
      for (int lag = 0; lag <= max_lag; ++lag) {
        int n_valid = N - lag;

        for (int i = 0; i < n_valid; ++i) {
          double dr2 = squared_displacement(traj[i + lag], traj[i]);
          accumulator.accumulate(group_index, lag, dr2);
        }
      }

      if (logger) {
        print_progress_bar(idx + 1, total, 40,
                           "MSD (Group " + std::to_string(group_index) + ")",
                           "particles");
      }
    }
#endif
    ++group_index;
  }

  accumulator.finalize();
  return accumulator;
}

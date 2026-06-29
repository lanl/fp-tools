/**
 * -----------------------------------------------------------------------------
 * @file msd.h
 * @brief Mean Squared Displacement (MSD) Accumulator and Analysis
 *
 * This header defines the MSDAccumulator class and the 'compute_msd' function
 * for calculating the Mean Squared Displacement (MSD) from molecular dynamics
 * trajectories. MSD tracks how far particles move, on average, as a function
 * of time.
 *
 * MSDAccumulator collects and averages squared displacements between particle 
 * positions at various time lags. The results can be written to disk or plotted 
 * via an external Python script.
 *
 * Provided functionality includes:
 *  - MSDAccumulator: store and average squared displacements by lag
 *  - compute_msd: iterate over a trajectory and accumulate MSDs
 *  - write: output MSD data to disk
 *  - plot: generate a plot of MSD vs time lag
 *
 * @author Leah Hartman (leahghartman) (LANL, University of Michigan)
 * -----------------------------------------------------------------------------
 */

#pragma once

#include "base.h"
#include "logger.h"

#include <map>
#include <string>
#include <vector>

/**
 * @class MSDAccumulator
 * @brief Collects and averages squared displacements to compute the MSD curve.
 *
 * The MSDAccumulator maintains running totals of squared displacements between 
 * frames separated by different time lags. After all values are accumulated, 
 * the data can be finalized to obtain average MSD values per lag. These values 
 * can be retrieved, written to disk, or plotted externally.
 */
class MSDAccumulator {
public:
  /**
   * @brief Construct a new MSDAccumulator.
   * @param max_lag Maximum time lag to store in the accumulator.
   */
  MSDAccumulator(int max_lag, int num_groups = 1);
  
  /**
   * @brief Accumulate a squared displacement value at a given lag.
   * @param lag Time lag (frame separation index).
   * @param dr2 Squared displacement between positions at the lag.
   */
  void accumulate(int group, int lag, double dr2);

  /**
   * @brief Finalize the MSD by averaging accumulated displacements per lag.
   */
  void finalize();
  void clear();

  /**
   * @brief Write the finalized MSD data to a file.
   * @param output_file Path to the output file.
   * @param dt Time step between trajectory frames.
   */
  void write(const std::string& output_file, const std::vector<std::string>& group_labels, double dt) const;

  /**
   * @brief Plot the finalized MSD curve using an external script.
   * @param output_file Path to the MSD data file (used by plot script).
   * @param output_dir Directory to save the resulting plot.
   * @param format Image format for the plot (e.g., "png", "svg").
   */
  void plot(const std::vector<std::vector<std::string>>& input_file_groups,
            const std::vector<std::vector<std::string>>& label_groups,
            const std::string& output_dir,
            const std::string& format) const;
  //void plot(const std::string& output_file, const std::string& output_dir, const std::string& format) const;
  void merge_group(const MSDAccumulator& other, int group_index);
  
  /**
   * @brief Get the finalized MSD values.
   * @return Const reference to the vector of average MSD values per lag.
   */
  const std::vector<double>& get_msd() const;

  /**
   * @brief Get the number of samples per lag.
   * @return Const reference to the vector of sample counts.
   */
  const std::vector<int>& get_counts() const;
  
private:
  std::vector<std::vector<double>> msd_sum;
  std::vector<std::vector<int>> counts;
  bool finalized{false};            ///< Indicates whether MSDs are averaged.
};

/**
 * @brief Compute the Mean Squared Displacement (MSD) from a trajectory.
 *
 * Iterates over a subset of trajectory frames and accumulates squared
 * displacements for all frame pairs separated by lag ≤ max_lag. The result is 
 * an MSDAccumulator containing average MSD values.
 *
 * @param frames Vector of Frame structs representing the trajectory.
 * @param max_lag Maximum time lag to consider (optional; if -1, use full range).
 * @return A finalized MSDAccumulator object containing MSD vs. lag.
 */
MSDAccumulator compute_msd(const std::vector<Frame>& frames, 
                           int max_lag = -1, 
                           std::vector<std::vector<std::string>> groups = {},
                           Logger* logger = nullptr);


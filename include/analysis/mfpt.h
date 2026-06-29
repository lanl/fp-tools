//  TODO: Add unit tests for FPTAccumulator logic
//  TODO: Add OpenMP support
//  TODO: Fix inconsistent MFPT behavior
//  TODO: Ensure this is actually comparable to Jerome's results
//  TODO: Refactor write_fpt() to support subdirectory creation
//  TODO: Finish documentation

/**
 * -----------------------------------------------------------------------------
 * @file mfpt.h
 * @brief First Passage Time (FPT) Accumulator and MFPT Analysis
 *
 * This header defines the FPTAccumulator class and related functions for
 * computing and analyzing First Passage times (FPTs) and Mean First Passage
 * Times (MFPTs) from molecular dynamics trajectories. The FPTAccumulator 
 * collects first-passage times -- the time it takes for particles to reach a 
 * specified displacement threshold (radius bin) and computes their averages
 * over multiple starting times.
 *
 * The resulting MFPT data can then be written to disk and plotted using an 
 * external Python script. This class does not include a plot() method because
 * raw FPT distributions are never plotted directly.
 *
 * Provided functions include:
 *  - compute_fpt: process a trajectory and calculate FPTs
 *  - write_fpt: write FPT/MFPT data to disk
 *  - plot_mfpt: plot the MFPT data file
 *  - compute_dofr: compute the logarithmic derivative D(r) of MFPT vs radius
 *  - plot_dofr: plot the computed D(r) curve
 *
 * @author Leah Hartman (leahghartman) (LANL, University of Michigan)
 * -----------------------------------------------------------------------------
 */

#pragma once

#include "base.h"

#include <map>
#include <string>
#include <vector>

/**
 * @class FPTAccumulator
 * @brief Collects and averages First Passage Times (FPTs) for different radii.
 *
 * FPTAccumulator tracks the times at which particles first reach a given 
 * displacement (radius bin) from their initial position. It accumulates
 * these times across multiple initial times and particles, computes average
 * first-passage times, and exposes them for writing or further analysis.
 *
 * Note: this class intentionally does not provide a plot() method because the 
 * raw individual FPT distributions are not plotted directly -- only the 
 * summarized MFPT data is plotted externally after writing to a file.
 */
class FPTAccumulator {
public:
  /**
   * @brief Construct a new FPTAccumulator.
   * @param num_bins Number of bins for the FPT histogram (corresponding to radius intervals).
   */
  FPTAccumulator(int num_bins);
  
  /**
   * @brief Accumulate a first-passage time into the specified bin.
   * @param bin_index Index of the radius bin.
   * @param fpt First-passage time to accumulate.
   */
  void accumulate(int bin_index, double fpt);
  
  /**
   * @brief Finalize the average first-passage time values per bin.
   * @param dt Time step between trajectory frames.
   */
  void finalize(double dt);

  /**
   * @brief Get the averaged first-passage time histogram.
   * @return A const reference to the vector of average FPTs per time bin.
   */
  const std::vector<double>& get_fpt_dist() const;

  /**
   * @brief Get the total number of first-passage events accumulated.
   * @return Total count of samples across all bins.
   */
  int get_total_counts() const;

  const std::vector<double> get_fpt_sum() const;
  const std::vector<int> get_counts() const;
  
private:
  std::vector<double> fpt_sum{};    ///< Sum of first-passage times per bin.
  std::vector<int> counts{};        ///< Number of samples accumulated per bin.
  std::vector<double> fpt_dist{};   ///< Averaged first-passage times per bin.
  int total_counts{0};              ///< Total number of first-passage samples.
  bool finalized{false};            ///< Indicates whether averages are finalized.
  double r_min;
  double r_max;
};

/**
 * @brief Compute First Passage Time (FPT) distributions over a trajectory.
 *
 * This function scans through the trajectory and builds FPT histograms for a set
 * of radius thresholds. It does not compute the mean itself — that is done in 
 * the FPTAccumulator’s finalize() method.
 *
 * @param frames Vector of trajectory frames.
 * @param dr Radius bin width.
 * @param r_max Maximum radius up to which FPTs will be calculated.
 * @param dt Time step between trajectory frames.
 * @return A vector of FPTAccumulator objects (one per radius bin).
 */
std::vector<FPTAccumulator> compute_fpt(const std::vector<Frame>& frames,
                             double dr,
                             double r_max,
                             double dt,
                             const std::string& species);

/**
 * @brief Write the FPT distributions and MFPT data to disk.
 *
 * This function writes:
 * - Individual FPT distributions per radius bin to files ftp_<rad>.out.
 * - A summary file mfpt.dat with average MFPT vs. radius.
 *
 * @param output_file Path to the mfpt.dat file.
 * @param fpt_accums Vector of FPTAccumulator objects.
 * @param dr Radius bin width.
 * @param dt Time step between trajectory frames.
 */
void write_fpt(const std::string& mfpt_file,
               const std::string& fpt_dir,
               const std::vector<FPTAccumulator>& fpt_accums,
               double dr,
               double dt,
               double box_volume,
               size_t total_atoms,
               const std::string& species);

/**
 * @brief Plot the mean first-passage time (MFPT) as a function of distance.
 *
 * Reads the summary file mfpt.dat and generates a plot of MFPT vs. radius.
 *
 * @param summary_file Path to mfpt.dat containing the average MFPT data.
 * @param output_dir Directory to save the plot image.
 * @param format Image file format (e.g., "png", "svg").
 */
void plot_mfpt(const std::string& summary_file,
                       const std::string& output_dir,
                       const std::string& format);

/**
 * @brief Compute the logarithmic derivative D(r) of the MFPT curve.
 *
 * This estimates the scaling exponent D(r) from the relation T(r) ~ r^D(r).
 *
 * @param mfpt_file Path to mfpt.dat containing the average MFPT data.
 * @param output_file Path where the computed D(r) data will be written.
 */
void compute_dofr(const std::string& mfpt_file, const std::string& output_file);

/**
 * @brief Plot the scaling exponent D(r) vs. distance.
 *
 * Reads the computed D(r) file and generates a plot of D(r) vs. radius.
 *
 * @param log_derivative_file Path to the D(r) data file.
 * @param output_dir Directory to save the plot image.
 * @param format Image file format (e.g., "png", "svg").
 */
void plot_dofr(const std::string& log_derivative_file,
                         const std::string& output_dir,
                         const std::string& format);


#pragma once

#include "base.h"
#include "logger.h"

#include <string>
#include <unordered_map>
#include <vector>

class RDF {
public:
  // Constructors
  RDF(double dr, double r_max);

  // Main accumulation function; user provides pairs specifying type combinations
  void calculate(const std::vector<Frame>& frames,
                 const std::vector<std::vector<std::string>>& pairs, 
                 Logger* logger);

  // Normalize accumulated histograms to produce g(r)
  void normalize(double box_volume, size_t num_frames,
                 const std::unordered_map<std::string, size_t>& atoms_per_type);

  // Write RDF data files for each pair
  void write(const std::string& output_dir, double& min_box_length) const;

  // Optional: plot RDF curves (calls external python script)
  void plot(const std::vector<std::vector<std::string>>& input_file_groups,
                          const std::vector<std::vector<std::string>>& label_groups,
                          const std::string& output_dir,
                          const std::string& format) const;

private:
  double dr_{};
  double r_max_{};
  size_t num_bins_{};
  double ws_radius_{};

  std::vector<double> r_vals_;
  std::vector<std::vector<double>> histogram_;
  std::vector<std::vector<double>> rdf_;
  std::vector<std::vector<std::string>> pairs_;
};

// Factory functions to compute RDF by specifying dr or num_bins
RDF build_rdf(const std::vector<Frame>& frames,
                           double dr,
                           double r_max,
                           double box_volume, 
                           std::vector<std::vector<std::string>> pairs,
                           Logger* logger);

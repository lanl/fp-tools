/*
* -----------------------------------------------------------------------------
* Configuration Structures
* -------------------------
* Defines the data structures used to represent the user input configuration
* for trajectory analysis. Supports input, output, and various analysis modules
* including MSD, RDF, FKT, FSKT, MFPT, and D(r).
*
* These structs are designed to be populated by parsing a YAML configuration file.
*
* Authors: leahghartman
* -----------------------------------------------------------------------------
*/

#pragma once

#include <string>
#include <vector>
#include <array>

// --- Input Configuration ---
struct InputConfig {
  std::string file{};           // Path to the input trajectory file; required
  std::string vel_file{};
  std::string coord_type{};     // Coordinate type, e.g. "wrapped" or "unwrapped"; optional
  int start_frame{0};           // Starting frame index (inclusive); default 0
  int end_frame{-1};            // Ending frame index (inclusive); default -1 (last frame) 
  int frame_interval{1};        // Interval to sample frames, e.g. 1 means every frame
};

// --- Output Configuration ---
struct OutputConfig {
  std::string path{"output"};  //
};

struct PropertiesConfig {
  std::string type{""};
  std::array<double, 3> box_lengths{0.0, 0.0, 0.0};
  double temp{};
  double density{};
  double dt{};
};

struct Curve {
  std::vector<std::string> types;  // Atom types in this curve
  std::string label;               // Label for plotting
};

// --- MSD Configuration ---
struct MSDConfig {
  bool enabled{false};
  double dt{};
  int max_lag{-1};            // Maximum lag time (frames) for MSD calculation; -1 means full length
  std::vector<std::vector<std::string>> groups {};
  std::string plot_format{};  // Format for output plot files, e.g. "png"; empty disables plotting
  std::vector<std::vector<Curve>> curves{};
};

// --- RDF Configuration ---
struct RDFConfig {
  bool enabled{false};
  double dr{};
  double r_max{};
  int num_bins{};
  std::vector<std::vector<std::string>> pairs{};
  std::vector<std::vector<Curve>> curves{};
  std::string plot_format{};
};

// --- MFPT Configuration ---
struct MFPTConfig {
    bool enabled{false};
    double dt{};
    double dr{};
    double r_max{};
    std::string plot_format{};
    std::vector<std::string> species{};
};

struct CustomCorr {
  std::string name;
  std::vector<std::string> files;
  std::vector<int> columns;
};

struct CorrConfig {
    bool enabled {};
    double dt{};
    std::string n_list{};

    bool normalize {};
    bool subtract_mean {};
    std::string zero_pad {};
    std::string plot_format {};

    std::vector<std::string> builtins {};
    std::array<double, 3> n{0.0, 0.0, 0.0};
    std::vector<CustomCorr> custom {};
};

// --- Container for ALL analysis configurations ---
struct AnalysisConfig {
  MSDConfig msd{};    // MSD parameters
  RDFConfig rdf{};    // RDF parameters
  MFPTConfig mfpt{};  // MFPT parameters
  CorrConfig corr{};
};

// --- Top-level configurationi struct aggregating all sections ---
struct Config {
  InputConfig input{};        // Input file and frame selection options
  OutputConfig output{};      // Output verbosity and status options
  PropertiesConfig properties{};
  AnalysisConfig analysis{};  // Analysis-specific parameters
};

// Parser function that reads from a TOML file and returns a populated Config struct
Config parse_config(const std::string& filename);

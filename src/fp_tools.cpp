// TODO: add a warning in case the user tries to specify more frames than are

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <string_view>
#include <unistd.h>
#include <vector>

// Core data structures and configuration parsing
#include "traj_prov.h"
#include "config.h"

// Input trajectory parsers
#include "binary.h"
#include "lammps.h"
#include "vasp.h"
#include "xyz.h"

// Analysis modules
#include "corr.h"
#include "mfpt.h"
#include "msd.h"
#include "rdf.h"

// Configuration parser and serializer for C++
#include "toml.hpp"

// Log files
#include "logger.h"

namespace fs = std::filesystem;
char hostname[1024];

std::string get_timestamp_string() {
  std::time_t t = std::time(nullptr);
  char buf[100];
  std::strftime(buf, sizeof(buf), "%F %T", std::localtime(&t));
  return std::string(buf);
}

// std::string infer_format(const std::string& filename) {
//     std::ifstream f(filename);
//     if (!f.is_open()) {
//         throw std::runtime_error("Cannot open trajectory file: " + filename);
//     }
//
//     std::string first_line;
//     std::getline(f, first_line);
//     f.close();
//
//     char header[16] = {0};
//     f.read(header, 16);
//     std::streamsize bytes_read = file.gcount();
//     file.close(); // Cleanly close immediately after inspection
//
//     if (bytes_read < 4) return TrajectoryType::UNKNOWN;
//
//     // --- CHECK FOR BINARY DCD (FORTRAN UNFORMATTED) ---
//     // Your explorer script proved the first 4 bytes are exactly: 04 00 00 00
//     if (header[0] == 0x04 && header[1] == 0x00 && header[2] == 0x00 && header[3] == 0x00) {
//         return "binary";
//     }
//     
//     if (filename.find("XDATCAR") != std::string::npos) return "vasp";
//     if (first_line.find("ITEM: TIMESTEP") != std::string::npos) return "lammps";
//
//     char* end;
//     const long n = std::strtol(first_line.c_str(), &end, 10);
//     if (end != first_line.c_str() && n > 0) return "xyz";
// }


std::string infer_format(const std::string& filename) {
    std::ifstream f(filename, std::ios::binary);
    if (!f.is_open()) {
        throw std::runtime_error("Cannot open trajectory file: " + filename);
    }

    int first_byte = f.peek();
    
    if (first_byte == 0x04) {
        f.close();
        return "binary";
    }

    std::string first_line;
    std::getline(f, first_line);
    f.close();

    if (filename.find("XDATCAR") != std::string::npos) return "vasp";
    if (first_line.find("ITEM: TIMESTEP") != std::string::npos) return "lammps";

    char* end;
    const long n = std::strtol(first_line.c_str(), &end, 10);
    if (end != first_line.c_str() && n > 0) return "xyz";

    return "unknown"; // Safe fallback if no match is found
}


int main(int argc, char **argv) {
    // Check for required configuration call format
    if (argc < 2) {
        std::cerr << "Usage: fptools <config.toml>\n";
        return 1;
    }

    double full_elapsed = 0.0;
    auto full_start = std::chrono::steady_clock::now();

    std::string config_file{argv[1]};
    Config config = parse_config(config_file);

    // Name all of the possible analysis directories, but don't create them unless
    // asked
    std::string output_dir = config.output.path;
    std::string msd_dir = output_dir + "/msd";
    std::string rdf_dir = output_dir + "/rdf";
    std::string mfpt_dir = output_dir + "/mfpt";
    std::string corr_dir = output_dir + "/corr";

    // Create main output directory
    if (!fs::exists(output_dir) && !fs::create_directories(output_dir)) {
        std::cerr << "Error: Could not create output directory: " << output_dir
                  << "\n";
        return 1;
    }

    // Create a log file that will contain everything the terminal does
    Logger logger(output_dir + "/run.log");
    gethostname(hostname, 1024);

    logger.banner("0.1.0");

    logger.sys_info(hostname, std::filesystem::current_path(),
                    get_timestamp_string());
    logger.newline();

    // Parsing config
    logger.arrow("Parsing configuration file: " + config_file);

    // Infer trajectory format
    std::string format = infer_format(config.input.file);
    logger.arrow("Detected input trajectory format: " + format);
    logger.newline();

    // Instantiate trajectory parser
    logger.section("Trajectory Information");
    logger.newline();
    logger.arrow("Opening and parsing trajectory file: " + config.input.file);

    // ---------------------------------------------------
    // ---  Instantiate appropriate trajectory parser  ---
    // ---------------------------------------------------
    std::unique_ptr<TrajectoryProvider> traj{};

    if (format == "lammps") {
        traj = std::make_unique<LAMMPSTrajectory>(config.input.file);
    } else if (format == "xyz") {
        traj = std::make_unique<XYZTrajectory>(config.input.file);
    } else if (format == "binary") {
        traj = std::make_unique<BinaryTrajectory>(config.input.file);
    } else if (format == "vasp") {
        traj = std::make_unique<VASPTrajectory>(config.input.file);
    }

    // -----------------------------------------------------
    // ---  Load frames according to config information  ---
    // -----------------------------------------------------
    std::vector<Frame> frames{};
    int frame_count{0};
    Frame frame{};
    std::size_t loaded_frames{0};

    int frame_span =
        (config.input.end_frame >= 0)
            ? config.input.end_frame - config.input.start_frame + 1
            : std::numeric_limits<int>::max(); // fallback if no end frame

    int expected_frames = (frame_span == std::numeric_limits<int>::max())
                            ? -1
                            : (frame_span + config.input.frame_interval - 1) /
                                  config.input.frame_interval;

    bool use_config_box = false;
    double box_volume = 1.0;
    while (traj->parse_snapshot(frame)) {
        if (config.input.end_frame >= 0 && frame_count > config.input.end_frame)
            break;

        if (frame_count >= config.input.start_frame &&
            ((frame_count - config.input.start_frame) %
                config.input.frame_interval ==
            0)) {
            
            if (loaded_frames == 0) {
                if (frame.is_triclinic) {
                    throw std::runtime_error(
                        "Triclinic simulation boxes are not yet supported.\n"
                        "Please use an orthogonal box or contact the developers.");
                }

                // Check if the file provided box bounds
                const bool file_has_bounds = 
                    (frame.box_bounds[0][1] - frame.box_bounds[0][0]) > 0.0;
                
                if (!file_has_bounds) {
                    const auto& bl = config.properties.box_lengths;
                    if (bl[0] <= 0.0 || bl[1] <= 0.0 || bl[2] <= 0.0) {
                        throw std::runtime_error(
                            "Box bounds are not present in the trajectory file "
                            "and no box_lengths were specified in [system].\n"
                            "Add box_lengths = [lx, ly, lz] to your input deck.");
                    }
                    use_config_box = true;
                    logger.arrow("Box bounds set from [system].box_lengths.");
                }
            }

            // Override box bounds from config if needed (apply to every frame)
            if (use_config_box) {
                const auto& bl = config.properties.box_lengths;
                frame.box_bounds[0][0] = 0.0; frame.box_bounds[0][1] = bl[0];
                frame.box_bounds[1][0] = 0.0; frame.box_bounds[1][1] = bl[1];
                frame.box_bounds[2][0] = 0.0; frame.box_bounds[2][1] = bl[2];
            }

            // Compute box volume from first frame
            if (loaded_frames == 0) {
                box_volume = (frame.box_bounds[0][1] - frame.box_bounds[0][0])
                           * (frame.box_bounds[1][1] - frame.box_bounds[1][0])
                           * (frame.box_bounds[2][1] - frame.box_bounds[2][0]);
            }

            frames.push_back(frame);
            ++loaded_frames;

            if (expected_frames > 0 &&
                (loaded_frames % 10 == 0 || loaded_frames == expected_frames)) {
                print_progress_bar(loaded_frames, expected_frames, 40, "Loading frames",
                                    "frames");
            }
        }
        ++frame_count;
    }

    // if (!frames.empty()) {
    //   const auto &bounds = frames[0].box_bounds;
    //   double lx = bounds[0][1] - bounds[0][0];
    //   double ly = bounds[1][1] - bounds[1][0];
    //   double lz = bounds[2][1] - bounds[2][0];
    //   box_volume = lx * ly * lz;
    // }

    // logger.newline();
    logger.arrow("Trajectory parsing completed successfully!");
    logger.newline();
    logger.traj_info(config.input, frames);

    // ------------------------------------
    // ---  Perform requested analyses  ---
    // ------------------------------------

    bool ran_msd = false;
    double msd_dt = 0.0;
    int msd_max_lag = 0;
    int msd_num_groups = 0;
    int msd_num_frames = static_cast<int>(frames.size());
    double msd_elapsed = 0.0;
    std::vector<std::vector<std::string>> msd_groups;

    if (config.analysis.msd.enabled) {
        ran_msd = true;
        msd_dt = config.properties.dt;
        msd_max_lag = config.analysis.msd.max_lag;
        msd_groups = config.analysis.msd.groups;
        msd_num_groups = static_cast<int>(msd_groups.size());
        double dt = config.properties.dt;
        int max_lag = config.analysis.msd.max_lag;
        std::string plot_format = config.analysis.msd.plot_format;
        const auto &groups = config.analysis.msd.groups;
        const auto &curve_sets = config.analysis.msd.curves;

        std::string msd_fig_dir = msd_dir + "/fig";

        fs::create_directories(msd_dir);
        fs::create_directories(msd_fig_dir);

        logger.msd_info(dt, max_lag, groups);
        auto msd_start = std::chrono::steady_clock::now();

        MSDAccumulator msd_accum(max_lag);
        msd_accum = compute_msd(frames, max_lag, groups, &logger);

        // Build group labels like "1_2_3"
        std::vector<std::string> group_labels;
        for (const auto &group : groups) {
            std::ostringstream oss;
            for (size_t i = 0; i < group.size(); ++i) {
                if (i > 0)
                    oss << "-";
                oss << group[i];
            }
            group_labels.push_back(oss.str());
        }

        // Write .dat files with label-based names
        logger.arrow("Writing MSD data to " + msd_dir + "...");
        msd_accum.write(msd_dir + "/msd", group_labels, config.properties.dt);

        std::vector<std::vector<std::string>> input_file_groups;
        std::vector<std::vector<std::string>> label_groups;

        for (const auto &curve_group : curve_sets) {
            std::vector<std::string> file_group;
            std::vector<std::string> label_group;

            for (const auto &curve : curve_group) {
                std::string types_label;

                for (size_t i = 0; i < curve.types.size(); ++i) {
                    if (i > 0)
                        types_label += "-";
                    types_label += curve.types[i];
                }
                std::string filename = msd_dir + "/msd_" + types_label + ".dat";
                file_group.push_back(filename);
                label_group.push_back(curve.label);
            }
            input_file_groups.push_back(file_group);
            label_groups.push_back(label_group);
        }

        // Plot results
        if (!plot_format.empty()) {
            logger.arrow("Generating MSD plot(s) (" + plot_format + ")...");
            msd_accum.plot(input_file_groups, label_groups, msd_dir, plot_format);
            logger.arrow("Plot(s) successfully saved to " + msd_dir + "/!");
        }
        auto msd_end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = msd_end - msd_start;
        msd_elapsed = elapsed_seconds.count();
        logger.newline();
    }

    // Radial Density Function (RDF) calculation
    std::vector<std::vector<std::string>> rdf_pairs = config.analysis.rdf.pairs;
    int rdf_num_frames = static_cast<int>(frames.size());
    double rdf_dt = config.properties.dt;
    double rdf_elapsed = 0.0;
    double rdf_dr = config.analysis.rdf.dr;
    int rdf_num_bins = config.analysis.rdf.num_bins;
    double rdf_r_max = config.analysis.rdf.r_max;
    bool ran_rdf = false;

    if (config.analysis.rdf.enabled) {
        ran_rdf = true;
        double dr = config.analysis.rdf.dr;
        int num_bins = config.analysis.rdf.num_bins;
        double r_max = config.analysis.rdf.r_max;
        const auto &pairs = config.analysis.rdf.pairs;
        const auto &curve_sets = config.analysis.rdf.curves;
        std::string plot_format = config.analysis.rdf.plot_format;

        double xlo = frames[0].box_bounds[0][0], xhi = frames[0].box_bounds[0][1];
        double ylo = frames[0].box_bounds[1][0], yhi = frames[0].box_bounds[1][1];
        double zlo = frames[0].box_bounds[2][0], zhi = frames[0].box_bounds[2][1];
        double lx = xhi - xlo, ly = yhi - ylo, lz = zhi - zlo;
        double min_box_length = std::min<double>({lx, ly, lx});

        std::string rdf_fig_dir = rdf_dir + "/fig";
        fs::create_directories(rdf_dir);
        fs::create_directories(rdf_fig_dir);

        logger.rdf_info(dr, num_bins, r_max, pairs);
        auto rdf_start = std::chrono::steady_clock::now();

        // Ensure atom types exist or default to "1"
        bool has_missing_types =
            frames[0].atoms.empty() || frames[0].atoms[0].type.empty();
        if (has_missing_types) {
            for (auto &frame : const_cast<std::vector<Frame> &>(frames)) {
                for (auto &atom : frame.atoms) {
                    atom.type = "1";
                }
            }
        }

        // Compute RDF
        RDF rdf_accum = RDF(dr, r_max);

        rdf_accum = build_rdf(frames, dr, r_max, box_volume, pairs, &logger);

        logger.arrow("Writing RDF data to " + rdf_dir + "...");
        rdf_accum.write(rdf_dir, min_box_length);

        // Plot RDF results if requested
        if (!plot_format.empty() && !curve_sets.empty()) {
            logger.arrow("Generating RDF plot(s) (" + plot_format + ")...");

            std::vector<std::vector<std::string>> input_file_groups;
            std::vector<std::vector<std::string>> label_groups;

            for (const auto &curve_group : curve_sets) {
                std::vector<std::string> file_group;
                std::vector<std::string> label_group;

                for (const auto &curve : curve_group) {
                    std::ostringstream oss_label;
                    for (size_t i = 0; i < curve.types.size(); ++i) {
                        if (i > 0)
                            oss_label << "_";
                        oss_label << curve.types[i];
                    }
                    std::string filename = rdf_dir + "/rdf-" + oss_label.str() + ".dat";
                    file_group.push_back(filename);
                    label_group.push_back(curve.label);
                }

                input_file_groups.push_back(file_group);
                label_groups.push_back(label_group);
            }
            rdf_accum.plot(input_file_groups, label_groups, rdf_dir, plot_format);
            logger.arrow("Plot(s) successfully saved to " + rdf_dir + "/!");
        }

        auto rdf_end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = rdf_end - rdf_start;
        rdf_elapsed = elapsed_seconds.count();
        logger.newline();
    }

    bool ran_mfpt = false;
    double mfpt_elapsed = 0.0;
    int mfpt_num_radii = 0;
    int mfpt_num_frames = static_cast<int>(frames.size());
    double mfpt_dt = config.analysis.mfpt.dt;
    double mfpt_dr = config.analysis.mfpt.dr;
    double mfpt_r_max = config.analysis.mfpt.r_max;

    if (config.analysis.mfpt.enabled) {
        ran_mfpt = true;

        // Extract parameters
        double dr = config.analysis.mfpt.dr;
        double r_max = config.analysis.mfpt.r_max;
        double dt = config.analysis.mfpt.dt;
        std::string plot_format = config.analysis.mfpt.plot_format;

        std::vector<std::string> multi_species_targets = config.analysis.mfpt.species;

        // Helper to generate unique filenames based on species
        auto get_species_suffix = [](const std::string& sp) {
            if (sp.empty()) return std::string("_all");
            std::string suffix = "_" + sp;
            return suffix;
        };

        // Info + timer
        logger.mfpt_info(dr, r_max, dt);
        auto mfpt_start = std::chrono::steady_clock::now();

        for (const std::string& target_species : multi_species_targets) {
       
            std::string suffix = get_species_suffix(target_species);
            logger.arrow("--- Starting MFPT Analysis for species: " + target_species + " ---");

            // Filter the atom count
            size_t species_atom_count = 0;
            if (target_species.empty() || target_species == "all") {
                species_atom_count = frames[0].atoms.size();
            } else {
                for (const auto& atom : frames[0].atoms) {
                    if (atom.type == target_species) {
                        species_atom_count++;
                    }
                }
            }

            if (species_atom_count == 0) {
                logger.arrow("Skipping species " + target_species + " (0 atoms found)");
                continue;
            }

            // Dynamically name outputs so nothing gets overwritten
            std::string fpt_dir = mfpt_dir + "/fpt";
            std::string mfpt_file = mfpt_dir + "/mfpt" + suffix + ".dat";
            std::string dr_file = output_dir + "/D_r" + suffix + ".dat";

            fs::create_directories(fpt_dir);

            // Print details for this specific run
            logger.mfpt_info(dr, r_max, dt);

            // 1. Compute MFPT for the current species subset
            std::vector<FPTAccumulator> fpt_accums = compute_fpt(frames, dr, r_max, dt, target_species);
        
            // 2. Write the distribution files (fpt_X_species.out) and summary map (mfpt_species.dat)
            logger.arrow("Writing MFPT data to " + mfpt_file + "...");
            write_fpt(mfpt_file, fpt_dir, fpt_accums, dr, dt, box_volume, species_atom_count, target_species);
    
            // 3. Extract the Log-Derivative D(r) from the file we just wrote
            logger.arrow("Computing D(r) from MFPT...");
            compute_dofr(mfpt_file, dr_file);

            // 4. Generate Python plots specifically tailored to this species run
            // logger.arrow("Plotting D(r)...");
            // plot_dofr(dr_file, mfpt_dir, "svg");
            //
            // if (!plot_format.empty()) {
            //     logger.arrow("Generating MFPT plot(s) (" + plot_format + ")...");
            //     plot_mfpt(mfpt_file, mfpt_dir, plot_format);
            //     logger.arrow("Plot(s) successfully saved to " + mfpt_dir + "/!");
            // }
        
            logger.newline();
        }

        // // Compute and write MFPT data
        // std::vector<FPTAccumulator> fpt_accums = compute_fpt(frames, dr, r_max, dt, target_species);
        // logger.arrow("Writing MFPT data to " + mfpt_dir + "...");
        // write_fpt(mfpt_file, fpt_dir, fpt_accums, dr, dt, box_volume,
        //           frames[0].atomnum, target_species);
        //
        // // Compute D(r)
        // logger.arrow("Computing D(r) from MFPT...");
        // compute_dofr(mfpt_file, dr_file);
        //
        // // logger.arrow("Plotting D(r)...");
        // // plot_dofr(dr_file, mfpt_dir, "svg");
        //
        // // if (!plot_format.empty()) {
        // //   logger.arrow("Generating MFPT plot(s) (" + plot_format + ")...");
        // //   plot_mfpt(mfpt_file, mfpt_dir, plot_format);
        // //   logger.arrow("Plot(s) successfully saved to " + mfpt_dir + "/!");
        // // }

        auto mfpt_end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = mfpt_end - mfpt_start;
        double mfpt_elapsed = elapsed_seconds.count();
        mfpt_num_radii = static_cast<int>(mfpt_r_max / mfpt_dr);

        logger.newline();
    }

    if (config.analysis.corr.enabled) {
        bool normalize{config.analysis.corr.normalize};
        bool subtract_mean{config.analysis.corr.subtract_mean};
        std::string zero_pad{config.analysis.corr.zero_pad};
        std::string plot_format{config.analysis.corr.plot_format};
        std::vector<std::string> builtins{config.analysis.corr.builtins};
        std::vector<CustomCorr> custom{config.analysis.corr.custom};

        double Lx = frames[0].box_bounds[0][1] - frames[0].box_bounds[0][0];
        double Ly = frames[0].box_bounds[1][1] - frames[0].box_bounds[1][0];
        double Lz = frames[0].box_bounds[2][1] - frames[0].box_bounds[2][0];
        constexpr double TWO_PI = 6.28318530717958647692;

        // --- RESOLVE Q-VECTORS (File List vs. Single Vector) ---
        std::vector<Vector3D> q_vectors;
        std::string n_list_path = config.analysis.corr.n_list;

        if (!n_list_path.empty()) {
            // Scenario A: User provided a file path containing multiple n-vectors
            std::ifstream n_file(n_list_path);
            if (!n_file.is_open()) {
                throw std::runtime_error("Could not open q-vector index file: " + n_list_path);
            }

            std::string line;
            while (std::getline(n_file, line)) {
                if (line.empty() || line[0] == '#') continue;

                std::istringstream iss(line);
                double nx, ny, nz;
                if (iss >> nx >> ny >> nz) {
                    Vector3D q {
                        (TWO_PI / Lx) * nx,
                        (TWO_PI / Ly) * ny,
                        (TWO_PI / Lz) * nz
                    };
                    q_vectors.push_back(q);
                }
            }
            n_file.close();

            if (q_vectors.empty()) {
                throw std::runtime_error("The q-vector index file was empty or incorrectly formatted: " + n_list_path);
            }
            logger.arrow("Loaded " + std::to_string(q_vectors.size()) + " q-vectors from file for shell-averaging.");
        } else {
            // Scenario B: Fallback to the single 'n' vector specified directly in the config
            double nx = config.analysis.corr.n[0];
            double ny = config.analysis.corr.n[1];
            double nz = config.analysis.corr.n[2];

            Vector3D q_vector{
                (TWO_PI / Lx) * nx,
                (TWO_PI / Ly) * ny,
                (TWO_PI / Lz) * nz
            };

            q_vectors.push_back(q_vector);
            logger.arrow("No q-vector list file provided. Using single vector n = [" 
                        + std::to_string(static_cast<int>(nx)) + ", " 
                        + std::to_string(static_cast<int>(ny)) + ", " 
                        + std::to_string(static_cast<int>(nz)) + "]");
        }

        std::vector<std::string> custom_names;
        for (const auto &c : custom) {
            custom_names.push_back(c.name);
        }

        fs::create_directories(corr_dir);
        std::string corr_file{corr_dir + "/corr.dat"};

        logger.corr_info(builtins, custom_names, zero_pad, normalize,
                         subtract_mean);

        CorrAccumulator corr_accum =
            compute_correlations(config.analysis.corr, frames, q_vectors);

        logger.arrow("Writing correlation data to " + corr_file);
        corr_accum.write(corr_dir, config.analysis.corr.dt);

        // Optionally plot correlation results using configured plot format
        if (!plot_format.empty()) {
            logger.arrow("Generating correlation plot (" + plot_format + ")...");
            corr_accum.plot(corr_dir, plot_format);
            logger.arrow("Plot successfully saved to " + corr_dir + "/!");
        }
        logger.newline();
    }

    auto full_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = full_end - full_start;

    logger.section("Summary");
    logger.newline();
    if (ran_msd) {
        logger.msd_summ(msd_groups, msd_num_frames, msd_dt, msd_elapsed);
    }
    if (ran_rdf) {
        logger.rdf_summ(rdf_pairs, rdf_num_frames, rdf_dt, rdf_elapsed, rdf_dr,
                        rdf_num_bins, rdf_r_max);
    }
    if (ran_mfpt) {
        logger.mfpt_summ(mfpt_num_radii, mfpt_num_frames, mfpt_dr, mfpt_r_max,
                        mfpt_dt, mfpt_elapsed);
    }

    std::ostringstream oss;
    oss << "[ TOTAL ] Total runtime: " << std::fixed << std::setprecision(2)
        << elapsed_seconds.count() << " seconds";
    logger.note(oss.str());
    logger.newline();
    return 0;
}

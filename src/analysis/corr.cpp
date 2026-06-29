
#include "corr.h"
#include "base.h"
#include "config.h"

#include <algorithm>
#include <cmath>
#include <fftw3.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <complex>

CorrAccumulator::CorrAccumulator() = default;
CorrAccumulator::~CorrAccumulator() = default;

// Helper: subtract mean in-place
void subtract_mean(std::vector<double> &data) {
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    for (auto &val : data)
        val -= mean;
}

// Helper: zero pad to nearest power of 2
size_t next_pow2(size_t n) {
    size_t pow2 = 1;
    while (pow2 < n)
        pow2 <<= 1;
    return pow2;
}

void CorrAccumulator::merge(const CorrAccumulator &other) {
    for (const auto &[key, val] : other.results) {
        results[key] = val;
    }
}

// Helper: read a specific column from file (1-based indexing)
std::vector<double> read_column(const std::string &path, size_t col_index) {
    std::ifstream infile(path);
    if (!infile.is_open()) {
        throw std::runtime_error("Could not open file: " + path);
    }

    std::vector<double> column_data;
    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#')
            continue;
        std::istringstream iss(line);
        double val;
        size_t current_col = 1;
        while (iss >> val) {
            if (current_col == col_index) {
                column_data.push_back(val);
                break;
            }
            ++current_col;
        }
    }
    return column_data;
}

// Core FFT-based correlation
std::vector<double> correlate(const std::vector<double> &a,
                              const std::vector<double> &b,
                              const CorrConfig &cfg) {
    size_t N = a.size();

    // Pad to avoid circular correlation stuff
    size_t pad_size = cfg.zero_pad == "next_pow2" ? next_pow2(2 * N) : 2 * N;

    // Zero-padded copies
    std::vector<double> A(pad_size, 0.0);
    std::vector<double> B(pad_size, 0.0);
    std::copy(a.begin(), a.end(), A.begin());
    std::copy(b.begin(), b.end(), B.begin());

    // Allocate FFTW buffers
    fftw_complex *fft_a =
        (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (pad_size / 2 + 1));
    fftw_complex *fft_b =
        (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (pad_size / 2 + 1));
    fftw_complex *fft_out =
        (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (pad_size / 2 + 1));

    // Create plans
    fftw_plan plan_a =
        fftw_plan_dft_r2c_1d(pad_size, A.data(), fft_a, FFTW_ESTIMATE);
    fftw_plan plan_b =
        fftw_plan_dft_r2c_1d(pad_size, B.data(), fft_b, FFTW_ESTIMATE);
    fftw_plan plan_back = 
        fftw_plan_dft_c2r_1d(pad_size, fft_out, A.data(), FFTW_ESTIMATE);

    // Forward transforms
    fftw_execute(plan_a);
    fftw_execute(plan_b);

    // For autocorrelations this is |FFT_a|^2
    for (size_t i = 0; i < pad_size / 2 + 1; ++i) {
        double re = fft_a[i][0] * fft_b[i][0] + fft_a[i][1] * fft_b[i][1];
        double im = fft_a[i][1] * fft_b[i][0] - fft_a[i][0] * fft_b[i][1];
        fft_out[i][0] = re;
        fft_out[i][1] = im;
    }

    // Inverse transform, the result goes back into A
    fftw_execute(plan_back);

    // Build the result vector and normalize by pad_size (FFTW doesn't normalize)
    std::vector<double> result(N);
    for (size_t i = 0; i < N; ++i) {
        result[i] = (A[i] / static_cast<double>(pad_size));
    }

    fftw_destroy_plan(plan_a);
    fftw_destroy_plan(plan_b);
    fftw_destroy_plan(plan_back);
    fftw_free(fft_a);
    fftw_free(fft_b);
    fftw_free(fft_out);

    return result;
}

CorrAccumulator compute_builtin(const CorrConfig &cfg,
                                const std::vector<Frame> &frames,
                                const std::vector<Vector3D>& q_vectors) {
    CorrAccumulator acc;
    size_t num_atoms = frames[0].atoms.size();
    size_t num_frames = frames.size();

    // Create a modified config that turns off normalization during component calculation
    CorrConfig internal_cfg = cfg;
    internal_cfg.normalize = false;
    internal_cfg.subtract_mean = false;


    // --- VACF CALCULATION (INDEP OF Q-STUFF) --------------------------------
    if (std::find(cfg.builtins.begin(), cfg.builtins.end(), "vacf") != cfg.builtins.end()) {
        std::vector<double> vacf(num_frames, 0.0);

        // Temporary vector to accumulate results from each atom
        std::vector<double> vacf_accum(num_frames, 0.0);

        for (size_t i = 0; i < num_atoms; ++i) {
            // Extract velocity components for atom i
            std::vector<double> vx(num_frames), vy(num_frames), vz(num_frames);
            for (size_t t = 0; t < num_frames; ++t) {
                vx[t] = frames[t].atoms[i].vx;
                vy[t] = frames[t].atoms[i].vy;
                vz[t] = frames[t].atoms[i].vz;
            }

            // Compute autocorrelations per component (raw, unnormalized)
            std::vector<double> c_x = correlate(vx, vx, internal_cfg);
            std::vector<double> c_y = correlate(vy, vy, internal_cfg);
            std::vector<double> c_z = correlate(vz, vz, internal_cfg);

            // Sum components to get total VACF for this atom
            for (size_t lag = 0; lag < num_frames; ++lag) {
                vacf_accum[lag] += c_x[lag] + c_y[lag] + c_z[lag];
            }
        }

        // Average over all atoms
        for (size_t lag = 0; lag < num_frames; ++lag) {
            vacf[lag] = vacf_accum[lag] / static_cast<double>(num_atoms);
        }

        // Actually normalize correctly at the end
        if (cfg.normalize && vacf[0] != 0.0) {
            const double norm = vacf[0];
            for (auto& v : vacf) v /= norm;
        }
        acc.store_result("vacf", vacf);
    }

    if (q_vectors.empty()) return acc;
    std::map<double, std::vector<Vector3D>> q_shells;
    for (const auto& q : q_vectors) {
        double q_mag = std::sqrt(q.x * q.x + q.y * q.y + q.z * q.z);
        
        bool found_shell = false;
        for (auto& [mag, shell] : q_shells) {
            if (std::abs(mag - q_mag) < 1e-3) { // Handle floating-point tolerance
                shell.push_back(q);
                found_shell = true;
                break;
            }
        }
        if (!found_shell) {
            q_shells[q_mag] = {q};
        }
    }

    // Vectors to hold the static structure factor curve points across shells
    std::vector<double> ssf_q_axis;
    std::vector<double> ssf_values;

    // LOOP THROUGH UNIQUE INDEPENDENT Q-SHELLS FOR COLLECTIVE PROPERTIES
    for (const auto& [q_mag, shell_vectors] : q_shells) {
        size_t num_q_in_shell = shell_vectors.size();
        
        // Formulate a clean string label for this specific shell coordinate
        std::string q_str = std::to_string(q_mag);
        // Truncate trailing string zeros to keep file names neat
        q_str.erase(q_str.find_last_not_of('0') + 1, std::string::npos);
        if (q_str.back() == '.') q_str.pop_back();

        // CALCULATOR 2: ISF / DSF CURVES FOR THIS SHELL
        if (std::find(cfg.builtins.begin(), cfg.builtins.end(), "dsf") != cfg.builtins.end()) {
            std::vector<double> F_qt_shell_avg(num_frames, 0.0);

            for (const auto& q_vec : shell_vectors) {
                std::vector<double> R(num_frames, 0.0);
                std::vector<double> I(num_frames, 0.0);

                for (size_t t = 0; t < num_frames; ++t) {
                    for (size_t j = 0; j < num_atoms; ++j) {
                        double phase = q_vec.x * frames[t].atoms[j].x +
                                       q_vec.y * frames[t].atoms[j].y +
                                       q_vec.z * frames[t].atoms[j].z;
                        R[t] += std::cos(phase);
                        I[t] += std::sin(phase);
                    }
                }

                std::vector<double> c_RR = correlate(R, R, internal_cfg);
                std::vector<double> c_II = correlate(I, I, internal_cfg);

                for (size_t lag = 0; lag < num_frames; ++lag) {
                    F_qt_shell_avg[lag] += ((c_RR[lag] + c_II[lag]) / static_cast<double>(num_atoms));
                }
            }

            // Average across vectors inside this particular magnitude group
            for (auto& val : F_qt_shell_avg) val /= static_cast<double>(shell_vectors.size());

            // Save static structure factor data points for our curve compilation
            ssf_q_axis.push_back(q_mag);
            ssf_values.push_back(F_qt_shell_avg[0]);

            // Dynamic Structure Factor S(q, w) calculation via real-to-complex FFT
            size_t num_frequencies = num_frames / 2 + 1;
            fftw_complex *dsf_complex = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * num_frequencies);
            std::vector<double> F_qt_copy = F_qt_shell_avg; 
            fftw_plan dsf_plan = fftw_plan_dft_r2c_1d(num_frames, F_qt_copy.data(), dsf_complex, FFTW_ESTIMATE);
            fftw_execute(dsf_plan);

            std::vector<double> S_qw(num_frequencies, 0.0);
            for (size_t w = 0; w < num_frequencies; ++w) S_qw[w] = dsf_complex[w][0];

            fftw_destroy_plan(dsf_plan);
            fftw_free(dsf_complex);
            acc.store_result("dsf_q_" + q_str, S_qw);

            if (cfg.normalize && F_qt_shell_avg[0] != 0.0) {
                const double norm = F_qt_shell_avg[0];
                for (auto &v : F_qt_shell_avg) v /= norm;
            }
            acc.store_result("isf_q_" + q_str, F_qt_shell_avg);

        }

        // CALCULATOR 3: LONGITUDINAL & TRANSVERSE CURRENT CORRELATIONS
        if (std::find(cfg.builtins.begin(), cfg.builtins.end(), "curr") != cfg.builtins.end()){
            std::vector<double> C_L_shell_avg(num_frames, 0.0);
            std::vector<double> C_T_shell_avg(num_frames, 0.0);

            for (const auto& q_vec : shell_vectors) {
                Vector3D q_hat = {0.0, 0.0, 0.0};
                if (q_mag > 0.0) {
                    q_hat = { q_vec.x / q_mag, q_vec.y / q_mag, q_vec.z / q_mag };
                }

                std::vector<double> R_x(num_frames, 0.0); std::vector<double> I_x(num_frames, 0.0);
                std::vector<double> R_y(num_frames, 0.0); std::vector<double> I_y(num_frames, 0.0);
                std::vector<double> R_z(num_frames, 0.0); std::vector<double> I_z(num_frames, 0.0);
                    
                for (size_t t = 0; t < num_frames; ++t) {
                    for (size_t j = 0; j < num_atoms; ++j) {
                        double phase = q_vec.x * frames[t].atoms[j].x +
                                       q_vec.y * frames[t].atoms[j].y +
                                       q_vec.z * frames[t].atoms[j].z;

                        double cos_p = std::cos(phase);
                        double sin_p = std::sin(phase);

                        R_x[t] += frames[t].atoms[j].vx * cos_p;
                        R_y[t] += frames[t].atoms[j].vy * cos_p;
                        R_z[t] += frames[t].atoms[j].vz * cos_p;

                        I_x[t] += frames[t].atoms[j].vx * sin_p;
                        I_y[t] += frames[t].atoms[j].vy * sin_p;
                        I_z[t] += frames[t].atoms[j].vz * sin_p;
                    }
                }

                std::vector<double> R_L(num_frames, 0.0);   std::vector<double> I_L(num_frames, 0.0);
                std::vector<double> R_T_x(num_frames, 0.0); std::vector<double> I_T_x(num_frames, 0.0);
                std::vector<double> R_T_y(num_frames, 0.0); std::vector<double> I_T_y(num_frames, 0.0);
                std::vector<double> R_T_z(num_frames, 0.0); std::vector<double> I_T_z(num_frames, 0.0);
                
                for (size_t t = 0; t < num_frames; ++t) {
                    R_L[t] = R_x[t] * q_hat.x + R_y[t] * q_hat.y + R_z[t] * q_hat.z;
                    I_L[t] = I_x[t] * q_hat.x + I_y[t] * q_hat.y + I_z[t] * q_hat.z;

                    R_T_x[t] = R_x[t] - R_L[t] * q_hat.x;
                    R_T_y[t] = R_y[t] - R_L[t] * q_hat.y;
                    R_T_z[t] = R_z[t] - R_L[t] * q_hat.z;

                    I_T_x[t] = I_x[t] - I_L[t] * q_hat.x;
                    I_T_y[t] = I_y[t] - I_L[t] * q_hat.y;
                    I_T_z[t] = I_z[t] - I_L[t] * q_hat.z;
                }

                std::vector<double> c_RL = correlate(R_L, R_L, internal_cfg);
                std::vector<double> c_IL = correlate(I_L, I_L, internal_cfg);
                std::vector<double> c_RTx = correlate(R_T_x, R_T_x, internal_cfg);
                std::vector<double> c_RTy = correlate(R_T_y, R_T_y, internal_cfg);
                std::vector<double> c_RTz = correlate(R_T_z, R_T_z, internal_cfg);
                std::vector<double> c_ITx = correlate(I_T_x, I_T_x, internal_cfg);
                std::vector<double> c_ITy = correlate(I_T_y, I_T_y, internal_cfg);
                std::vector<double> c_ITz = correlate(I_T_z, I_T_z, internal_cfg);

                for (size_t lag = 0; lag < num_frames; ++lag) {
                    C_L_shell_avg[lag] += (c_RL[lag] + c_IL[lag]) / static_cast<double>(num_atoms);
                    C_T_shell_avg[lag] += (c_RTx[lag] + c_ITx[lag] + 
                                       c_RTy[lag] + c_ITy[lag] + 
                                       c_RTz[lag] + c_ITz[lag]) / static_cast<double>(num_atoms);
                }
            }

            // Normalization shell-averaging step
            for (size_t lag = 0; lag < num_frames; ++lag) {
                C_L_shell_avg[lag] /= static_cast<double>(num_q_in_shell);
                C_T_shell_avg[lag] /= static_cast<double>(num_q_in_shell);
            }
               
            if (cfg.normalize) {
                if (C_L_shell_avg[0] != 0.0) {
                    double norm_L = C_L_shell_avg[0];
                    for (auto& v : C_L_shell_avg) v /= norm_L;
                }
                if (C_T_shell_avg[0] != 0.0) {
                    double norm_T = C_T_shell_avg[0];
                    for (auto& v : C_T_shell_avg) v /= norm_T;
                }
            }
            acc.store_result("cl_t_q_" + q_str, C_L_shell_avg);
            acc.store_result("ct_t_q_" + q_str, C_T_shell_avg);

            // Longitudinal & Transverse Frequency Spectra
            size_t num_frequencies = num_frames / 2 + 1;
            fftw_complex *curr_complex = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * num_frequencies);

            std::vector<double> C_L_t_copy = C_L_shell_avg; 
            fftw_plan plan_L = fftw_plan_dft_r2c_1d(num_frames, C_L_t_copy.data(), curr_complex, FFTW_ESTIMATE);
            fftw_execute(plan_L);
            std::vector<double> C_L_w(num_frequencies, 0.0);
            for (size_t w = 0; w < num_frequencies; ++w) C_L_w[w] = curr_complex[w][0];
            fftw_destroy_plan(plan_L);       
                
            std::vector<double> C_T_t_copy = C_T_shell_avg; 
            fftw_plan plan_T = fftw_plan_dft_r2c_1d(num_frames, C_T_t_copy.data(), curr_complex, FFTW_ESTIMATE);
            fftw_execute(plan_T);
            std::vector<double> C_T_w(num_frequencies, 0.0);
            for (size_t w = 0; w < num_frequencies; ++w) C_T_w[w] = curr_complex[w][0];
            fftw_destroy_plan(plan_T);   

            fftw_free(curr_complex);

            acc.store_result("cl_w_q_" + q_str, C_L_w);
            acc.store_result("ct_w_q_" + q_str, C_T_w);
        }
    }

    // PHASE 3: STORE ENTIRE STATIC STRUCTURE FACTOR DATA SET
    if (!ssf_q_axis.empty()) {
        acc.store_result("ssf_q_axis", ssf_q_axis);
        acc.store_result("ssf_values", ssf_values);
    }

    return acc;
}

CorrAccumulator compute_correlations(const CorrConfig &cfg,
                                     const std::vector<Frame> &frames,
                                     const std::vector<Vector3D>& q_vectors) {
    CorrAccumulator acc;

    // --- Built-in functions (vacf, vcf) ---
    if (!cfg.builtins.empty()) {
        CorrAccumulator builtin_acc = compute_builtin(cfg, frames, q_vectors);
        acc.merge(builtin_acc);
    }

    // --- Custom file/column pair correlations ---
    for (const auto &pair : cfg.custom) {
        std::vector<double> a, b;

        if (pair.files.size() == 1) {
            // Autocorrelation: use same file and column for both
            a = read_column(pair.files[0], pair.columns[0]);
            b = a;
        } else {
            a = read_column(pair.files[0], pair.columns[0]);
            b = read_column(pair.files[1], pair.columns[1]);
        }

        if (cfg.subtract_mean) {
            subtract_mean(a);
            // Only subtract mean from b if it's not a reference to a
            if (&a != &b)
                subtract_mean(b);
        }
        auto result = correlate(a, b, cfg);
        acc.store_result(pair.name, result);
    }
    return acc;
}

void CorrAccumulator::write(const std::string &output_dir, double dt) const {
    namespace fs = std::filesystem;
    fs::create_directories(output_dir);

    // 1. Process and isolate the global static structure factor curve first
    if (results.count("ssf_values") && results.count("ssf_q_axis")) {
        std::ofstream ssf_out(output_dir + "/ssf_curve.dat");
        ssf_out << "# q_magnitude S(q)\n";
        const auto& q_ax = results.at("ssf_q_axis");
        const auto& vals = results.at("ssf_values");
        for (size_t i = 0; i < vals.size(); ++i) {
            ssf_out << q_ax[i] << " " << vals[i] << "\n";
        }
    }

    // 2. Dynamically organize the rest into clean subdirectories
    for (const auto &[name, data] : results) {
        // Skip global curve arrays
        if (name == "ssf_q_axis" || name == "ssf_values") continue; 

        // Check if this variable belongs to a specific q-shell
        size_t q_pos = name.find("_q_");
        if (q_pos != std::string::npos) {
            std::string prop_name = name.substr(0, q_pos);
            std::string q_val_str = name.substr(q_pos + 3);

            // Create a subdirectory path for this specific magnitude
            std::string shell_dir = output_dir + "/q_magnitude_" + q_val_str;
            fs::create_directories(shell_dir);

            std::string outpath = shell_dir + "/" + prop_name + ".dat";
            std::ofstream out(outpath);
            if (!out.is_open()) continue;

            // Frequency spectrum format handler (DSF, current frequencies)
            if (prop_name == "dsf" || prop_name == "cl_w" || prop_name == "ct_w") {
                size_t sym_length = (data.size() - 1) * 2; 
                double d_omega = (2.0 * M_PI) / (static_cast<double>(sym_length) * dt);
                out << "# omega spectral_density\n";
                for (size_t i = 0; i < data.size(); ++i) {
                    out << i * d_omega << " " << data[i] << "\n";
                }
            } 
            else {
                // Standard time correlation format handler (ISF, current decays)
                out << "# time correlation_value\n";
                for (size_t i = 0; i < data.size(); ++i) {
                    out << i * dt << " " << data[i] << "\n";
                }
            }
        } 
        else {
            // Fallback pathway for completely q-independent properties (like VACF)
            std::string outpath = output_dir + "/" + name + ".dat";
            std::ofstream out(outpath);
            if (!out.is_open()) continue;

            out << "# time correlation_value\n";
            for (size_t i = 0; i < data.size(); ++i) {
                out << i * dt << " " << data[i] << "\n";
            }
        }
    }
}

void CorrAccumulator::plot(const std::string &output_dir,
                           const std::string &format) const {
    namespace fs = std::filesystem;
    fs::path abs_output_dir = fs::absolute(output_dir);
    const std::string script_path = std::string(SOURCE_DIR) + "/src/fp_plot.py";

    std::string command = "python3 " + script_path + " corr " +
                        abs_output_dir.string() + " " + format;

    int result = std::system(command.c_str());
    if (result != 0) {
        std::cerr << "Warning: Failed to run correlation plotting script.\n";
    }
}

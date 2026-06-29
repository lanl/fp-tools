#pragma once

#include "base.h"
#include "config.h"

#include <fftw3.h>
#include <map>
#include <string>
#include <vector>
#include <unordered_map>

/**
 * @class CorrAccumulator
 * @brief Computes built-in and custom correlation functions via FFTW.
 */
class CorrAccumulator {
public:
    // Constructor
    CorrAccumulator();
    ~CorrAccumulator();

    void add_builtin(const std::string &name, const std::string &input_dir);
    void add_custom(const std::string &name,
                  const std::vector<std::string> &files,
                  const std::vector<int> &columns,
                  const std::string &input_dir);

    void accumulate(const std::string &zero_pad = "next_pow2",
                  bool subtract_mean = true, bool normalize = true);

    // Optional: write g(t) data to files
    void write(const std::string &output_dir, double dt) const;

    // Optional: plot using Python
    void plot(const std::string &output_dir, const std::string &format) const;

    std::unordered_map<std::string, std::vector<double>> results;

    void store_result(const std::string &name, const std::vector<double> &data) {
        results[name] = data;
    }
    void merge(const CorrAccumulator &other);

private:
    // Correlation computation
    std::vector<double> correlate(const std::vector<double> &a,
                                const std::vector<double> &b,
                                const std::string &pad_method,
                                bool subtract_mean, bool normalize);

    // FFT utilities
    size_t get_padded_size(size_t n, const std::string &method) const;

    // I/O utilities
    std::vector<double> read_column(const std::string &filename,
                                  int column_index);
    void save_result(const std::string &name, const std::vector<double> &g_t);

    // Results
    std::vector<std::string> result_names;
    std::vector<std::vector<double>> g_t_list;

    // FFT workspace (reused)
    fftw_complex *fft_a = nullptr;
    fftw_complex *fft_b = nullptr;
    fftw_complex *ifft_result = nullptr;
    fftw_plan plan_a;
    fftw_plan plan_b;
    fftw_plan plan_ifft;

    // FFTW buffer size
    size_t fft_size = 0;

    // Cleanup
    void cleanup_fftw();
};

CorrAccumulator compute_correlations(const CorrConfig &cfg,
                                     const std::vector<Frame> &frames,
                                     const std::vector<Vector3D>& q_vectors);

CorrAccumulator compute_dsf(const CorrConfig& cfg, 
                            const std::vector<Frame>& frames,
                            const std::vector<Vector3D>& q_vectors);

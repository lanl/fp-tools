#pragma once

#include "base.h"
#include "config.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

class Logger {
public:
    explicit Logger(const std::string& filepath,
                    bool enable_console=true);           // console on/off
    ~Logger();

    void banner(const std::string& version);
    void kv     (const std::string& key,
                 const std::string& val,
                 int pad = 24);                          // “  – key … = value”
    void kv(const std::string& key, double val, int pad);
    void note   (const std::string& msg);                // plain info line
    void newline();
    void raw(const std::string& text);
    void status_line(const std::string& key, const std::string& value);
    void arrow(const std::string& message);

    void section(const std::string& title);

    void rdf_info(double dr, int num_bins, double r_max, std::vector<std::vector<std::string>> pairs);
    void msd_info(double dt, int max_lag, const std::vector<std::vector<std::string>>& groups);
    void corr_info(const std::vector<std::string>& builtins,
                       const std::vector<std::string>& custom_names,
                       const std::string& zero_pad,
                       bool normalize,
                       bool subtract_mean);
    void mfpt_info(double dr, double r_max, double dt);

    void bullet(const std::string& text);

    void msd_summ(const std::vector<std::vector<std::string>>& groups, int num_frames, double dt, double elapsed_seconds);
    void rdf_summ(const std::vector<std::vector<std::string>>& pairs,
                      int num_frames,
                      double dt,
                      double elapsed_seconds,
                      std::optional<double> dr,
                      std::optional<int> num_bins,
                      double r_max);
    void traj_info(const InputConfig& config, const std::vector<Frame>& frames);

    void sys_info(const std::string& host,
                  const std::string& wd,
                  const std::string& ts);
    void list(const std::string& label, const std::vector<std::string>& items);
    void list(const std::string& label, const std::vector<std::vector<std::string>>& list_of_lists); 

    void mfpt_summ(int num_radii, int num_frames, double dr, double r_max, double dt, double elapsed_seconds);

private:
    std::ofstream out;
    bool   to_console;

    void write_both(const std::string& txt, bool bold=false);
};

void print_progress_bar(std::size_t current, std::size_t total, int bar_width = 50, const std::string& prefix = "", const std::string& unit_label = "");

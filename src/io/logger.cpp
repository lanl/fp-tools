#include "logger.h"
#include "base.h"
#include "config.h"
#include <iomanip>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// ---------- ctor / dtor ----------------------------------------------------
Logger::Logger(const std::string &fp, bool enable_console)
    : to_console(enable_console) {
  out.open(fp);
  if (!out.is_open())
    throw std::runtime_error("Cannot open log file: " + fp);
}
Logger::~Logger() {
  if (out.is_open())
    out.close();
}

// ---------- internal helper -----------------------------------------------
void Logger::write_both(const std::string &txt, bool bold) {
  if (to_console) {
    if (bold)
      std::cout << "\033[1m";
    std::cout << txt;
    if (bold)
      std::cout << "\033[0m";
  }
  out << txt;
}

// ---------- public API -----------------------------------------------------
void Logger::banner(const std::string &version) {
  const std::string line(60, '-');
  write_both("\n" + line + "\n");
  write_both("             Fp‑Tools – First Passage Toolkit\n", true);
  write_both("                      version " + version + "\n");
  write_both("              Los Alamos National Lab, 2025\n");
  write_both(line + "\n\n");
}

void Logger::kv(const std::string &key, const std::string &val, int pad) {
  std::ostringstream oss;
  oss << "  - " << std::left << std::setw(pad) << key << "= " << val << "\n";
  write_both(oss.str());
}

// Overload for numeric values
void Logger::kv(const std::string &key, double val, int pad) {
  std::ostringstream oss;
  oss << "  - " << std::left << std::setw(pad) << key << "= " << std::fixed
      << std::setprecision(2) << val << "\n";
  write_both(oss.str());
}

void Logger::note(const std::string &msg) { write_both(msg + "\n"); }

void Logger::newline() { write_both("\n"); }

void Logger::sys_info(const std::string &host, const std::string &wd,
                      const std::string &ts) {

  note("  Start time: " + ts);
  note("  Running on host: " + host);
  note("  Current working directory: " + wd);
}

static void print_colored(const std::string &txt, const char *ansi_color) {
  std::cout << ansi_color << txt << "\033[0m"; // reset after text
}

void Logger::arrow(const std::string &msg) {
  const std::string line = "▶ " + msg + "\n";
  out << line; // logfile (no color)

  if (to_console)
    print_colored(line, "\033[92m");
}

void Logger::bullet(const std::string &text) {
  std::cout << "\t• " << text << "\n";
}

void Logger::section(const std::string &title) {
  const int width = 60; // choose your width (matching banner)
  std::string line(width, '=');

  int padding = (width - static_cast<int>(title.size())) / 2;
  if (padding < 0)
    padding = 0;

  out << line << "\n";
  out << std::string(padding, ' ') << title << "\n";
  out << line << "\n";

  std::cout << line << "\n";
  std::cout << std::string(padding, ' ') << title << "\n";
  std::cout << line << "\n";
}

void Logger::rdf_info(double dr, int num_bins, double r_max,
                      std::vector<std::vector<std::string>> pairs) {
    section("Radial Distribution Function (RDF or g(r))");
    newline();
    kv("r_max", r_max, 12);
    kv("dr", std::to_string(dr), 12);

    std::ostringstream pair_str;
    for (const auto &p : pairs) {
        if (p.size() == 1) {
            pair_str << "[" << p[0] << "]";
        } else {
            pair_str << "[" << p[0] << ", " << p[1] << "] ";
        }
    }
    kv("pairs", pair_str.str(), 12);
    newline();

    arrow("Beginning RDF computation for specified atomic pairs...");
}

void Logger::msd_info(double dt, int max_lag,
                      const std::vector<std::vector<std::string>> &groups) {
  section("Mean Squared Displacement (MSD)");
  newline();
  kv("dt", dt, 9);
  kv("max_lag", max_lag, 9);

  // Format groups into a string like: [2], [3]
  std::ostringstream oss;
  for (size_t i = 0; i < groups.size(); ++i) {
    oss << "[";
    for (size_t j = 0; j < groups[i].size(); ++j) {
      oss << groups[i][j];
      if (j < groups[i].size() - 1)
        oss << ", ";
    }
    oss << "]";
    if (i < groups.size() - 1)
      oss << ", ";
  }
  kv("Group(s)", oss.str(), 9);

  newline();
  arrow("Beginning MSD computation...");
}

void Logger::mfpt_info(double dr, double r_max, double dt) {
  section("Mean Free Passage Time (MFPT) and D(r)");
  newline();

  kv("dt", dt, 6);
  kv("dr", dr, 6);
  kv("r_max", r_max, 6);
  newline();
  arrow("Beginning MFPT/D(r) calculation...");
}

void Logger::corr_info(const std::vector<std::string> &builtins,
                       const std::vector<std::string> &custom_names,
                       const std::string &zero_pad, bool normalize,
                       bool subtract_mean) {
  section("Correlation Function(s)");
  newline();

  kv("Zero padding", zero_pad);
  kv("Normalize", normalize ? "true" : "false");
  kv("Subtract mean", subtract_mean ? "true" : "false");
  if (!builtins.empty()) {
    kv("Built-in functions:", " ");
    for (const auto &name : builtins) {
      bullet(name);
    }
    newline();
  }

  if (!custom_names.empty()) {
    kv("Custom correlation pairs:", " ");
    for (const auto &name : custom_names) {
      bullet(name);
    }
    newline();
  }

  arrow("Beginning correlation computation...");
}

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

void print_progress_bar(std::size_t current, std::size_t total, int bar_width,
                        const std::string &prefix,
                        const std::string &unit_label) {
  if (total == 0)
    return;

  static const char *fractional_blocks[] = {"",  "▏", "▎", "▍", "▌",
                                            "▋", "▊", "▉", "█"};
  static const char spinner_chars[] = {'|', '/', '-', '\\'};
  static int spinner_index = 0;

  float progress = static_cast<float>(current) / total;
  std::size_t full_blocks = static_cast<std::size_t>(progress * bar_width);
  float fractional = (progress * bar_width) - full_blocks;
  std::size_t fractional_index = static_cast<std::size_t>(fractional * 8);

  std::string bar;

  // Add full blocks (each visually 1 char)
  for (std::size_t i = 0; i < full_blocks; ++i)
    bar += "█";

  // Add fractional block (also visually 1 char or empty)
  if (fractional_index > 0)
    bar += fractional_blocks[fractional_index];

  // Calculate visual length of the bar so far
  std::size_t bar_length = full_blocks + (fractional_index > 0 ? 1 : 0);

  // Pad with spaces to ensure fixed width
  while (bar_length < static_cast<std::size_t>(bar_width)) {
    bar += " ";
    ++bar_length;
  }

  // Green color codes
  const std::string green = "\033[1;32m";
  const std::string reset = "\033[0m";

  // Build and print line
  std::cout << "\r\033[K";
  std::cout << green << "▶ " << prefix << " |" << bar << "| "
            << std::setw(3) << static_cast<int>(progress * 100.0f) << "% ("
            << current << "/" << total << " " << unit_label << ") " << green
            << reset << std::flush;

  if (current == total) std::cout << std::endl;

}

void Logger::msd_summ(const std::vector<std::vector<std::string>> &groups,
                      int num_frames, double dt, double elapsed_seconds) {
  std::ostringstream oss;

  oss << "[ MSD ] Finished in " << std::fixed << std::setprecision(2)
      << elapsed_seconds << " seconds → " << groups.size() << " group"
      << (groups.size() == 1 ? "" : "s") << ", " << num_frames
      << " frames, dt = " << dt << ".\n";

  std::cout << oss.str();
}

void Logger::rdf_summ(const std::vector<std::vector<std::string>> &pairs,
                      int num_frames, double dt, double elapsed_seconds,
                      std::optional<double> dr, std::optional<int> num_bins,
                      double r_max) {
  std::ostringstream oss;

  oss << "[ RDF ] Finished in " << std::fixed << std::setprecision(2)
      << elapsed_seconds << " seconds → " << pairs.size() << " pair"
      << (pairs.size() == 1 ? "" : "s") << ", " << num_frames << " frames"
      << ", dt = " << dt;

  if (dr.has_value()) {
    oss << ", dr = " << dr.value();
  } else if (num_bins.has_value()) {
    oss << ", " << num_bins.value() << " bins";
  }

  oss << ", r_max = " << r_max << ".\n";
  std::cout << oss.str();
}

void Logger::mfpt_summ(int num_radii, int num_frames, double dr, double r_max,
                       double dt, double elapsed_seconds) {
  std::ostringstream oss;

  oss << "[ MFPT ] Finished in " << std::fixed << std::setprecision(2)
      << elapsed_seconds << " seconds → " << num_radii << " radii (dr = " << dr
      << ", r_max = " << r_max << "), " << num_frames << " frames, dt = " << dt
      << ".\n";

  std::cout << oss.str();
}

void Logger::traj_info(const InputConfig &config,
                       const std::vector<Frame> &frames) {
  if (frames.empty()) {
    note("No frames available to summarize.");
    return;
  }

  const auto &bb = frames[0].box_bounds;

  std::ostringstream xbounds, ybounds, zbounds;
  xbounds << "[" << std::fixed << std::setprecision(2) << bb[0][0] << ", "
          << std::fixed << std::setprecision(2) << bb[0][1] << "]";
  ybounds << "[" << std::fixed << std::setprecision(2) << bb[1][0] << ", "
          << std::fixed << std::setprecision(2) << bb[1][1] << "]";
  zbounds << "[" << std::fixed << std::setprecision(2) << bb[2][0] << ", "
          << std::fixed << std::setprecision(2) << bb[2][1] << "]";

  double Lx = bb[0][1] - bb[0][0];
  double Ly = bb[1][1] - bb[1][0];
  double Lz = bb[2][1] - bb[2][0];
  std::ostringstream Lbox;
  Lbox << std::fixed << std::setprecision(2) << Lx << " x " << std::fixed
       << std::setprecision(2) << Ly << " x " << std::fixed
       << std::setprecision(2) << Lz;

  // Collect atom types
  std::set<std::string> types;
  for (const auto &atom : frames[0].atoms)
    types.insert(atom.type);

  std::ostringstream type_oss;
  size_t count = 0;
  for (const auto &t : types) {
    if (count++ > 0)
      type_oss << ", ";
    type_oss << t;
  }

  // Print aligned key–value pairs
  kv("Trajectory file ", config.file, 22);
  kv("Parsed frames ", static_cast<int>(frames.size()), 22);
  kv("Frame interval ", config.frame_interval, 22);
  kv("Box bounds (x) ", xbounds.str(), 22);
  kv("Box bounds (y) ", ybounds.str(), 22);
  kv("Box bounds (z) ", zbounds.str(), 22);
  kv("Box size (Lx, Ly, Lz) ", Lbox.str(), 22);
  kv("Atom types present ", type_oss.str(), 22);
  kv("Atoms per frame ", static_cast<int>(frames[0].atomnum), 22);
  newline();
}

void Logger::list(const std::string &label,
                  const std::vector<std::string> &items) {
  note("  - " + label + ": ");
  for (const auto &item : items) {
    note("  - " + item);
  }
}

void Logger::list(const std::string &label,
                  const std::vector<std::vector<std::string>> &list_of_lists) {
  note("  - " + label + ": ");
  for (const auto &group : list_of_lists) {
    std::ostringstream oss;
    oss << "    - [";
    for (size_t i = 0; i < group.size(); ++i) {
      oss << group[i];
      if (i < group.size() - 1)
        oss << ", ";
    }
    oss << "]";
    note(oss.str());
  }
}

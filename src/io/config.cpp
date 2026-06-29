#include "config.h"
#include "config_utils.h"
#include "toml.hpp"
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace tomlutil;

Config parse_config(const std::string &filename) {
    Config config{};
    toml::table root = toml::parse_file(filename);

    // Input
    if (const auto *input = get_table(root, "input")) {
        config.input.file = get_or(*input, "file", std::string{""});
        config.input.vel_file = get_or(*input, "vel_file", std::string{""});
        config.input.coord_type =
            get_or(*input, "coord_type", std::string{"unwrapped"});
        config.input.start_frame = get_or(*input, "start_frame", 0);
        config.input.end_frame = get_or(*input, "end_frame", -1);
        config.input.frame_interval = get_or(*input, "frame_interval", 1);
    }

    // Output
    if (const auto *output = get_table(root, "output")) {
        config.output.path = get_or(*output, "path", std::string{"output"});
    }

    // System
    if (const auto *properties = get_table(root, "properties")) {
        config.properties.type = get_or(*properties, "type", std::string{""});
        config.properties.temp = get_or(*properties, "temp", 0);
        config.properties.density = get_or(*properties, "density", 10);
        config.properties.dt = get_or(*properties, "dt", 0.01);

        if (const auto *box_arr = get_array(*properties, "box_lengths")) {
            size_t idx{0};
            for (const auto& val : *box_arr) {
                if (idx >= 3) break;
                if (const auto v = val.value<double>()) {
                    config.properties.box_lengths[idx++] = *v;
                }
            }
        }
    }

    // Analysis
    if (const auto *analysis = get_table(root, "analysis")) {

        if (const auto *msd = get_table(*analysis, "msd")) {
            config.analysis.msd.enabled = get_or(*msd, "enabled", false);
            config.analysis.msd.max_lag = get_or(*msd, "max_lag", -1);
            config.analysis.msd.plot_format =
                get_or(*msd, "plot_format", std::string{""});

            if (const auto *groups = get_array(*msd, "groups")) {
                for (const auto &group : *groups) {
                    std::vector<std::string> entry;
                    if (const auto *inner = group.as_array()) {
                        for (const auto &item : *inner) {
                            if (item.is_string()) {
                                entry.push_back(*item.value<std::string>());
                            }
                        }
                        if (!entry.empty()) {
                            config.analysis.msd.groups.push_back(entry);
                        } else {
                            std::cerr << "Warning: MSD group must contain at least one atom "
                                         "type. Skipping invalid entry.\n";
                        }
                    }
                }
            }

            if (config.analysis.msd.groups.empty()) {
                config.analysis.msd.groups.push_back({"all"});
            }

            if (const auto *curves_array = get_array(*msd, "curves")) {
                for (const auto &plot_entry : *curves_array) {
                    std::vector<Curve> plot_curves;
                    if (const auto *inner_array = plot_entry.as_array()) {
                        for (const auto &curve_entry : *inner_array) {
                            if (const auto *curve_table = curve_entry.as_table()) {
                                Curve curve;

                                // Parse 'types' field
                                if (const auto *types_array = get_array(*curve_table, "types")) {
                                    for (const auto &type_val : *types_array) {
                                        if (type_val.is_string()) {
                                            curve.types.push_back(*type_val.value<std::string>());
                                        }
                                    }
                                }

                                // Parse 'label' field
                                curve.label = get_or(*curve_table, "label", std::string{""});

                                if (!curve.types.empty()) {
                                    plot_curves.push_back(std::move(curve));
                                } else {
                                    std::cerr << "Warning: MSD curve must have at least one "
                                                 "type. Skipping.\n";
                                }
                            }
                        }
                    }

                    if (!plot_curves.empty()) {
                        config.analysis.msd.curves.push_back(std::move(plot_curves));
                    } else {
                        std::cerr << "Warning: MSD plot group must contain at least one "
                                     "curve. Skipping.\n";
                    }
                }
            }
        }
    

        if (const auto *rdf = get_table(*analysis, "rdf")) {
            config.analysis.rdf.enabled = get_or(*rdf, "enabled", false);
            config.analysis.rdf.dr = get_or(*rdf, "dr", -1.0);
            config.analysis.rdf.r_max = get_or(*rdf, "r_max", -1.0);
            config.analysis.rdf.num_bins = get_or(*rdf, "num_bins", -1);
            config.analysis.rdf.plot_format =
                get_or(*rdf, "plot_format", std::string{""});

            if (const auto *combos = get_array(*rdf, "pairs")) {
                for (const auto &pair : *combos) {
                    std::vector<std::string> combo;
                    if (const auto *inner = pair.as_array()) {
                        for (const auto &item : *inner) {
                            if (item.is_string()) 
                                combo.push_back(*item.value<std::string>());
                        }
                        if (combo.size() == 2) {
                            config.analysis.rdf.pairs.push_back(combo);
                        } else {
                            std::cerr << "Warning: RDF combination must contain exactly two "
                                         "atom types. Skipping invalid entry.\n";
                        }
                    }
                }
            }

            if (config.analysis.rdf.pairs.empty()) {
                config.analysis.rdf.pairs.push_back({"all"});
            }

            if (const auto *curves_array = get_array(*rdf, "curves")) {
                for (const auto &plot_entry : *curves_array) {
                    std::vector<Curve> plot_curves;
                    if (const auto *inner_array = plot_entry.as_array()) {
                        for (const auto &curve_entry : *inner_array) {
                            if (const auto *curve_table = curve_entry.as_table()) {
                                Curve curve;
                                if (const auto *types_array = get_array(*curve_table, "types")) {
                                    for (const auto &type_val : *types_array) {
                                        if (type_val.is_string()) {
                                            curve.types.push_back(*type_val.value<std::string>());
                                        }
                                    }
                                }
                                curve.label = get_or(*curve_table, "label", std::string{""});

                                if (!curve.types.empty()) {
                                    plot_curves.push_back(std::move(curve));
                                }
                            }
                        }
                    }
                    if (!plot_curves.empty()) {
                        config.analysis.rdf.curves.push_back(std::move(plot_curves));
                    }
                }
            }
        }

        if (const auto *mfpt = get_table(*analysis, "mfpt")) {
            config.analysis.mfpt.enabled = get_or(*mfpt, "enabled", false);
            config.analysis.mfpt.dt = get_or(*mfpt, "dt", 0.001);
            config.analysis.mfpt.dr = get_or(*mfpt, "dr", 0.01);
            config.analysis.mfpt.r_max = get_or(*mfpt, "r_max", 1.0);
            config.analysis.mfpt.plot_format =
                get_or(*mfpt, "plot_format", std::string{""});

            if (const auto *spec_node = get_array(*mfpt, "species")) {
                for (const auto &val : *spec_node) {
                    if (const auto str = val.value<std::string>()) {
                        config.analysis.mfpt.species.push_back(*str);
                    }
                }
            }
        }

        if (const auto *corr = get_table(*analysis, "correlations")) {
            config.analysis.corr.enabled = get_or(*corr, "enabled", false);
            config.analysis.corr.dt = get_or(*corr, "dt", 0.0);
            config.analysis.corr.normalize = get_or(*corr, "normalize", false);
            config.analysis.corr.zero_pad = get_or(*corr, "zero_pad", std::string{});
            config.analysis.corr.plot_format =
                get_or(*corr, "plot_format", std::string{});
            config.analysis.corr.subtract_mean =
                get_or(*corr, "subtract_mean", false);

            config.analysis.corr.n_list = get_or(*corr, "n_list", std::string{""});
        
            if (const auto *n_node = get_array(*corr, "n")) {
                size_t idx{0};
                for (const auto &val : *n_node) {
                    if (idx >= 3) break;
                    if (const auto d = val.value<double>()) {
                        config.analysis.corr.n[idx++] = *d;
                    }
                }
            }

            if (const auto *builtin_node = get_array(*corr, "builtin")) {
                for (const auto &val : *builtin_node) {
                    if (const auto str = val.value<std::string>()) {
                        config.analysis.corr.builtins.push_back(*str);
                    }
                }
            }

            if (const auto *customs = get_array(*corr, "custom")) {
                for (const auto &entry : *customs) {
                    if (!entry.is_table())
                        continue;
                    const auto &t = *entry.as_table();

                    CustomCorr c;
                    c.name = get_or(t, "name", std::string{});

                    if (const auto *files_arr = get_array(t, "files")) {
                        for (const auto &val : *files_arr) {
                            if (val.is_string()) {
                                c.files.push_back(*val.value<std::string>());
                            }
                        }
                    }
                    if (const auto *col_node = get_array(t, "columns")) {
                        if (const auto *arr = col_node->as_array()) {
                            for (const auto &val : *arr) {
                                if (const auto i = val.value<int64_t>()) {
                                    c.columns.push_back(static_cast<int>(*i));
                                }
                            }
                        }
                    }
                    config.analysis.corr.custom.push_back(c);
                }
            }
        }
    }
    return config;
}

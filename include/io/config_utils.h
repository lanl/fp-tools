#pragma once

#include "toml.hpp"

#include <string>

namespace tomlutil {

template <typename T>
T get_or(const toml::table& tbl, std::string_view key, const T& def) {
  if (auto node = tbl.get(key)) {
    if (auto val = node->value<T>())
      return *val;
  }
  return def;
}

inline const toml::table* get_table(const toml::table& parent, std::string_view key) {
  if (auto node = parent.get(key)) {
    return node->as_table();
  }
  return nullptr;
}

inline const toml::array* get_array(const toml::table& parent, std::string_view key) {
  if (auto node = parent.get(key)) {
    return node->as_array();
  }
  return nullptr;
}

}  // namespace tomlutil


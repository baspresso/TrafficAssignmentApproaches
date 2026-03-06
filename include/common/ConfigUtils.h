#ifndef CONFIG_UTILS_H
#define CONFIG_UTILS_H

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace TrafficAssignment::ConfigUtils {

namespace fs = std::filesystem;

struct CliOptions {
  std::unordered_map<std::string, std::string> values;
  bool help_requested = false;
};

inline std::string ToLowerCopy(std::string text) {
  std::transform(
    text.begin(),
    text.end(),
    text.begin(),
    [](unsigned char c) { return static_cast<char>(std::tolower(c)); }
  );
  return text;
}

inline std::string NormalizeKey(std::string key) {
  key = ToLowerCopy(std::move(key));
  std::replace(key.begin(), key.end(), '-', '_');
  return key;
}

inline std::string TrimCopy(const std::string& text) {
  const auto first = text.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) {
    return "";
  }
  const auto last = text.find_last_not_of(" \t\r\n");
  return text.substr(first, last - first + 1);
}

inline bool ParseBool(const std::string& raw_value, const std::string& field_name) {
  const std::string value = ToLowerCopy(raw_value);
  if (value == "1" || value == "true" || value == "yes" || value == "on") {
    return true;
  }
  if (value == "0" || value == "false" || value == "no" || value == "off") {
    return false;
  }
  throw std::runtime_error(
    "Invalid boolean value for '" + field_name + "': '" + raw_value + "'"
  );
}

template <typename T>
T ParseNumber(const std::string& raw_value, const std::string& field_name) {
  std::stringstream stream(raw_value);
  T value {};
  stream >> value;
  if (stream.fail() || !stream.eof()) {
    throw std::runtime_error(
      "Invalid numeric value for '" + field_name + "': '" + raw_value + "'"
    );
  }
  return value;
}

inline std::optional<std::string> GetEnvValue(const char* name) {
  const char* value = std::getenv(name);
  if (value == nullptr || value[0] == '\0') {
    return std::nullopt;
  }
  return std::string(value);
}

inline CliOptions ParseCliOptions(int argc, char** argv) {
  CliOptions cli;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      cli.help_requested = true;
      continue;
    }
    if (arg.rfind("--", 0) != 0) {
      throw std::runtime_error("Unsupported argument format: '" + arg + "'");
    }

    arg = arg.substr(2);
    std::string key;
    std::string value;
    const auto eq_pos = arg.find('=');
    if (eq_pos != std::string::npos) {
      key = arg.substr(0, eq_pos);
      value = arg.substr(eq_pos + 1);
    } else {
      key = arg;
      if ((i + 1) < argc && std::string(argv[i + 1]).rfind("--", 0) != 0) {
        value = argv[++i];
      } else {
        value = "true";
      }
    }

    cli.values[NormalizeKey(key)] = value;
  }
  return cli;
}

inline fs::path FindProjectRoot() {
  fs::path current = fs::current_path();
  while (true) {
    if (fs::exists(current / "CMakeLists.txt") &&
        fs::exists(current / "data" / "TransportationNetworks")) {
      return current;
    }
    if (!current.has_parent_path() || current == current.parent_path()) {
      return fs::current_path();
    }
    current = current.parent_path();
  }
}

inline fs::path ResolvePath(const fs::path& raw_path, const fs::path& project_root) {
  if (raw_path.empty()) {
    return raw_path;
  }
  if (raw_path.is_absolute()) {
    return raw_path;
  }
  if (fs::exists(project_root / raw_path)) {
    return project_root / raw_path;
  }
  return fs::absolute(raw_path);
}

}  // namespace TrafficAssignment::ConfigUtils

#endif

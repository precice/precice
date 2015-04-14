// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_UTILS_PUBLISHER_HPP_
#define PRECICE_UTILS_PUBLISHER_HPP_

#include <string>

namespace precice {
namespace utils {
class Publisher {
public:
  struct ScopedPrefix {
    ScopedPrefix(std::string const& prefix);

    ~ScopedPrefix();

  private:
    std::string _prefix;
  };

  struct ScopedPublication {
    ScopedPublication(std::string const& filePath, std::string const& data);

    ~ScopedPublication();

  private:
    std::string _filePath;
  };

public:
  static void read(std::string const& filePath, std::string& data);

  static void write(std::string const& filePath, std::string const& data);

  static void remove(std::string const& filePath);

  static void setPrefix(std::string const& prefix);

  static std::string const& prefix();

private:
  static std::string _prefix;
};
} // namespace utils
} // namespace precice

#endif /* PRECICE_UTILS_PUBLISHER_HPP_ */

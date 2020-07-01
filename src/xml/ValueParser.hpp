#pragma once

#include <Eigen/Core>
#include <string>
#include "logging/Logger.hpp"
#include "utils/String.hpp"

namespace precice {
namespace xml {

static logging::Logger _log{"xml::ValueParser"};

inline void readValueSpecific(const std::string &rawValue, double &value)
{
  try {
    if (rawValue.find('/') != std::string::npos) {
      std::string left  = rawValue.substr(0, rawValue.find('/'));
      std::string right = rawValue.substr(rawValue.find('/') + 1, rawValue.size() - rawValue.find('/') - 1);

      value = std::stod(left) / std::stod(right);
    } else {
      value = std::stod(rawValue);
    }
  } catch (...) {
    PRECICE_ERROR("String to Double error");
  }
}

inline void readValueSpecific(const std::string &rawValue, int &value)
{
  try {
    value = std::stoi(rawValue);
  } catch (...) {
    PRECICE_ERROR("String to Int error");
  }
}

inline void readValueSpecific(const std::string &rawValue, std::string &value)
{
  value = rawValue;
}

inline void readValueSpecific(const std::string &rawValue, bool &value)
{
  value = precice::utils::convertStringToBool(rawValue);
}

inline void readValueSpecific(const std::string &rawValue, Eigen::VectorXd &value)
{
  Eigen::VectorXd vec;

  std::string valueString(rawValue);
  bool        componentsLeft = true;
  int         i              = 0;
  while (componentsLeft) {
    std::string tmp1(rawValue);
    // erase entries before i-th entry
    for (int j = 0; j < i; j++) {
      if (tmp1.find(';') != std::string::npos) {
        tmp1.erase(0, tmp1.find(';') + 1);
      } else {
        componentsLeft = false;
      }
    }
    // if we are not in the last vector component...
    if (tmp1.find(';') != std::string::npos) {
      // ..., erase entries after i-th entry
      tmp1.erase(tmp1.find(';'), tmp1.size());
    }

    if (componentsLeft) {
      vec.conservativeResize(vec.rows() + 1);
      vec(vec.rows() - 1) = std::stod(tmp1);
    }
    i++;
  }

  value = vec;
}

} // namespace xml
} // namespace precice

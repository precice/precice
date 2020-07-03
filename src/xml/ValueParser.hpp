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

void readValueSpecific(const std::string &rawValue, Eigen::VectorXd &value);

} // namespace xml
} // namespace precice

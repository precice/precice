#pragma once

#include <Eigen/Core>
#include <string>
#include "utils/String.hpp"

namespace precice {
namespace xml {

void readValueSpecific(const std::string &rawValue, double &value);

void readValueSpecific(const std::string &rawValue, int &value);

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

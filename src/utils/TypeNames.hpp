/*
 * Returns pretty type names for various types.
 * Used for generating the XML documentation.
 * They are inlined, so that they don't need to be defined in a seperate cpp file.
 */

#pragma once

#include <Eigen/Core>
#include <string>

namespace precice {
namespace utils {

inline std::string getTypeName(const double &var)
{
  return "float";
}

inline std::string getTypeName(const std::string &var)
{
  return "string";
}

inline std::string getTypeName(const bool &var)
{
  return "boolean";
}

inline std::string getTypeName(const int &var)
{
  return "integer";
}

inline std::string getTypeName(Eigen::VectorXd const &var)
{
  return "vector";
}

} // namespace utils
} // namespace precice

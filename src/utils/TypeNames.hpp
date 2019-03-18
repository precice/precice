/*
 * Returns pretty type names for various types.
 * Used for generating the XML documentation.
 * They are inlined, so that they don't need to be defined in a seperate cpp file.
 */

#pragma once

#include <string>
#include <Eigen/Core>

namespace precice {
namespace utils {

inline std::string getTypeName(const double&)
{
  return "float";
}

inline std::string getTypeName(const std::string&)
{
  return "string";
}

inline std::string getTypeName(const bool&)
{
  return "boolean";
}

inline std::string getTypeName(const int&)
{
  return "integer";
}

inline std::string getTypeName(Eigen::VectorXd const &)
{
  return "vector";
}

}}

/*
 * Returns pretty type names for various types.
 * Used for generating the XML documentation.
 * They are inlined, so that they don't need to be defined in a separate cpp file.
 */

#pragma once

#include <Eigen/Core>
#include <string>

namespace precice::utils {

<<<<<<< HEAD
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
=======
template <typename T>
constexpr std::string_view getTypeName(const T &)
{
  using U = std::remove_cv_t<std::remove_reference_t<T>>;

  if constexpr (std::is_same_v<U, double>)
    return "float";
  else if constexpr (std::is_same_v<U, std::string>)
    return "string";
  else if constexpr (std::is_same_v<U, bool>)
    return "boolean";
  else if constexpr (std::is_same_v<U, int>)
    return "integer";
  else if constexpr (std::is_same_v<U, Eigen::VectorXd>)
    return "vector";
  else
    static_assert(!sizeof(U), "Unsupported type passed to getTypeName()");
>>>>>>> 5afb891e (Refactor typename getter logic)
}

} // namespace precice::utils

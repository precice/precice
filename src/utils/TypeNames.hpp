/*
 * Returns pretty type names for various types.
 * Used for generating the XML documentation.
 */

#pragma once

#include <Eigen/Core>
#include <string>

namespace precice::utils {

template <typename T>
constexpr std::string_view getTypeName(const T &)
{
  // if we move to cpp20 can switch this with std::remove_cvref_t
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
    static_assert(false, "Unsupported type passed to getTypeName()");
}

} // namespace precice::utils

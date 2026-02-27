/*
 * Returns pretty type names for various types.
 * Used for generating the XML documentation.
 *
 * Note that all are templated so we do not need
 * a separate .cpp file
 *
 * Due to constexpr and string_view, it allows
 * for plausible compile time eval and no allocs
 */

#pragma once

#include <Eigen/Core>
#include <string_view>

namespace precice::utils {

// primary template declaration
// (other explicit specifications will follow this)
template <typename T>
constexpr std::string_view getTypeName()
{
  // in the future as the project moves to cpp20
  // this can be replaced by concepts, but this
  // check is fine for now.
  static_assert(sizeof(T) == 0,
                "Unsupported type passed to getTypeName()");
  return {};
}

// value based forwarder
template <typename T>
constexpr std::string_view getTypeName(const T &)
{
  return getTypeName<T>();
}

template <>
constexpr std::string_view getTypeName<double>()
{
  return "float";
}

template <>
constexpr std::string_view getTypeName<std::string>()
{
  return "string";
}

template <>
constexpr std::string_view getTypeName<bool>()
{
  return "boolean";
}

template <>
constexpr std::string_view getTypeName<int>()
{
  return "integer";
}

template <>
constexpr std::string_view getTypeName<Eigen::VectorXd>()
{
  return "vector";
}

} // namespace precice::utils

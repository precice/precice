#pragma once

#include <Eigen/Core>
#include <cstdlib>

namespace precice::math {

constexpr double NUMERICAL_ZERO_DIFFERENCE = 1.0e-14;

/// Compares two Eigen::MatrixBase for equality up to tolerance
template <class DerivedA, class DerivedB>
constexpr bool equals(const Eigen::MatrixBase<DerivedA> &A,
                      const Eigen::MatrixBase<DerivedB> &B,
                      double                             tolerance = NUMERICAL_ZERO_DIFFERENCE)
{
  return A.isApprox(B, tolerance);
}

/// Compares two scalar (arithmetic) types
inline bool equals(double a, double b, double tolerance = NUMERICAL_ZERO_DIFFERENCE) noexcept
{
  auto d     = std::abs(a - b);
  auto scale = std::max(1.0, std::max(std::abs(a), std::abs(b)));
  return d <= tolerance * scale;
}

inline bool isZero(double a, double tolerance = NUMERICAL_ZERO_DIFFERENCE) noexcept
{
  return std::abs(a) < tolerance;
}

inline bool greater(double a, double b, double tolerance = NUMERICAL_ZERO_DIFFERENCE) noexcept
{
  return a > b + tolerance;
}

inline bool greaterEquals(double a, double b, double tolerance = NUMERICAL_ZERO_DIFFERENCE) noexcept
{
  return (a > b) || ::precice::math::equals(a, b, tolerance);
}

inline bool smaller(double a, double b, double tolerance = NUMERICAL_ZERO_DIFFERENCE) noexcept
{
  return a + tolerance < b;
}

inline bool smallerEquals(double a, double b, double tolerance = NUMERICAL_ZERO_DIFFERENCE) noexcept
{
  return (a < b) || ::precice::math::equals(a, b, tolerance);
}

} // namespace precice::math

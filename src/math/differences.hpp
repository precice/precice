#pragma once

#include <Eigen/Core>

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
template <class Scalar>
typename std::enable_if<std::is_arithmetic<Scalar>::value, bool>::type equals(const Scalar a, const Scalar b, const Scalar tolerance = NUMERICAL_ZERO_DIFFERENCE)
{
  auto d = std::abs(a - b);
  // Handle differences close to 0
  if (d < tolerance) {
    return true;
  }
  // Handles larger difference
  return (d <= tolerance * std::max(std::abs(a), std::abs(b)));
}

template <class Scalar>
typename std::enable_if<std::is_arithmetic<Scalar>::value, bool>::type isZero(const Scalar a, const Scalar tolerance = NUMERICAL_ZERO_DIFFERENCE)
{
  return std::abs(a) < tolerance;
}

template <class Scalar>
typename std::enable_if<std::is_arithmetic<Scalar>::value, bool>::type greater(Scalar A, Scalar B, Scalar tolerance = NUMERICAL_ZERO_DIFFERENCE)
{
  return A > B + tolerance;
}

template <class Scalar>
typename std::enable_if<std::is_arithmetic<Scalar>::value, bool>::type greaterEquals(Scalar A, Scalar B, Scalar tolerance = NUMERICAL_ZERO_DIFFERENCE)
{
  return A + tolerance >= B;
}

template <class Scalar>
typename std::enable_if<std::is_arithmetic<Scalar>::value, bool>::type smaller(Scalar A, Scalar B, Scalar tolerance = NUMERICAL_ZERO_DIFFERENCE)
{
  return A + tolerance < B;
}

template <class Scalar>
typename std::enable_if<std::is_arithmetic<Scalar>::value, bool>::type smallerEquals(Scalar A, Scalar B, Scalar tolerance = NUMERICAL_ZERO_DIFFERENCE)
{
  return A <= B + tolerance;
}

} // namespace precice::math

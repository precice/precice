#pragma once

#include <Eigen/Core>
#include "utils/assertion.hpp"

namespace precice {
namespace math {

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
  return std::abs(a - b) <= tolerance;
}

template <class DerivedA, class DerivedB>
bool oneGreater(const Eigen::MatrixBase<DerivedA> &A,
                const Eigen::MatrixBase<DerivedB> &B,
                double                             tolerance = math::NUMERICAL_ZERO_DIFFERENCE)
{
  PRECICE_ASSERT(A.rows() == B.rows(), "Matrices with different number of rows can't be compared.");
  PRECICE_ASSERT(A.cols() == B.cols(), "Matrices with different number of cols can't be compared.");

  return ((A - B).array() > tolerance).any();
}

template <class DerivedA, class DerivedB>
bool oneGreaterEquals(const Eigen::MatrixBase<DerivedA> &A,
                      const Eigen::MatrixBase<DerivedB> &B,
                      double                             tolerance = math::NUMERICAL_ZERO_DIFFERENCE)
{
  PRECICE_ASSERT(A.rows() == B.rows(), "Matrices with different number of rows can't be compared.");
  PRECICE_ASSERT(A.cols() == B.cols(), "Matrices with different number of cols can't be compared.");

  return ((A - B).array() >= -tolerance).any();
}

template <class DerivedA, class DerivedB>
bool allGreater(const Eigen::MatrixBase<DerivedA> &A,
                const Eigen::MatrixBase<DerivedB> &B,
                double                             tolerance = math::NUMERICAL_ZERO_DIFFERENCE)
{
  PRECICE_ASSERT(A.rows() == B.rows(), "Matrices with different number of rows can't be compared.");
  PRECICE_ASSERT(A.cols() == B.cols(), "Matrices with different number of cols can't be compared.");

  return ((A - B).array() > tolerance).all();
}

template <class DerivedA, class DerivedB>
bool allGreaterEquals(const Eigen::MatrixBase<DerivedA> &A,
                      const Eigen::MatrixBase<DerivedB> &B,
                      double                             tolerance = math::NUMERICAL_ZERO_DIFFERENCE)
{
  PRECICE_ASSERT(A.rows() == B.rows(), "Matrices with different number of rows can't be compared.");
  PRECICE_ASSERT(A.cols() == B.cols(), "Matrices with different number of cols can't be compared.");

  return ((A - B).array() >= tolerance).all();
}

template <class Scalar>
typename std::enable_if<std::is_arithmetic<Scalar>::value, bool>::type greater(Scalar A, Scalar B, Scalar tolerance = NUMERICAL_ZERO_DIFFERENCE)
{
  return A - B > tolerance;
}

template <class Scalar>
typename std::enable_if<std::is_arithmetic<Scalar>::value, bool>::type greaterEquals(Scalar A, Scalar B, Scalar tolerance = NUMERICAL_ZERO_DIFFERENCE)
{
  return A - B >= -tolerance;
}

template <class Scalar>
typename std::enable_if<std::is_arithmetic<Scalar>::value, bool>::type smaller(Scalar A, Scalar B, Scalar tolerance = NUMERICAL_ZERO_DIFFERENCE)
{
  return A - B < -tolerance;
}

template <class Scalar>
typename std::enable_if<std::is_arithmetic<Scalar>::value, bool>::type smallerEquals(Scalar A, Scalar B, Scalar tolerance = NUMERICAL_ZERO_DIFFERENCE)
{
  return A - B <= tolerance;
}

} // namespace math
} // namespace precice

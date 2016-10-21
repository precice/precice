#pragma once

#include <Eigen/Dense>
#include "tarch/la/DynamicVector.h"
#include "tarch/la/Vector.h"

namespace precice {
namespace math {

constexpr double NUMERICAL_ZERO_DIFFERENCE = 1.0e-14;

/// Compares two tarch::la::DynamicVectors for equality up to tolerance
template<class lScalar, class rScalar>
constexpr bool equals (const tarch::la::DynamicVector<lScalar>& lVector,
                       const tarch::la::DynamicVector<rScalar>& rVector,
                       double tolerance = NUMERICAL_ZERO_DIFFERENCE)
{
  return tarch::la::equals(lVector, rVector, tolerance);
}


/// Compares two tarch::la::Vectors for equality up to tolerance
template<int size, class lScalar, class rScalar>
constexpr bool equals (const tarch::la::Vector<size, lScalar>& lVector,
                       const tarch::la::Vector<size, rScalar>& rVector,
                       double tolerance = NUMERICAL_ZERO_DIFFERENCE)
{
  return tarch::la::equals(lVector, rVector, tolerance);
}

/// Compares a tarch::la::Vector and an Eigen vector for equality
template<int size, class Scalar, class Derived>
constexpr bool equals (const tarch::la::Vector<size, Scalar>& A,
                       const Eigen::MatrixBase<Derived>& B,
                       double tolerance = NUMERICAL_ZERO_DIFFERENCE)
{
  return tarch::la::equals(A,
                           static_cast<tarch::la::Vector<size, Scalar>>(B),
                           tolerance);
}

/// Compares an Eigen vector and a tarch::la::Vector for equality
template<int size, class Scalar, class Derived>
constexpr bool equals (const Eigen::MatrixBase<Derived>& A,
                       const tarch::la::Vector<size, Scalar>& B,
                       double tolerance = NUMERICAL_ZERO_DIFFERENCE)
{
  return tarch::la::equals(static_cast<tarch::la::Vector<size, Scalar>>(A),
                           B,
                           tolerance);
}

/// Compares an Eigen vector and a tarch::la::DynamicVector for equality
template<class Scalar, class Derived>
constexpr bool equals (const Eigen::MatrixBase<Derived>& A,
                       const tarch::la::DynamicVector<Scalar>& B,
                       double tolerance = NUMERICAL_ZERO_DIFFERENCE)
{
  return tarch::la::equals(static_cast<tarch::la::DynamicVector<Scalar>>(A),
                           B,                           
                           tolerance);
}

template<class Scalar, class Derived>
constexpr bool equals (const tarch::la::DynamicVector<Scalar>& A,
                       const Eigen::MatrixBase<Derived>& B,
                       double tolerance = NUMERICAL_ZERO_DIFFERENCE)
{
  return equals(B, A, tolerance);
}


/// Compares two Eigen::MatrixBase for equality up to tolerance
template <typename DerivedA, typename DerivedB>
constexpr bool equals (const Eigen::MatrixBase<DerivedA>& A,
                       const Eigen::MatrixBase<DerivedB>& B,
                       double tolerance = NUMERICAL_ZERO_DIFFERENCE)
{
  return A.isApprox(B, tolerance);
}

// template <class Scalar, int RowsA, int RowsB, int Cols, int Options, int MaxRowsA, int MaxRowsB, int MaxCols>
// constexpr bool equals (const Eigen::Matrix<Scalar, RowsA, Cols, Options, MaxRowsA, MaxCols>& A,
//                        const Eigen::Matrix<Scalar, RowsB, Cols, Options, MaxRowsB, MaxCols>& B,
//                        double tolerance = NUMERICAL_ZERO_DIFFERENCE)
// {
//   return A.isApprox(B, tolerance);
// }


// Compares two scalar types for equality up to tolerance
// template<class A, class B>
// bool equals(const A& a, const B& b, const double tolerance = NUMERICAL_ZERO_DIFFERENCE)
// {
//   return boost::math::relative_difference(a, b) <= tolerance;
// }

bool equals(const double a, const double b, const double tolerance = NUMERICAL_ZERO_DIFFERENCE);

template<class DerivedA, class DerivedB>
bool oneGreater(const Eigen::MatrixBase<DerivedA> & A,
                const Eigen::MatrixBase<DerivedB> & B,
                double tolerance = math::NUMERICAL_ZERO_DIFFERENCE)
{
  assertion(A.rows() == B.rows(), "Matrices with different number of rows can't be compared.");
  assertion(A.rows() == B.rows(), "Matrices with different number of cols can't be compared.");
  for (int col = 0; col < A.cols(); ++col)
    for (int row = 0; row < A.rows(); ++row)
      if ( A(row, col) - B(row, col) > tolerance )
        return true;
  
  return false;
}

template<class DerivedA, class DerivedB>
bool oneGreaterEquals(const Eigen::MatrixBase<DerivedA> & A,
                      const Eigen::MatrixBase<DerivedB> & B,
                      double tolerance = math::NUMERICAL_ZERO_DIFFERENCE)
{
  assertion(A.rows() == B.rows(), "Matrices with different number of rows can't be compared.");
  assertion(A.rows() == B.rows(), "Matrices with different number of cols can't be compared.");
  for (int col = 0; col < A.cols(); ++col)
    for (int row = 0; row < A.rows(); ++row) {
      if ( A(row, col) - B(row, col) >= -tolerance )
        return true;
    }
  return false;
}


template<class DerivedA, class DerivedB>
bool allGreater(const Eigen::MatrixBase<DerivedA> & A,
                const Eigen::MatrixBase<DerivedB> & B,
                double tolerance = math::NUMERICAL_ZERO_DIFFERENCE)
{
  assertion(A.rows() == B.rows(), "Matrices with different number of rows can't be compared.");
  assertion(A.rows() == B.rows(), "Matrices with different number of cols can't be compared.");
  for (int col = 0; col < A.cols(); ++col)
    for (int row = 0; row < A.rows(); ++row)
      if ( A(row, col) - B(row, col) <= tolerance )
        return false;
  
  return true;
}

template<class DerivedA, class DerivedB>
bool allGreaterEquals(const Eigen::MatrixBase<DerivedA> & A,
                      const Eigen::MatrixBase<DerivedB> & B,
                      double tolerance = math::NUMERICAL_ZERO_DIFFERENCE)
{
  assertion(A.rows() == B.rows(), "Matrices with different number of rows can't be compared.");
  assertion(A.rows() == B.rows(), "Matrices with different number of cols can't be compared.");
  for (int col = 0; col < A.cols(); ++col)
    for (int row = 0; row < A.rows(); ++row)
      if ( A(row, col) - B(row, col) < -tolerance )
        return false;
  
  return true;
}

template<class Scalar>
typename std::enable_if<std::is_arithmetic<Scalar>::value, bool>
::type greater(Scalar A, Scalar B, Scalar tolerance = NUMERICAL_ZERO_DIFFERENCE)
{
  return A - B > tolerance;
}

template<class Scalar>
typename std::enable_if<std::is_arithmetic<Scalar>::value, bool>
::type greaterEquals(Scalar A, Scalar B, Scalar tolerance = NUMERICAL_ZERO_DIFFERENCE)
{
  return A - B >= -tolerance;
}


template<class Scalar>
typename std::enable_if<std::is_arithmetic<Scalar>::value, bool>
::type smaller(Scalar A, Scalar B, Scalar tolerance = NUMERICAL_ZERO_DIFFERENCE)
{
  return A - B < -tolerance;
}


template<class Scalar>
typename std::enable_if<std::is_arithmetic<Scalar>::value, bool>
::type smallerEquals(Scalar A, Scalar B, Scalar tolerance = NUMERICAL_ZERO_DIFFERENCE)
{
  return A - B <= -tolerance;
}


}
}

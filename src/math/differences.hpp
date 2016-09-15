#pragma once

#include <boost/math/special_functions/relative_difference.hpp>
#include "tarch/la/DynamicVector.h"
#include "tarch/la/Vector.h"

namespace precice {
namespace math {

constexpr double NUMERICAL_ZERO_DIFFERENCE = 1.0e-14;

/// Compares two scalar types for equality up to tolerance
template<class A, class B>
constexpr bool equals(const A a, const B b, const double tolerance = NUMERICAL_ZERO_DIFFERENCE)
{
  return boost::math::relative_difference(a, b) <= tolerance;
}

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


} }

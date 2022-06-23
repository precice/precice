#pragma once

#include "math/constants.hpp"
#include "math/differences.hpp"
#include "math/la.hpp"

namespace precice {
namespace math {

/// Return the sign, one of {-1, 0, 1}
inline int sign(double number)
{
  if (greater(number, 0.0)) {
    return 1;
  } else if (greater(0.0, number)) {
    return -1;
  }
  return 0;
}

/// Computes the power of an integer (x^iexp) using recursion, which is much faster than std::pow(x, iexp)
template <int iexp, typename T>
inline constexpr T pow_int(const T x)
{
  static_assert(iexp >= 0, "Exponent must be an integer greater zero.");

  if (iexp == 0)
    return static_cast<T>(1.);
  else
    // exponentiation by squaring
    return ((iexp % 2 == 1) ? x * pow_int<iexp / 2>(x * x) : pow_int<iexp / 2>(x * x));
}

} // namespace math
} // namespace precice

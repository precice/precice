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

} // namespace math
} // namespace precice

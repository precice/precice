#pragma once

#include "math/differences.hpp"
#include "math/constants.hpp"

namespace precice {
namespace math {

inline int sign (double number)
{
  if ( greater(number, 0.0) ) {
    return 1;
  }
  else if ( greater(0.0, number) ) {
    return -1;
  }
  return 0;
}

}
}

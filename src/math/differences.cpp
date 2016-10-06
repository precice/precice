#include "differences.hpp"

#include <boost/math/special_functions/relative_difference.hpp>

namespace precice {
namespace math {

bool equals(const double a, const double b, const double tolerance)
{
  return boost::math::relative_difference(a, b) <= tolerance;
}


}}

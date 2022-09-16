#include "AbsoluteConvergenceMeasure.hpp"
#include "logging/LogMacros.hpp"
#include "math/differences.hpp"

namespace precice::cplscheme::impl {

AbsoluteConvergenceMeasure::AbsoluteConvergenceMeasure(double convergenceLimit)
    : _convergenceLimit(convergenceLimit)
{
  PRECICE_ASSERT(math::greater(_convergenceLimit, 0.0),
                 "Absolute convergence limit has to be greater than zero!");
}

} // namespace precice::cplscheme::impl

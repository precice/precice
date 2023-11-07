#include "AbsoluteOrRelativeConvergenceMeasure.hpp"
#include "logging/LogMacros.hpp"

namespace precice::cplscheme::impl {

AbsoluteOrRelativeConvergenceMeasure::AbsoluteOrRelativeConvergenceMeasure(double absLimit, double relLimit)
    : _convergenceLimitPercent(relLimit),
      _convergenceLimit(absLimit)

{
  PRECICE_ASSERT(math::greater(_convergenceLimit, 0.0),
                 "Absolute convergence limit has to be greater than zero!");
  PRECICE_ASSERT(math::greater(_convergenceLimitPercent, 0.0) && math::greaterEquals(1.0, _convergenceLimitPercent),
                 "Relative convergence limit has to be in (0;1] !");
}
} // namespace precice::cplscheme::impl

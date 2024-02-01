#include "ResidualRelativeConvergenceMeasure.hpp"
#include "logging/LogMacros.hpp"
#include "math/differences.hpp"

namespace precice::cplscheme::impl {

ResidualRelativeConvergenceMeasure::ResidualRelativeConvergenceMeasure(double convergenceLimitPercent)
    : _convergenceLimitPercent(convergenceLimitPercent)
{
  PRECICE_ASSERT(math::greater(_convergenceLimitPercent, 0.0) && math::greaterEquals(1.0, _convergenceLimitPercent),
                 "Relative convergence limit has to be in ]0;1] !");
}
} // namespace precice::cplscheme::impl

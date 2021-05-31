#include "RelativeConvergenceMeasure.hpp"
#include "logging/LogMacros.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

RelativeConvergenceMeasure::RelativeConvergenceMeasure(double convergenceLimitPercent)
    : _convergenceLimitPercent(convergenceLimitPercent)
{
  PRECICE_ASSERT(math::greater(_convergenceLimitPercent, 0.0) && math::greaterEquals(1.0, _convergenceLimitPercent),
                 "Relative convergence limit has to be in ]0;1] !");
}
} // namespace impl
} // namespace cplscheme
} // namespace precice

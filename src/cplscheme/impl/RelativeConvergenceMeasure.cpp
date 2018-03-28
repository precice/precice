#include "RelativeConvergenceMeasure.hpp"
#include "utils/Globals.hpp"

namespace precice
{
namespace cplscheme
{
namespace impl
{

logging::Logger RelativeConvergenceMeasure::
    _log("cplscheme::RelativeConvergenceMeasure");

RelativeConvergenceMeasure::RelativeConvergenceMeasure(
    double convergenceLimitPercent)
    : ConvergenceMeasure(),
      _convergenceLimitPercent(convergenceLimitPercent),
      _normDiff(0.0),
      _norm(0.0),
      _isConvergence(false)
{
  CHECK(math::greater(_convergenceLimitPercent, 0.0) && math::greaterEquals(1.0, _convergenceLimitPercent),
        "Relative convergence limit has to be in ]0;1] !");
}
}
}
} // namespace precice, cplscheme, impl

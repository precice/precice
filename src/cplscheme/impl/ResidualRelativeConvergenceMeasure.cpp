#include "ResidualRelativeConvergenceMeasure.hpp"
#include "math/math.hpp"
#include "utils/Globals.hpp"

namespace precice
{
namespace cplscheme
{
namespace impl
{

logging::Logger ResidualRelativeConvergenceMeasure::
    _log("cplscheme::ResidualRelativeConvergenceMeasure");

ResidualRelativeConvergenceMeasure::ResidualRelativeConvergenceMeasure(
    double convergenceLimitPercent)
    : ConvergenceMeasure(),
      _convergenceLimitPercent(convergenceLimitPercent),
      _isFirstIteration(true),
      _normFirstResidual(std::numeric_limits<double>::max()),
      _normDiff(0.0),
      _isConvergence(false)
{
  CHECK(math::greater(_convergenceLimitPercent, 0.0) && math::greaterEquals(1.0, _convergenceLimitPercent),
        "Relative convergence limit has to be in ]0;1] !");
}
}
}
} // namespace precice, cplscheme, impl

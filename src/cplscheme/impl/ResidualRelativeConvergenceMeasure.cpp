#include "ResidualRelativeConvergenceMeasure.hpp"
#include "math/math.hpp"

namespace precice
{
namespace cplscheme
{
namespace impl
{


ResidualRelativeConvergenceMeasure::ResidualRelativeConvergenceMeasure(double convergenceLimitPercent)
  : _convergenceLimitPercent(convergenceLimitPercent)
{
  P_CHECK(math::greater(_convergenceLimitPercent, 0.0) && math::greaterEquals(1.0, _convergenceLimitPercent),
        "Relative convergence limit has to be in ]0;1] !");
}
}
}
} // namespace precice, cplscheme, impl

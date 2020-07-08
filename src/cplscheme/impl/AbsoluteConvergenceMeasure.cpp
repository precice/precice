#include "AbsoluteConvergenceMeasure.hpp"
#include "logging/LogMacros.hpp"
#include "math/differences.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

AbsoluteConvergenceMeasure::AbsoluteConvergenceMeasure(double convergenceLimit)
    : _convergenceLimit(convergenceLimit)
{
  PRECICE_CHECK(not math::greaterEquals(0.0, _convergenceLimit),
                "Absolute convergence limit has to be greater than zero!");
}

} // namespace impl
} // namespace cplscheme
} // namespace precice

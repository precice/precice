#include "AbsoluteConvergenceMeasure.hpp"
#include "math/math.hpp"

namespace precice
{
namespace cplscheme
{
namespace impl
{

AbsoluteConvergenceMeasure::AbsoluteConvergenceMeasure(double convergenceLimit)
  : _convergenceLimit(convergenceLimit)
{
  P_CHECK(not math::greaterEquals(0.0, _convergenceLimit),
        "Absolute convergence limit has to be greater than zero!");
}

}}} // namespace precice, cplscheme, impl

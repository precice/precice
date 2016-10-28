#include "AbsoluteConvergenceMeasure.hpp"
#include "math/math.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

logging::Logger AbsoluteConvergenceMeasure::
   _log ( "precice::cplscheme::AbsoluteConvergenceMeasure" );

AbsoluteConvergenceMeasure:: AbsoluteConvergenceMeasure
(
   double convergenceLimit )
:
   ConvergenceMeasure (),
   _convergenceLimit ( convergenceLimit ),
   _normDiff(0.0),
   _isConvergence ( false )
{
   CHECK ( ! math::greaterEquals(0.0, _convergenceLimit),
           "Absolute convergence limit has to be greater than zero!" );
}

}}} // namespace precice, cplscheme

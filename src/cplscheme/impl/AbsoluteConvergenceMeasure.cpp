#include "AbsoluteConvergenceMeasure.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

tarch::logging::Log AbsoluteConvergenceMeasure::
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
   preciceCheck ( ! tarch::la::greaterEquals(0.0, _convergenceLimit),
                  "AbsoluteConvergenceLimit()", "Absolute convergence limit "
                  << "has to be greater than zero!" );
}

}}} // namespace precice, cplscheme

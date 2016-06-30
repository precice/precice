#include "ResidualRelativeConvergenceMeasure.hpp"
#include "utils/Globals.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

logging::Logger ResidualRelativeConvergenceMeasure::
   _log ( "precice::cplscheme::ResidualRelativeConvergenceMeasure" );


ResidualRelativeConvergenceMeasure:: ResidualRelativeConvergenceMeasure
(
   double convergenceLimitPercent )
:
   ConvergenceMeasure (),
   _convergenceLimitPercent ( convergenceLimitPercent ),
   _isFirstIteration ( true ),
   _normFirstResidual ( std::numeric_limits<double>::max() ),
   _normDiff(0.0),
   _isConvergence ( false )
{
   preciceCheck ( tarch::la::greater(_convergenceLimitPercent, 0.0)
                  && tarch::la::greaterEquals(1.0, _convergenceLimitPercent),
                  "ResidualRelativeConvergenceMeasure()", "Relative convergence limit "
                  << "has in ]0;1] !" );
}

}}} // namespace precice, cplscheme, impl

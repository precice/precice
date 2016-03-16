// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "RelativeConvergenceMeasure.hpp"
#include "utils/Globals.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

logging::Logger RelativeConvergenceMeasure::
   _log ( "precice::cplscheme::RelativeConvergenceMeasure" );


RelativeConvergenceMeasure:: RelativeConvergenceMeasure
(
   double convergenceLimitPercent )
:
   ConvergenceMeasure (),
   _convergenceLimitPercent ( convergenceLimitPercent ),
   _normDiff(0.0),
   _norm(0.0),
   _isConvergence ( false )
{
   preciceCheck ( tarch::la::greater(_convergenceLimitPercent, 0.0)
                  && tarch::la::greaterEquals(1.0, _convergenceLimitPercent),
                  "RelativeConvergenceLimit()", "Relative convergence limit "
                  << "has in ]0;1] !" );
}

//void RelativeConvergenceMeasure:: startMeasurement ()
//{
//   _sumSquaredDifferences = 0.0;
//   _sumSquaredNew = 0.0;
//   _isConvergence = false;
//}

//void RelativeConvergenceMeasure:: finishMeasurement ()
//{
//   double twoNormDifferences = std::sqrt ( _sumSquaredDifferences );
//   double twoNormNew = std::sqrt ( _sumSquaredNew );
//   _isConvergence = tarch::la::greaterEquals (
//         twoNormNew * _convergenceLimitPercent, twoNormDifferences );
//   preciceInfo ( "finsihMeasurement()", "Relative convergence measure: "
//                 << "two-norm differences = " << twoNormDifferences
//                 << ", convergence limit = "
//                 << twoNormNew * _convergenceLimitPercent
//                 << ", convergence = " << _isConvergence );
//}

}}} // namespace precice, cplscheme, impl

// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
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

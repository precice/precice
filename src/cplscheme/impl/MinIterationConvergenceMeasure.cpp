// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "MinIterationConvergenceMeasure.hpp"
#include "utils/Globals.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

logging::Logger MinIterationConvergenceMeasure::
  _log("precice::cplscheme::MinIterationConvergenceMeasure");

MinIterationConvergenceMeasure:: MinIterationConvergenceMeasure
(
  int minimumIterationCount )
:
  ConvergenceMeasure(),
  _minimumIterationCount(minimumIterationCount),
  _currentIteration(0),
  _isConvergence(false)
{}

void MinIterationConvergenceMeasure:: newMeasurementSeries()
{
  _currentIteration = 0;
  _isConvergence = false;
}

}}} // namespace precice, cplscheme, impl

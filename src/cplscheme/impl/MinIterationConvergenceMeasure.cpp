#include "MinIterationConvergenceMeasure.hpp"
#include "utils/Globals.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

logging::Logger MinIterationConvergenceMeasure::
  _log("cplscheme::MinIterationConvergenceMeasure");

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

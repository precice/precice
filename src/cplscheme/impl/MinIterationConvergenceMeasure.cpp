#include "MinIterationConvergenceMeasure.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

MinIterationConvergenceMeasure::MinIterationConvergenceMeasure(
    int minimumIterationCount)
    : ConvergenceMeasure(),
      _minimumIterationCount(minimumIterationCount)
{
}

void MinIterationConvergenceMeasure::newMeasurementSeries()
{
  _currentIteration = 0;
  _isConvergence    = false;
}
} // namespace impl
} // namespace cplscheme
} // namespace precice

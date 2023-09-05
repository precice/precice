#include <algorithm>

#include "cplscheme/CouplingScheme.hpp"
#include "logging/LogMacros.hpp"
#include "math/bspline.hpp"
#include "math/differences.hpp"
#include "mesh/Data.hpp"
#include "time/Time.hpp"
#include "time/Waveform.hpp"

namespace precice::time {

Waveform::Waveform(const int degree)
    : _degree(degree)
{
  PRECICE_ASSERT(Time::MIN_WAVEFORM_DEGREE <= _degree && _degree <= Time::MAX_WAVEFORM_DEGREE);
}

int Waveform::getDegree() const
{
  return _degree;
}

time::Storage &Waveform::timeStepsStorage()
{
  return _timeStepsStorage;
}

const time::Storage &Waveform::timeStepsStorage() const
{
  return _timeStepsStorage;
}

Eigen::VectorXd Waveform::sample(double normalizedDt) const
{
  const int usedDegree = computeUsedDegree(_degree, _timeStepsStorage.nTimes());

  PRECICE_ASSERT(math::equals(this->_timeStepsStorage.maxStoredNormalizedDt(), time::Storage::WINDOW_END), this->_timeStepsStorage.maxStoredNormalizedDt()); // sampling is only allowed, if a window is complete.

  if (_degree == 0) {
    return this->_timeStepsStorage.getSampleAtOrAfter(normalizedDt).values;
  }

  PRECICE_ASSERT(usedDegree >= 1);

  const auto data = _timeStepsStorage.getTimesAndValues();

  return math::bspline::interpolateAt(data.first, data.second, usedDegree, normalizedDt);
}

int Waveform::computeUsedDegree(int requestedDegree, int numberOfAvailableSamples) const
{
  int usedDegree = -1;
  PRECICE_ASSERT(requestedDegree <= 3);
  if (requestedDegree == 0 || numberOfAvailableSamples < 2) {
    usedDegree = 0;
  } else if (requestedDegree == 1 || numberOfAvailableSamples < 3) {
    usedDegree = 1;
  } else if (requestedDegree == 2 || numberOfAvailableSamples < 4) {
    usedDegree = 2;
  } else if (requestedDegree == 3 || numberOfAvailableSamples < 5) {
    usedDegree = 3;
  } else {
    PRECICE_ASSERT(false); // not supported
  }
  return usedDegree;
}

} // namespace precice::time

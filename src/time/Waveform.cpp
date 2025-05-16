#include <algorithm>

#include "cplscheme/CouplingScheme.hpp"
#include "logging/LogMacros.hpp"
#include "math/differences.hpp"
#include "mesh/Data.hpp"
#include "time/Waveform.hpp"

namespace precice::time {

Waveform::Waveform(const int degree)
{
  _timeStepsStorage.setInterpolationDegree(degree);
}

time::Storage &Waveform::timeStepsStorage()
{
  return _timeStepsStorage;
}

const time::Storage &Waveform::timeStepsStorage() const
{
  return _timeStepsStorage;
}

SampleResult Waveform::sample(double time) const
{
  return _timeStepsStorage.sample(time);
}
} // namespace precice::time

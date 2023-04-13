#include "time/Waveform.hpp"
#include <algorithm>
#include <unsupported/Eigen/Splines>
#include "cplscheme/CouplingScheme.hpp"
#include "logging/LogMacros.hpp"
#include "math/differences.hpp"
#include "time/Time.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice::time {

Waveform::Waveform(
    const int interpolationOrder)
    : _interpolationOrder(interpolationOrder), _storage(0, interpolationOrder)
{
  PRECICE_ASSERT(Time::MIN_INTERPOLATION_ORDER <= _interpolationOrder && _interpolationOrder <= Time::MAX_INTERPOLATION_ORDER);
}

int Waveform::getInterpolationOrder() const
{
  return _interpolationOrder;
}

void Waveform::initialize(const Eigen::VectorXd &values)
{
  PRECICE_ASSERT(_storage.nTimes() == 0);
  _storage.initialize(values);
  PRECICE_ASSERT(_interpolationOrder >= Time::MIN_INTERPOLATION_ORDER);
}

void Waveform::store(const Eigen::VectorXd &values, double normalizedDt)
{
  if (math::equals(_storage.maxStoredNormalizedDt(), time::Storage::WINDOW_END)) { // reached end of window and trying to write new data from next window. Clearing window first.
    _storage.clear(false);
  }
  PRECICE_ASSERT((_storage.nTimes() == 0) || (values.size() == _storage.nDofs()), "Size of new data does not match size of existing data", values.size(), _storage.nDofs());
  _storage.setValuesAtTime(normalizedDt, values);
}

Eigen::VectorXd Waveform::sample(double normalizedDt)
{
  PRECICE_ASSERT(math::equals(this->_storage.maxStoredNormalizedDt(), time::Storage::WINDOW_END), this->_storage.maxStoredNormalizedDt()); // sampling is only allowed, if a window is complete.

  return this->_storage.sampleAt(normalizedDt);
}

void Waveform::moveToNextWindow()
{
  _storage.move();
}

} // namespace precice::time

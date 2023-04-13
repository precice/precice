#include "time/Storage.hpp"
#include <unsupported/Eigen/Splines>
#include "cplscheme/CouplingScheme.hpp"
#include "math/differences.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/assertion.hpp"

namespace precice::time {

const double Storage::WINDOW_START = 0.0;

const double Storage::WINDOW_END = 1.0;

Storage::Storage(int extrapolationOrder, int interpolationOrder)
    : _sampleStorage{}, _extrapolationOrder(extrapolationOrder), _interpolationOrder(interpolationOrder)
{
}

void Storage::initialize(Eigen::VectorXd values)
{
  _sampleStorage.emplace_back(std::make_pair(WINDOW_START, values));
  _sampleStorage.emplace_back(std::make_pair(WINDOW_END, values));
}

Eigen::VectorXd Storage::getValuesAtTime(double time)
{
  for (auto &sample : _sampleStorage) {
    if (math::equals(sample.first, time)) {
      return sample.second;
    }
  }
  PRECICE_ASSERT(false, "no values found!", time, getTimes());
}

void Storage::setValuesAtTime(double time, Eigen::VectorXd values, bool mustOverwriteExisting)
{
  PRECICE_ASSERT(math::smallerEquals(WINDOW_START, time), "Setting values outside of valid range!");
  PRECICE_ASSERT(math::smallerEquals(time, WINDOW_END), "Sampling outside of valid range!");
  if (!mustOverwriteExisting) {
    PRECICE_ASSERT(math::smaller(maxStoredNormalizedDt(), time), maxStoredNormalizedDt(), time, "Trying to overwrite existing values or to write values with a time that is too small. Please use clear(), if you want to reset the storage.");
    _sampleStorage.emplace_back(std::make_pair(time, values));
  } else {
    // check that key "time" exists.
    auto sample = std::find_if(_sampleStorage.begin(), _sampleStorage.end(), [&time](const auto &s) { return math::equals(s.first, time); });
    PRECICE_ASSERT(sample != _sampleStorage.end(), time, "Key does not exist, cannot overwrite values.");
    // overwrite values at "time"
    for (auto &sample : _sampleStorage) {
      if (math::equals(sample.first, time)) {
        sample.second = values;
        return;
      }
    }
    PRECICE_ASSERT(false, "unreachable!");
  }
}

double Storage::maxStoredNormalizedDt()
{
  if (_sampleStorage.size() == 0) {
    return -1; // invalid return
  } else {
    return _sampleStorage.back().first;
  }
}

int Storage::nTimes()
{
  return _sampleStorage.size();
}

int Storage::nDofs()
{
  PRECICE_ASSERT(_sampleStorage.size() > 0);
  return _sampleStorage[0].second.size();
}

void Storage::move()
{
  PRECICE_ASSERT(nTimes() > 0);
  auto valuesAtBeginning = getValuesAtEnd();
  auto valuesAtEnd       = computeExtrapolation();
  _sampleStorage.clear();
  _sampleStorage.emplace_back(std::make_pair(WINDOW_START, valuesAtBeginning));
  _sampleStorage.emplace_back(std::make_pair(WINDOW_END, valuesAtEnd));
}

void Storage::clear(bool keepWindowStart)
{
  if (keepWindowStart) {
    PRECICE_ASSERT(nTimes() > 0, "Storage does not contain any data!");
    Eigen::VectorXd keep = getValuesAtBeginning(); // we keep data at _storageDict[0.0]
    _sampleStorage.clear();
    _sampleStorage.emplace_back(std::make_pair(WINDOW_START, keep));
  } else {
    _sampleStorage.clear();
  }
}

Eigen::VectorXd Storage::getValuesAtOrAfter(double before)
{
  auto sample = std::find_if(_sampleStorage.begin(), _sampleStorage.end(), [&before](const auto &s) { return math::greaterEquals(s.first, before); });
  PRECICE_ASSERT(sample != _sampleStorage.end(), "no values found!");

  return sample->second;
}

Eigen::VectorXd Storage::getTimes()
{
  auto times = Eigen::VectorXd(nTimes());
  for (int i = 0; i < times.size(); i++) {
    times[i] = _sampleStorage[i].first;
  }
  return times;
}

std::pair<Eigen::VectorXd, Eigen::MatrixXd> Storage::getTimesAndValues()
{
  auto times  = Eigen::VectorXd(nTimes());
  auto values = Eigen::MatrixXd(nDofs(), nTimes());
  for (int i = 0; i < times.size(); i++) {
    times[i]      = _sampleStorage[i].first;
    values.col(i) = _sampleStorage[i].second;
  }
  return std::make_pair(times, values);
}

Eigen::VectorXd Storage::computeExtrapolation()
{
  if (_extrapolationOrder == 0 || _extrapolationOrder == cplscheme::CouplingScheme::UNDEFINED_EXTRAPOLATION_ORDER) {
    return getValuesAtEnd(); // use values at end of window as initial guess for next
  } else if (_extrapolationOrder == 1) {
    return 2 * getValuesAtEnd() - getValuesAtBeginning(); // use linear extrapolation from window at beginning and end of window.
  }
  PRECICE_UNREACHABLE("Invalid _extrapolationOrder")
}

Eigen::VectorXd Storage::getValuesAtBeginning()
{
  PRECICE_ASSERT(math::equals(_sampleStorage.front().first, WINDOW_START), _sampleStorage.front().first);
  return _sampleStorage.front().second;
}

Eigen::VectorXd Storage::getValuesAtEnd()
{
  PRECICE_ASSERT(math::equals(_sampleStorage.back().first, WINDOW_END), _sampleStorage.back().first);
  return _sampleStorage.back().second;
}

// helper function to compute x(t) from given data (x0,t0), (x1,t1), ..., (xn,tn) via B-spline interpolation (implemented using Eigen).
Eigen::VectorXd Storage::bSplineInterpolationAt(double t, Eigen::VectorXd ts, Eigen::MatrixXd xs, int splineDegree)
{
  // organize data in columns. Each column represents one sample in time.
  PRECICE_ASSERT(xs.cols() == ts.size());
  const int ndofs = xs.rows(); // number of dofs. Each dof needs it's own interpolant.

  Eigen::VectorXd interpolated(ndofs);

  const int splineDimension = 1;

  for (int i = 0; i < ndofs; i++) {
    auto spline     = Eigen::SplineFitting<Eigen::Spline<double, splineDimension>>::Interpolate(xs.row(i), splineDegree, ts);
    interpolated[i] = spline(t)[0]; // get component of spline associated with xs.row(i)
  }

  return interpolated;
}

Eigen::VectorXd Storage::sampleAt(double normalizedDt)
{
  const int usedOrder = computeUsedOrder(_interpolationOrder, nTimes());

  PRECICE_ASSERT(math::equals(this->maxStoredNormalizedDt(), time::Storage::WINDOW_END), this->maxStoredNormalizedDt()); // sampling is only allowed, if a window is complete.

  if (_interpolationOrder == 0) {
    return this->getValuesAtOrAfter(normalizedDt);
  }

  PRECICE_ASSERT(usedOrder >= 1);

  auto data = getTimesAndValues();
  return bSplineInterpolationAt(normalizedDt, data.first, data.second, usedOrder);
}

int Storage::computeUsedOrder(int requestedOrder, int numberOfAvailableSamples)
{
  int usedOrder = -1;
  PRECICE_ASSERT(requestedOrder <= 3);
  if (requestedOrder == 0 || numberOfAvailableSamples < 2) {
    usedOrder = 0;
  } else if (requestedOrder == 1 || numberOfAvailableSamples < 3) {
    usedOrder = 1;
  } else if (requestedOrder == 2 || numberOfAvailableSamples < 4) {
    usedOrder = 2;
  } else if (requestedOrder == 3 || numberOfAvailableSamples < 5) {
    usedOrder = 3;
  } else {
    PRECICE_ASSERT(false); // not supported
  }
  return usedOrder;
}

} // namespace precice::time

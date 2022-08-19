#include "time/Waveform.hpp"
#include <algorithm>
#include <eigen3/unsupported/Eigen/Splines>
#include "cplscheme/CouplingScheme.hpp"
#include "logging/LogMacros.hpp"
#include "time/Time.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice {
namespace time {

Waveform::Waveform(
    const int interpolationOrder)
    : _interpolationOrder(interpolationOrder)
{
  PRECICE_ASSERT(not _storageIsInitialized);
  PRECICE_ASSERT(Time::MIN_INTERPOLATION_ORDER <= _interpolationOrder && _interpolationOrder <= Time::MAX_INTERPOLATION_ORDER);
}

int Waveform::getInterpolationOrder() const
{
  return _interpolationOrder;
}

void Waveform::initialize(const Eigen::VectorXd &values)
{
  int storageSize;
  _timeStepsStorage = std::map<double, Eigen::VectorXd>{};
  PRECICE_ASSERT(_interpolationOrder >= Time::MIN_INTERPOLATION_ORDER);
  _timeStepsStorage[0.0] = Eigen::VectorXd(values);
  _timeStepsStorage[1.0] = Eigen::VectorXd(values);
  _numberOfStoredSamples = 2; // the first sample is automatically initialized as zero and stored.
  _storageIsInitialized  = true;
  PRECICE_ASSERT(this->valuesSize() == values.size());
}

void Waveform::store(const Eigen::VectorXd &values)
{
  PRECICE_ASSERT(_storageIsInitialized);
  double normalizedDtForEnd = 1.0; // use dt associated with end of window
  this->store(values, normalizedDtForEnd);
}

void Waveform::store(const Eigen::VectorXd &values, double normalizedDt)
{
  PRECICE_ASSERT(_storageIsInitialized);
  // dt has to be in interval (0.0, 1.0]
  PRECICE_ASSERT(normalizedDt > 0.0); // cannot override value at beginning of window. It is locked!
  PRECICE_ASSERT(normalizedDt <= 1.0);
  this->storeAt(values, normalizedDt);
}

double Waveform::maxStoredDt()
{
  double maxDt = -1;
  for (auto timeStep : _timeStepsStorage) {
    if (timeStep.first > maxDt) {
      maxDt = timeStep.first;
    }
  }
  PRECICE_ASSERT(maxDt >= 0);
  return maxDt;
}

void Waveform::storeAt(const Eigen::VectorXd values, double dt)
{
  if (maxStoredDt() < 1.0) { // did not reach end of window yet, so dt has to strictly increase
    PRECICE_ASSERT(dt > maxStoredDt());
  } else {                                         // reached end of window and trying to write new data from next window. Clearing window first.
    Eigen::VectorXd keep = _timeStepsStorage[0.0]; // we keep data at _timeStepsStorage[0.0]
    _timeStepsStorage.clear();
    _timeStepsStorage[0.0] = keep;
    _numberOfStoredSamples = 1;
  }
  PRECICE_ASSERT(values.size() == this->valuesSize(), values.size(), this->valuesSize());
  this->_timeStepsStorage[dt] = Eigen::VectorXd(values);
  _numberOfStoredSamples++;
}

// helper function to compute x(t) from given data (x0,t0), (x1,t1), ..., (xn,tn) via third degree B-spline interpolation (implemented using Eigen). Alternatives are cubic spline interpolation (see https://en.wikipedia.org/wiki/Spline_interpolation#Algorithm_to_find_the_interpolating_cubic_spline), but this is not offered by Eigen or Boost.
Eigen::VectorXd bSplineInterpolationAt(double t, Eigen::VectorXd ts, Eigen::MatrixXd xs, int splineDegree)
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

Eigen::VectorXd Waveform::sample(double normalizedDt)
{
  PRECICE_ASSERT(_storageIsInitialized);
  PRECICE_ASSERT(normalizedDt >= 0, "Sampling outside of valid range!");
  PRECICE_ASSERT(normalizedDt <= 1, "Sampling outside of valid range!");

  const int usedOrder = computeUsedOrder(_interpolationOrder, _numberOfStoredSamples);

  PRECICE_ASSERT(maxStoredDt() == 1.0); // sampling is only allowed, if a window is complete.

  // @TODO: Improve efficiency: Check whether key = normalizedDt is in _timeStepsStorage. If yes, just get value and return. No need for interpolation.

  if (_interpolationOrder == 0) {
    // @TODO: Remove constant interpolation in preCICE v3.0? Usecase is unclear and does not generalize well with BSpline interpolation. It's also not 100% clear what to do at the jump.
    // constant interpolation = just use sample at the end of the window: x(dt) = x^t
    // At beginning of window use result from last window x(0) = x^(t-1)
    return Eigen::VectorXd(this->_timeStepsStorage[findTimeAfter(normalizedDt)]);
  }

  PRECICE_ASSERT(usedOrder >= 1);

  /** @TODO for higher-order quadratic interpolation there are several possibilities:
   * 1. Use data from this window and last window. Then we do not need to consider any samples from subcycling
   * 2. Use data from this window. Requires at least polynomial degree p substeps in window to create BSpline.
   * 3. Use data from this window, but perform a p-th least squares fit. If we don't do subcycling the system is underdetermined, if we do 2 substeps this option is identical to option 2. If we do 3 or more substeps we will get a least squares fit. Important: Might lead to discontinuities at the window boundary!
   **/

  auto timesAscending = getTimesAscending();
  auto nTimes         = timesAscending.size();
  auto nDofs          = this->_timeStepsStorage[0.0].size();
  PRECICE_ASSERT(timesAscending[0] == 0.0);
  PRECICE_ASSERT(timesAscending[nTimes - 1] == 1.0);
  Eigen::MatrixXd dataAscending(nDofs, nTimes);
  int             i = 0;
  for (int i = 0; i < nTimes; i++) {
    dataAscending.col(i) = this->_timeStepsStorage[timesAscending[i]];
  }
  return bSplineInterpolationAt(normalizedDt, timesAscending, dataAscending, usedOrder);
}

void Waveform::moveToNextWindow()
{
  PRECICE_ASSERT(_storageIsInitialized);
  auto initialGuess = this->sample(1.0); // use value at end of window as initial guess for next
  _timeStepsStorage.clear();
  _timeStepsStorage[0.0] = Eigen::VectorXd(initialGuess);
  _timeStepsStorage[1.0] = Eigen::VectorXd(initialGuess); // initial guess is always constant extrapolation
  _numberOfStoredSamples = 2;
}

int Waveform::valuesSize()
{
  PRECICE_ASSERT(_storageIsInitialized);
  return _timeStepsStorage[0.0].size();
}

int Waveform::computeUsedOrder(int requestedOrder, int numberOfAvailableSamples)
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

Eigen::VectorXd Waveform::getTimesAscending()
{
  // create std::vector with all keys
  std::vector<double> keys;
  for (auto timeStep : _timeStepsStorage) {
    keys.push_back(timeStep.first);
  }

  // sort vector
  std::sort(keys.begin(), keys.end());

  // copy data into Eigen::VectorXd to return
  auto times = Eigen::VectorXd(keys.size());
  for (int i = 0; i < keys.size(); i++) {
    times[i] = keys[i];
  }
  return times;
}

// @todo improve efficiency by using some ordered list to store times.
double Waveform::findTimeAfter(double normalizedDt)
{
  PRECICE_ASSERT(0.0 <= normalizedDt);
  PRECICE_ASSERT(normalizedDt <= 1.0);

  double timeAfter = 1.0;

  for (auto timeStep : _timeStepsStorage) {
    if (normalizedDt <= timeStep.first && timeStep.first <= timeAfter) { // current timeStep is after normalizedDt and earlier than current timeAfter
      timeAfter = timeStep.first;
    }
  }

  return timeAfter;
}

} // namespace time
} // namespace precice

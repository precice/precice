#include "time/Waveform.hpp"
#include <algorithm>
#include <eigen3/unsupported/Eigen/Splines>
#include "cplscheme/CouplingScheme.hpp"
#include "logging/LogMacros.hpp"
#include "math/differences.hpp"
#include "time/Time.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice::time {

Waveform::Waveform(
    const int interpolationOrder)
    : _interpolationOrder(interpolationOrder)
{
  PRECICE_ASSERT(_timeStepsStorage.size() == 0);
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
}

void Waveform::store(const Eigen::VectorXd &values, double normalizedDt)
{
  PRECICE_ASSERT(_timeStepsStorage.size() > 0);
  // dt has to be in interval (0.0, 1.0]
  PRECICE_ASSERT(normalizedDt > 0.0); // cannot override value at beginning of window. It is locked!
  PRECICE_ASSERT(normalizedDt <= 1.0);

  if (math::equals(maxStoredNormalizedDt(), 1.0)) { // reached end of window and trying to write new data from next window. Clearing window first.
    Eigen::VectorXd keep = _timeStepsStorage[0.0];  // we keep data at _timeStepsStorage[0.0]
    _timeStepsStorage.clear();
    _timeStepsStorage[0.0] = keep;
  } else { // did not reach end of window yet, so dt has to strictly increase
    PRECICE_ASSERT(normalizedDt > maxStoredNormalizedDt(), normalizedDt, maxStoredNormalizedDt());
  }
  PRECICE_ASSERT(values.size() == _timeStepsStorage[0.0].size());
  this->_timeStepsStorage[normalizedDt] = Eigen::VectorXd(values);
}

double Waveform::maxStoredNormalizedDt()
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

// helper function to compute x(t) from given data (x0,t0), (x1,t1), ..., (xn,tn) via B-spline interpolation (implemented using Eigen).
Eigen::VectorXd bSplineInterpolationAt(double t, Eigen::VectorXd ts, Eigen::MatrixXd xs, int splineDegree)
{
  // organize data in columns. Each column represents one sample in time.
  PRECICE_ASSERT(xs.cols() == ts.size());
  const int ndofs = xs.rows(); // number of dofs. Each dof needs its own interpolant.

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
  PRECICE_ASSERT(_timeStepsStorage.size() > 0);
  PRECICE_ASSERT(normalizedDt >= 0, "Sampling outside of valid range!");
  PRECICE_ASSERT(normalizedDt <= 1, "Sampling outside of valid range!");

  const int usedOrder = computeUsedOrder(_interpolationOrder, _timeStepsStorage.size());

  PRECICE_ASSERT(math::equals(maxStoredNormalizedDt(), 1.0), maxStoredNormalizedDt()); // sampling is only allowed, if a window is complete.

  // @TODO: Improve efficiency: Check whether key = normalizedDt is in _timeStepsStorage. If yes, just get value and return. No need for interpolation.

  if (_interpolationOrder == 0) {
    // @TODO: Remove constant interpolation in preCICE v3.0? Usecase is unclear and does not generalize well with BSpline interpolation. It's also not 100% clear what to do at the jump.
    // constant interpolation = just use sample at the end of the window: x(dt) = x^t
    // At beginning of window use result from last window x(0) = x^(t-1)
    return Eigen::VectorXd(this->_timeStepsStorage[findTimeAfter(normalizedDt)]);
  }

  PRECICE_ASSERT(usedOrder >= 1);

  auto timesAscending = getTimesAscending();
  auto nTimes         = timesAscending.size();
  auto nDofs          = this->_timeStepsStorage[0.0].size();
  PRECICE_ASSERT(math::equals(timesAscending[0], 0.0));
  PRECICE_ASSERT(math::equals(timesAscending[nTimes - 1], 1.0));
  Eigen::MatrixXd dataAscending(nDofs, nTimes);
  int             i = 0;
  for (int i = 0; i < nTimes; i++) {
    dataAscending.col(i) = this->_timeStepsStorage[timesAscending[i]];
  }
  return bSplineInterpolationAt(normalizedDt, timesAscending, dataAscending, usedOrder);
}

void Waveform::moveToNextWindow()
{
  PRECICE_ASSERT(_timeStepsStorage.size() > 0);
  auto initialGuess = this->sample(maxStoredNormalizedDt()); // use value at end of window as initial guess for next
  _timeStepsStorage.clear();
  _timeStepsStorage[0.0] = Eigen::VectorXd(initialGuess);
  _timeStepsStorage[1.0] = Eigen::VectorXd(initialGuess); // initial guess is always constant extrapolation
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

} // namespace precice::time

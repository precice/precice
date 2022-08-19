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

  // need to provide 2d coordinates (t,x) to make all points unique (0,0,1) is forbidden ((0,0), (0.5,0), (1,1)) is allowed
  const int splineDimension = 2;
  for (int i = 0; i < ndofs; i++) {
    Eigen::MatrixXd data(2, xs.cols());
    for (int j = 0; j < xs.cols(); j++) {
      data.row(0)[j] = ts[j];
      data.row(1)[j] = xs.row(i)[j];
    }
    // create know vector
    Eigen::VectorXd knots(ts.size());
    double          tmin = ts[0];
    double          tmax = ts[ts.size() - 1];
    for (int j = 0; j < xs.cols(); j++) {
      knots[j] = (ts[j] - tmin) / (tmax - tmin); // scale knots to interval [0,1]. Use ts as reference point.
    }
    double tScaled  = (t - tmin) / (tmax - tmin); // scale sampling time to interval [0,1];
    auto   spline   = Eigen::SplineFitting<Eigen::Spline<double, splineDimension>>::Interpolate(data, splineDegree, knots);
    interpolated[i] = spline(tScaled)[1]; // get component of spline associated with xs.row(i)
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
    // constant interpolation = just use sample at the end of the window: x(dt) = x^t
    // At beginning of window use result from last window x(0) = x^(t-1)
    return Eigen::VectorXd(this->_timeStepsStorage[findTimeAfter(normalizedDt)]);
  }

  PRECICE_ASSERT(usedOrder >= 1);

  /** @TODO for quadratic interpolation there are several possibilities:
   * 1. Use data from this window and last window. Then we do not need to consider any samples from subcycling
   * 2. Use data from this window. Requires at least 2 substeps in window. Note: For 3 or more substeps we will directly apply cubic spline interpolation, because this is the natural choice and piecewise quadratic spline interpolation is very academic.
   * 3. Use data from this window, but perform a least squares fit. If we don't do subcycling the system is underdetermined, if we do 2 substeps this option is identical to option 2. If we do 3 or more substeps we will get a least squares fit. Important: Might lead to discontinuities at the window boundary!
   **/

  double          low      = 0.0;
  double          high     = 1.0;
  auto            timesVec = getTimesAscending(low, high);
  Eigen::VectorXd timesAscending(timesVec.size());
  Eigen::MatrixXd dataAscending(this->_timeStepsStorage[0.0].size(), timesVec.size());
  int             i = 0;
  for (int i = 0; i < timesVec.size(); i++) {
    timesAscending[i]    = timesVec[i];
    dataAscending.col(i) = this->_timeStepsStorage[timesVec[i]];
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

std::vector<double> Waveform::getTimesAscending(double low, double high)
{
  // create std::vector with all keys
  std::vector<double> times;
  for (auto timeStep : _timeStepsStorage) {
    times.push_back(timeStep.first);
  }
  std::vector<double> timesFiltered;
  std::copy_if(times.begin(), times.end(), std::back_inserter(timesFiltered), [low](double e) { return e >= low; }); // remove elements below lower bound
  times = timesFiltered;
  timesFiltered.clear();
  std::copy_if(times.begin(), times.end(), std::back_inserter(timesFiltered), [high](double e) { return e <= high; }); // remove elements above higher bound
  times = timesFiltered;
  timesFiltered.clear();
  std::sort(times.begin(), times.end()); // sort vector
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

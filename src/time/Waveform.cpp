#include "time/Waveform.hpp"
#include <algorithm>
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

void Waveform::clearTimeStepsStorage()
{
  // see https://stackoverflow.com/a/8234813, erase elements from map we are iterating over
  for (auto timeStep = _timeStepsStorage.begin(); timeStep != _timeStepsStorage.end();) {
    if (timeStep->first > 0.0) {
      timeStep = _timeStepsStorage.erase(timeStep);
    } else {
      ++timeStep;
    }
  }
}

int Waveform::maxNumberOfStoredWindows()
{
  PRECICE_ASSERT(_storageIsInitialized);
  if (_interpolationOrder > 0) {
    return _interpolationOrder;
  } else { // also for zeroth order interpolation we store one window (this window)
    return 1;
  }
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
  } else {                                // reached end of window and trying to write new data from next window. Clearing window first.
    auto startValues = this->sample(0.0); // use value at beginning of this window (= end of last window)
    this->clearTimeStepsStorage();        // clear storage for data with dt > 0.0 before redoing this window
  }
  PRECICE_ASSERT(values.size() == this->valuesSize(), values.size(), this->valuesSize());
  this->_timeStepsStorage[dt] = Eigen::VectorXd(values);
}

Eigen::VectorXd Waveform::sample(double normalizedDt)
{
  PRECICE_ASSERT(_storageIsInitialized);
  PRECICE_ASSERT(normalizedDt >= 0, "Sampling outside of valid range!");
  PRECICE_ASSERT(normalizedDt <= 1, "Sampling outside of valid range!");

  const int usedOrder = computeUsedOrder(_interpolationOrder, _numberOfStoredSamples);

  PRECICE_ASSERT(maxStoredDt() == 1.0); // sampling is only allowed, if a window is complete.

  // @todo do we need to explicitly differentiate between piecewise and non-piecewise interplation below?
  if (usedOrder == 0) {
    // constant interpolation = just use sample at the end of the window: x(dt) = x^t
    // At beginning of window use result from last window x(0) = x^(t-1)
    return Eigen::VectorXd(this->_timeStepsStorage[findTimeAfter(normalizedDt)]);
  }
  Eigen::VectorXd interpolatedValue;
  if (usedOrder == 1) {
    // linear interpolation inside window: x(dt) = dt * x^t + (1-dt) * x^(t-1)
    PRECICE_ASSERT(_numberOfStoredSamples > 1);
    interpolatedValue = this->_timeStepsStorage[1.0] * normalizedDt;        // = dt * x^t
    interpolatedValue += this->_timeStepsStorage[0.0] * (1 - normalizedDt); // = dt * x^t + (1-dt) * x^(t-1)
    return interpolatedValue;
  }
  PRECICE_ASSERT(usedOrder == 2);
  // quadratic interpolation inside window: x(dt) = x^t * (dt^2 + dt)/2 + x^(t-1) * (1-dt^2)+ x^(t-2) * (dt^2-dt)/2
  interpolatedValue = this->_timeStepsStorage[1.0] * (normalizedDt + 1) * normalizedDt * 0.5;
  interpolatedValue += this->_timeStepsStorage[0.0] * (1 - normalizedDt * normalizedDt);
  interpolatedValue += this->_timeStepsStorage[-1.0] * (normalizedDt - 1) * normalizedDt * 0.5;
  return interpolatedValue;
}

void Waveform::moveToNextWindow()
{
  PRECICE_ASSERT(_storageIsInitialized);
  PRECICE_ASSERT(maxNumberOfStoredWindows() <= 2); // other options are currently not implemented or supported.

  if (maxNumberOfStoredWindows() == 2) {
    _numberOfStoredSamples        = 3;
    this->_timeStepsStorage[-1.0] = Eigen::VectorXd(this->_timeStepsStorage[0.0]); // store values from past windows
  }
  auto initialGuess = this->sample(1.0); // use value at end of window as initial guess for next
  this->clearTimeStepsStorage();         // clear storage before entering next window
  _timeStepsStorage[0.0] = Eigen::VectorXd(initialGuess);
  _timeStepsStorage[1.0] = Eigen::VectorXd(initialGuess); // initial guess is always constant extrapolation
}

int Waveform::valuesSize()
{
  PRECICE_ASSERT(_storageIsInitialized);
  return _timeStepsStorage[0.0].size();
}

int Waveform::computeUsedOrder(int requestedOrder, int numberOfAvailableSamples)
{
  int usedOrder = -1;
  if (requestedOrder == 0) {
    usedOrder = 0;
  } else if (requestedOrder == 1) {
    if (numberOfAvailableSamples < 2) {
      usedOrder = 0;
    } else {
      usedOrder = 1;
    }
  } else if (requestedOrder == 2) {
    if (numberOfAvailableSamples < 2) {
      usedOrder = 0;
    } else if (numberOfAvailableSamples < 3) {
      usedOrder = 1;
    } else {
      usedOrder = 2;
    }
  } else {
    PRECICE_ASSERT(false);
  }
  return usedOrder;
}

// @todo improve efficiency by using some ordered list to store times.
double Waveform::findTimeBefore(double normalizedDt)
{
  PRECICE_ASSERT(0.0 <= normalizedDt);
  PRECICE_ASSERT(normalizedDt <= 1.0);

  double timeBefore = 0.0;

  for (auto timeStep : _timeStepsStorage) {
    if (timeBefore <= timeStep.first && timeStep.first <= normalizedDt) { // current timeStep is before normalizedDt and later than current timeBefore
      timeBefore = timeStep.first;
    }
  }

  return timeBefore;
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

#include "time/Waveform.hpp"
#include <algorithm>
#include "cplscheme/CouplingScheme.hpp"
#include "logging/LogMacros.hpp"
#include "time/Time.hpp"
#include "utils/EigenHelperFunctions.hpp"

#include <eigen3/unsupported/Eigen/Splines>
#include <cmath>

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

// helper function to compute x(t) from given data (x0,t0), (x1,t1), ..., (xn,tn) via third degree B-spline interpolation (implemented using Eigen). Alternatives are cubic spline interpolation (see https://en.wikipedia.org/wiki/Spline_interpolation#Algorithm_to_find_the_interpolating_cubic_spline), but this is not offered by Eigen or Boost.
Eigen::VectorXd bSplineInterpolationAt(double t, Eigen::VectorXd ts, Eigen::MatrixXd xs, int splineDegree)
{
  std::cout << "use degree:" << splineDegree << std::endl;
  // organize data in columns. Each column represents one sample in time.
  PRECICE_ASSERT(xs.cols() == ts.size());
  const int ndofs = xs.rows(); // number of dofs. Each dof needs it's own interpolant.

  Eigen::VectorXd interpolated(ndofs);


  const int splineDimension = 2;

  for (int i = 0; i < ndofs; i++) {
    // need to provide 2d coordinates (t,x) to make all points unique (0,0,1) is forbidden ((0,0), (0.5,0), (1,1)) is allowed
    Eigen::MatrixXd data(2, xs.cols());

    int method = 1;
    // only needed for method (3)
    Eigen::VectorXd chordLength(xs.cols()-1);
    double chordLengthSum = 0;

    for (int j = 0; j < xs.cols(); j++) {
      data.row(0)[j] = ts[j];
      data.row(1)[j] = xs.row(i)[j];
      if(method == 3) {
        if(j>0){
          double dts = ts[j]-ts[j-1];
          double dxs = xs.row(i)[j] - xs.row(i)[j-1];
          chordLength[j-1] = sqrt(dts*dts + dxs*dxs);
          chordLengthSum += chordLength[j-1];
        }
      }
    }

    // create know vector
    Eigen::VectorXd knots(ts.size());
    double          tmin = ts[0];
    double          tmax = ts[ts.size() - 1];
    for (int j = 0; j < xs.cols(); j++) {
      if(method == 1) {
        knots[j] = (ts[j] - tmin) / (tmax - tmin); // scale knots to interval [0,1]. Use ts as reference point.
      }
      if(method == 2) {  // uniform (see https://pages.mtu.edu/%7Eshene/COURSES/cs3621/NOTES/INT-APP/PARA-uniform.html)
        knots[j] = ((double)j)/(xs.cols()-1);
      }
      if(method == 3) {  //  chord length (see https://pages.mtu.edu/%7Eshene/COURSES/cs3621/NOTES/INT-APP/PARA-chord-length.html)
        if (j==0) {
          knots[j] = 0;
        } else if (j==xs.cols()-1) {
          knots[j] = 1;
        } else {
          knots[j] = 0;
        }
        for (int k = 0; k < j; k++) {
          knots[j] += chordLength[k];
        }
        knots[j] /= chordLengthSum;
      }
    }

    double tScaled  = (t - tmin) / (tmax - tmin); // scale sampling time to interval [0,1];
    auto   spline   = Eigen::SplineFitting<Eigen::Spline<double, splineDimension>>::Interpolate(data, splineDegree, knots);

    std::cout << "(t,x) = (" << spline(tScaled)[0] << "," << spline(tScaled)[1] << ")" << std::endl;

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

  bool   usePastWindow = true;
  double low           = 0.0;
  double high          = 1.0;
  if (usePastWindow) { // this just helps us to use quadratic interpolation without subcycling, as well. Use-Case unclear.
    low = -1.0;  // If we remove this, we can also assume that all our samples live on t in [0,1]. This simplifies the knot vector in the BSpline.
  }
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
  PRECICE_ASSERT(maxNumberOfStoredWindows() <= 2); // other options are currently not implemented or supported.

  if (maxNumberOfStoredWindows() == 2) {
    _numberOfStoredSamples = 3;
    // @TODO do we really want to store values from old windows? Use-case unclear and initialization difficult.
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
  // @todo simplify?
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
  } else if (requestedOrder == 3) {
    if (numberOfAvailableSamples < 2) {
      usedOrder = 0;
    } else if (numberOfAvailableSamples < 3) {
      usedOrder = 1;
    } else if (numberOfAvailableSamples < 4) {
      usedOrder = 2;
    } else {
      usedOrder = 3;
    }
  } else {
    PRECICE_ASSERT(false);
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

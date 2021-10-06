#include "time/Waveform.hpp"
#include <algorithm>
#include "logging/LogMacros.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice {
namespace time {

const int Waveform::UNDEFINED_INTERPOLATION_ORDER = -1;

Waveform::Waveform(
    const int initializedNumberOfData,
    const int extrapolationOrder,
    const int interpolationOrder)
    : _extrapolationOrder(extrapolationOrder)
{
  /**
     * Reserve storage depending on required extrapolation order. Extrapolation happens in-place. Therefore, for zeroth
     * order extrapolation we need one column (to read from and write to), for first order two, for second order three. 
     * Note that extrapolationOrder = 0 is an exception, since we want to always work with at least two samples. One at
     * the beginning and one at the end of the time window. Therefore, we use 2 samples for zeroth and first order
     * extrapolation.
     */
  PRECICE_ASSERT(extrapolationOrder >= 0);
  PRECICE_ASSERT(interpolationOrder >= 0 || interpolationOrder == UNDEFINED_INTERPOLATION_ORDER);
  const int initializedNumberOfSamples = std::max({2, _extrapolationOrder + 1, interpolationOrder + 1});
  _timeWindows                         = Eigen::MatrixXd::Zero(initializedNumberOfData, initializedNumberOfSamples);
  _numberOfValidSamples                = 1; // we assume that upon creation the first sample is always valid.
  PRECICE_ASSERT(numberOfSamples() == initializedNumberOfSamples);
  PRECICE_ASSERT(numberOfData() == initializedNumberOfData);
}

void Waveform::store(const Eigen::VectorXd &data)
{
  int columnID = 0;
  PRECICE_ASSERT(_timeWindows.cols() > columnID, numberOfSamples(), columnID);
  PRECICE_ASSERT(data.size() == numberOfData(), data.size(), numberOfData());
  this->_timeWindows.col(columnID) = data;
}

Eigen::VectorXd Waveform::sample(double normalizedDt, int order)
{
  PRECICE_ASSERT(normalizedDt >= 0, "Sampling outside of valid range!");
  PRECICE_ASSERT(normalizedDt <= 1, "Sampling outside of valid range!");
  // @ todo: Add more logic here? Set order in constructor; keep track of time windows inside class. See https://github.com/precice/precice/pull/1004/files?file-filters%5B%5D=.cmake&file-filters%5B%5D=.hpp#r642223767.
  return this->interpolateData(order, normalizedDt);
}

void Waveform::moveToNextWindow()
{
  auto initialGuess = extrapolateData();
  utils::shiftSetFirst(this->_timeWindows, initialGuess); // archive old samples and store initial guess
  if (_numberOfValidSamples < numberOfSamples()) {        // together with the initial guess the number of valid samples increases
    _numberOfValidSamples++;
  }
}

int Waveform::numberOfSamples()
{
  return _timeWindows.cols();
}

int Waveform::numberOfValidSamples()
{
  return _numberOfValidSamples;
}

int Waveform::numberOfData()
{
  return _timeWindows.rows();
}

const Eigen::MatrixXd &Waveform::lastTimeWindows()
{
  return _timeWindows;
}

/**
 * @brief Computes which order may be used for extrapolation or interpolation.
 * 
 * Order of extrapolation or interpolation is determined by number of valid samples and maximum order defined by the user.
 * Example: If only two samples are available, the maximum order we may use is 1, even if the user demands order 2.
 *
 * @param requestedOrder Order requested by the user.
 * @param numberOfAvailableSamples Samples available for extrapolation or interpolation.
 * @return Order that may be used.
 */
static int computeUsedOrder(int requestedOrder, int numberOfAvailableSamples)
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

Eigen::VectorXd Waveform::extrapolateData()
{
  const int usedOrder = computeUsedOrder(_extrapolationOrder, _numberOfValidSamples);

  if (usedOrder == 0) {
    PRECICE_ASSERT(this->numberOfSamples() > 0);
    return this->_timeWindows.col(0);
  }
  Eigen::VectorXd extrapolatedValue;
  if (usedOrder == 1) { //timesteps is increased before extrapolate is called
    PRECICE_DEBUG("Performing first order extrapolation");
    PRECICE_ASSERT(this->numberOfSamples() > 1);
    extrapolatedValue = this->_timeWindows.col(0) * 2.0; // = 2*x^t
    extrapolatedValue -= this->_timeWindows.col(1);      // = 2*x^t - x^(t-1)
    return extrapolatedValue;
  }
  PRECICE_ASSERT(usedOrder == 2);
  // uses formula given in https://doi.org/10.1016/j.compstruc.2008.11.013, p.796, Algorithm line 1
  PRECICE_DEBUG("Performing second order extrapolation");
  PRECICE_ASSERT(this->numberOfSamples() > 2);
  extrapolatedValue = this->_timeWindows.col(0) * 2.5;  // = 2.5*x^t
  extrapolatedValue -= this->_timeWindows.col(1) * 2.0; // = 2.5*x^t - 2*x^(t-1)
  extrapolatedValue += this->_timeWindows.col(2) * 0.5; // = 2.5*x^t - 2*x^(t-1) + 0.5*x^(t-2)
  return extrapolatedValue;
}

Eigen::VectorXd Waveform::interpolateData(const int order, const double normalizedDt)
{
  const int usedOrder = computeUsedOrder(order, _numberOfValidSamples);

  PRECICE_ASSERT(normalizedDt >= 0, "Sampling outside of valid range!");
  PRECICE_ASSERT(normalizedDt <= 1, "Sampling outside of valid range!");
  if (usedOrder == 0) {
    // constant interpolation = just use sample at the end of the window: x(dt) = x^t
    PRECICE_ASSERT(this->numberOfSamples() > 0);
    return this->_timeWindows.col(0);
  }
  Eigen::VectorXd interpolatedValue;
  if (usedOrder == 1) {
    // linear interpolation inside window: x(dt) = dt * x^t + (1-dt) * x^(t-1)
    PRECICE_ASSERT(this->numberOfSamples() > 1);
    interpolatedValue = this->_timeWindows.col(0) * normalizedDt;        // = dt * x^t
    interpolatedValue += this->_timeWindows.col(1) * (1 - normalizedDt); // = dt * x^t + (1-dt) * x^(t-1)
    return interpolatedValue;
  }
  PRECICE_ASSERT(usedOrder == 2);
  // quadratic interpolation inside window: x(dt) = x^t * (dt^2 + dt)/2 + x^(t-1) * (1-dt^2)+ x^(t-2) * (dt^2-dt)/2
  interpolatedValue = this->_timeWindows.col(0) * (normalizedDt + 1) * normalizedDt * 0.5;
  interpolatedValue += this->_timeWindows.col(1) * (1 - normalizedDt * normalizedDt);
  interpolatedValue += this->_timeWindows.col(2) * (normalizedDt - 1) * normalizedDt * 0.5;
  return interpolatedValue;
}

} // namespace time
} // namespace precice

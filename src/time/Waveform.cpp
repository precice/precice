#include "time/Waveform.hpp"
#include <algorithm>
#include "logging/LogMacros.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice {
namespace time {

const int Waveform::UNDEFINED_INTERPOLATION_ORDER = -1;

Waveform::Waveform(
    const int valuesSize,
    const int extrapolationOrder,
    const int interpolationOrder)
    : _extrapolationOrder(extrapolationOrder), _interpolationOrder(interpolationOrder)
{
  /**
     * Reserve storage depending on required extrapolation order. Extrapolation happens in-place. Therefore, for zeroth
     * order extrapolation we need one column (to read from and write to), for first order two, for second order three. 
     * Note that extrapolationOrder = 0 is an exception, since we want to always work with at least two samples. One at
     * the beginning and one at the end of the time window. Therefore, we use 2 samples for zeroth and first order
     * extrapolation.
     */
  PRECICE_ASSERT(_extrapolationOrder >= 0);
  PRECICE_ASSERT(_interpolationOrder >= 0 || _interpolationOrder == UNDEFINED_INTERPOLATION_ORDER);
  const int sampleStorageSize = std::max({2, _extrapolationOrder + 1, interpolationOrder + 1});
  _timeWindowsStorage         = Eigen::MatrixXd::Zero(valuesSize, sampleStorageSize);
  _numberOfStoredSamples      = 1; // the first sample is automatically initialized as zero and stored.
  PRECICE_ASSERT(this->sizeOfSampleStorage() == sampleStorageSize);
  PRECICE_ASSERT(this->valuesSize() == valuesSize);
}

void Waveform::store(const Eigen::VectorXd &values)
{
  int columnID = 0;
  PRECICE_ASSERT(_timeWindowsStorage.cols() > columnID, sizeOfSampleStorage(), columnID);
  PRECICE_ASSERT(values.size() == this->valuesSize(), values.size(), this->valuesSize());
  this->_timeWindowsStorage.col(columnID) = values;
}

Eigen::VectorXd Waveform::sample(double normalizedDt)
{
  PRECICE_ASSERT(normalizedDt >= 0, "Sampling outside of valid range!");
  PRECICE_ASSERT(normalizedDt <= 1, "Sampling outside of valid range!");
  PRECICE_ASSERT(_interpolationOrder != UNDEFINED_INTERPOLATION_ORDER, "Sampling is only allowed, if Waveform is configured correspondingly");
  return this->interpolate(normalizedDt);
}

const Eigen::VectorXd Waveform::getInitialGuess()
{
  return _timeWindowsStorage.col(0);
}

void Waveform::moveToNextWindow()
{
  auto initialGuess = extrapolate();
  utils::shiftSetFirst(this->_timeWindowsStorage, initialGuess); // archive old samples and store initial guess
  if (_numberOfStoredSamples < sizeOfSampleStorage()) {          // together with the initial guess the number of stored samples increases
    _numberOfStoredSamples++;
  }
}

const Eigen::VectorXd Waveform::getInitialGuess()
{
  return _timeWindowsStorage.col(0);
}

int Waveform::sizeOfSampleStorage()
{
  return _timeWindowsStorage.cols();
}

int Waveform::valuesSize()
{
  return _timeWindowsStorage.rows();
}

/**
 * @brief Computes which order may be used for extrapolation or interpolation.
 * 
 * Order of extrapolation or interpolation is determined by number of stored samples and maximum order defined by the user.
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

Eigen::VectorXd Waveform::extrapolate()
{
  const int usedOrder = computeUsedOrder(_extrapolationOrder, _numberOfStoredSamples);

  if (usedOrder == 0) {
    PRECICE_ASSERT(_numberOfStoredSamples > 0);
    return _timeWindowsStorage.col(0);
  }
  Eigen::VectorXd extrapolatedValue;
  if (usedOrder == 1) { //timesteps is increased before extrapolate is called
    PRECICE_DEBUG("Performing first order extrapolation");
    PRECICE_ASSERT(_numberOfStoredSamples > 1);
    extrapolatedValue = _timeWindowsStorage.col(0) * 2.0; // = 2*x^t
    extrapolatedValue -= _timeWindowsStorage.col(1);      // = 2*x^t - x^(t-1)
    return extrapolatedValue;
  }
  PRECICE_ASSERT(usedOrder == 2);
  // uses formula given in https://doi.org/10.1016/j.compstruc.2008.11.013, p.796, Algorithm line 1
  PRECICE_DEBUG("Performing second order extrapolation");
  PRECICE_ASSERT(_numberOfStoredSamples > 2);
  extrapolatedValue = _timeWindowsStorage.col(0) * 2.5;  // = 2.5*x^t
  extrapolatedValue -= _timeWindowsStorage.col(1) * 2.0; // = 2.5*x^t - 2*x^(t-1)
  extrapolatedValue += _timeWindowsStorage.col(2) * 0.5; // = 2.5*x^t - 2*x^(t-1) + 0.5*x^(t-2)
  return extrapolatedValue;
}

Eigen::VectorXd Waveform::interpolate(const double normalizedDt)
{
  const int usedOrder = computeUsedOrder(_interpolationOrder, _numberOfStoredSamples);

  PRECICE_ASSERT(normalizedDt >= 0, "Sampling outside of valid range!");
  PRECICE_ASSERT(normalizedDt <= 1, "Sampling outside of valid range!");
  if (usedOrder == 0) {
    // constant interpolation = just use sample at the end of the window: x(dt) = x^t
    PRECICE_ASSERT(_numberOfStoredSamples > 0);
    return this->_timeWindowsStorage.col(0);
  }
  Eigen::VectorXd interpolatedValue;
  if (usedOrder == 1) {
    // linear interpolation inside window: x(dt) = dt * x^t + (1-dt) * x^(t-1)
    PRECICE_ASSERT(_numberOfStoredSamples > 1);
    interpolatedValue = this->_timeWindowsStorage.col(0) * normalizedDt;        // = dt * x^t
    interpolatedValue += this->_timeWindowsStorage.col(1) * (1 - normalizedDt); // = dt * x^t + (1-dt) * x^(t-1)
    return interpolatedValue;
  }
  PRECICE_ASSERT(usedOrder == 2);
  // quadratic interpolation inside window: x(dt) = x^t * (dt^2 + dt)/2 + x^(t-1) * (1-dt^2)+ x^(t-2) * (dt^2-dt)/2
  interpolatedValue = this->_timeWindowsStorage.col(0) * (normalizedDt + 1) * normalizedDt * 0.5;
  interpolatedValue += this->_timeWindowsStorage.col(1) * (1 - normalizedDt * normalizedDt);
  interpolatedValue += this->_timeWindowsStorage.col(2) * (normalizedDt - 1) * normalizedDt * 0.5;
  return interpolatedValue;
}

} // namespace time
} // namespace precice

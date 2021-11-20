#include "time/Waveform.hpp"
#include <algorithm>
#include "cplscheme/CouplingScheme.hpp"
#include "logging/LogMacros.hpp"
#include "time/Time.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice {
namespace time {

Waveform::Waveform(
    const int extrapolationOrder,
    const int interpolationOrder)
    : _extrapolationOrder(extrapolationOrder), _interpolationOrder(interpolationOrder)
{
  PRECICE_ASSERT(not _storageIsInitialized);
}

void Waveform::initialize(
    const int valuesSize)
{
  int sampleStorageSize  = std::max({2, _extrapolationOrder + 1, _interpolationOrder + 1});
  _timeWindowsStorage    = Eigen::MatrixXd::Zero(valuesSize, sampleStorageSize);
  _numberOfStoredSamples = 1; // the first sample is automatically initialized as zero and stored.
  _storageIsInitialized  = true;
  PRECICE_ASSERT(this->sizeOfSampleStorage() == sampleStorageSize);
  PRECICE_ASSERT(this->valuesSize() == valuesSize);
}

void Waveform::resizeData(int newValuesSize)
{
  _timeWindowsStorage = Eigen::MatrixXd::Zero(newValuesSize, sizeOfSampleStorage());
  PRECICE_ASSERT(valuesSize() == newValuesSize);
}

void Waveform::store(const Eigen::VectorXd &values)
{
  PRECICE_ASSERT(_storageIsInitialized);
  int columnID = 0;
  this->storeAt(values, columnID);
}

void Waveform::storeAt(const Eigen::VectorXd values, int columnID)
{
  PRECICE_ASSERT(_timeWindowsStorage.cols() > columnID, sizeOfSampleStorage(), columnID);
  PRECICE_ASSERT(values.size() == this->valuesSize(), values.size(), this->valuesSize());
  this->_timeWindowsStorage.col(columnID) = values;
}

Eigen::VectorXd Waveform::sample(double normalizedDt)
{
  PRECICE_ASSERT(_storageIsInitialized);
  PRECICE_ASSERT(normalizedDt >= 0, "Sampling outside of valid range!");
  PRECICE_ASSERT(normalizedDt <= 1, "Sampling outside of valid range!");
  PRECICE_ASSERT(_interpolationOrder != Time::UNDEFINED_INTERPOLATION_ORDER, "Sampling is only allowed, if Waveform is configured correspondingly");
  return this->interpolate(normalizedDt);
}

void Waveform::moveToNextWindow()
{
  PRECICE_ASSERT(_storageIsInitialized);
  if (_extrapolationOrder != Time::UNDEFINED_EXTRAPOLATION_ORDER) { // @todo: we should define the extrapolation data for all data (if the waveform is owned by the mesh::Data), another argument for configuring extrapolation order there.
    auto initialGuess = extrapolate();
    utils::shiftSetFirst(this->_timeWindowsStorage, initialGuess); // archive old samples and store initial guess
    if (_numberOfStoredSamples < sizeOfSampleStorage()) {          // together with the initial guess the number of stored samples increases
      _numberOfStoredSamples++;
    }
  }
}

const Eigen::VectorXd Waveform::getInitialGuess()
{
  PRECICE_ASSERT(_storageIsInitialized);
  return _timeWindowsStorage.col(0);
}

int Waveform::sizeOfSampleStorage()
{
  PRECICE_ASSERT(_storageIsInitialized);
  return _timeWindowsStorage.cols();
}

Eigen::VectorXd Waveform::getSample(int sampleID)
{
  //PRECICE_ASSERT(sampleID < _numberOfStoredSamples);  // @todo use this stricted assertion?
  PRECICE_ASSERT(sampleID < sizeOfSampleStorage());
  return _timeWindowsStorage.col(sampleID);
}

void Waveform::setExtrapolationOrder(int extrapolationOrder)
{
  _extrapolationOrder = extrapolationOrder;
}

int Waveform::valuesSize()
{
  PRECICE_ASSERT(_storageIsInitialized);
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
  PRECICE_ASSERT(_storageIsInitialized);

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

#include "time/Waveform.hpp"
#include <algorithm>
#include "logging/LogMacros.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice {
namespace time {

Waveform::Waveform(
    const int extrapolationOrder)
    : _extrapolationOrder(extrapolationOrder)
{
  PRECICE_ASSERT(not _storageIsInitialized);
}

void Waveform::initializeData(
    const int valuesSize)
{
  int sampleStorageSize  = std::max({2, _extrapolationOrder + 1});
  _timeWindowsStorage    = Eigen::MatrixXd::Zero(valuesSize, sampleStorageSize);
  _numberOfStoredSamples = 1; // the first sample is automatically initialized as zero and stored.
  _storageIsInitialized  = true;
  PRECICE_ASSERT(this->sizeOfSampleStorage() == sampleStorageSize);
  PRECICE_ASSERT(this->valuesSize() == valuesSize);
}

void Waveform::store(const Eigen::VectorXd &values)
{
  PRECICE_ASSERT(_storageIsInitialized);
  int columnID = 0;
  PRECICE_ASSERT(_timeWindowsStorage.cols() > columnID, sizeOfSampleStorage(), columnID);
  PRECICE_ASSERT(values.size() == this->valuesSize(), values.size(), this->valuesSize());
  this->_timeWindowsStorage.col(columnID) = values;
}

void Waveform::moveToNextWindow()
{
  PRECICE_ASSERT(_storageIsInitialized);
  auto initialGuess = extrapolate();
  utils::shiftSetFirst(this->_timeWindowsStorage, initialGuess); // archive old samples and store initial guess
  if (_numberOfStoredSamples < sizeOfSampleStorage()) {          // together with the initial guess the number of stored samples increases
    _numberOfStoredSamples++;
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

} // namespace time
} // namespace precice

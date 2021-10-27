#include "time/Waveform.hpp"
#include <algorithm>
#include "logging/LogMacros.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice {
namespace time {

Waveform::Waveform(
    const int dataCount,
    const int extrapolationOrder)
    : _extrapolationOrder(extrapolationOrder)
{
  /**
     * Reserve storage depending on required extrapolation order. Extrapolation happens in-place. Therefore, for zeroth
     * order extrapolation we need one column (to read from and write to), for first order two, for second order three. 
     * Note that extrapolationOrder = 0 is an exception, since we want to always work with at least two samples. One at
     * the beginning and one at the end of the time window. Therefore, we use 2 samples for zeroth and first order
     * extrapolation.
     */
  int sampleStorageSize = std::max({2, _extrapolationOrder + 1});
  _timeWindowsStorage                   = Eigen::MatrixXd::Zero(dataCount, sampleStorageSize);
  _numberOfStoredSamples          = 1; // we assume that upon creation the first sample is always valid.
  PRECICE_ASSERT(this->sizeOfSampleStorage() == sampleStorageSize);
  PRECICE_ASSERT(this->dataCount() == dataCount);
}

void Waveform::store(const Eigen::VectorXd &data)
{
  int columnID = 0;
  PRECICE_ASSERT(_timeWindowsStorage.cols() > columnID, sizeOfSampleStorage(), columnID);
  PRECICE_ASSERT(data.size() == dataCount(), data.size(), dataCount());
  this->_timeWindowsStorage.col(columnID) = data;
}

void Waveform::moveToNextWindow()
{
  auto initialGuess = extrapolateData();
  utils::shiftSetFirst(this->_timeWindowsStorage, initialGuess); // archive old samples and store initial guess
  if (_numberOfStoredSamples < sizeOfSampleStorage()) {     // together with the initial guess the number of valid samples increases
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

int Waveform::dataCount()
{
  return _timeWindowsStorage.rows();
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

} // namespace time
} // namespace precice

#include "Extrapolation.hpp"
#include <algorithm>
#include "cplscheme/CouplingScheme.hpp"
#include "logging/LogMacros.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice::cplscheme::impl {

Extrapolation::Extrapolation(
    const int extrapolationOrder)
    : _extrapolationOrder(extrapolationOrder)
{
  PRECICE_ASSERT(not _storageIsInitialized);
}

void Extrapolation::initialize(
    const int valuesSize)
{
  int sampleStorageSize  = std::max({_extrapolationOrder + 1});
  _timeWindowsStorage    = Eigen::MatrixXd::Zero(valuesSize, sampleStorageSize);
  _numberOfStoredSamples = 1; // the first sample is automatically initialized as zero and stored.
  _storageIsInitialized  = true;
  PRECICE_ASSERT(this->sizeOfSampleStorage() == sampleStorageSize);
  PRECICE_ASSERT(this->valuesSize() == valuesSize);
}

void Extrapolation::store(const Eigen::VectorXd &values)
{
  PRECICE_ASSERT(_storageIsInitialized);
  if (values.size() != this->valuesSize()) { // @todo work-around, but becomes unnecessary with https://github.com/precice/precice/pull/1639
    initialize(values.size());
  }
  PRECICE_ASSERT(values.size() == this->valuesSize(), values.size(), this->valuesSize());
  this->_timeWindowsStorage.col(0) = values;
}

void Extrapolation::moveToNextWindow()
{
  PRECICE_ASSERT(_storageIsInitialized);
  auto initialGuess = extrapolate();
  utils::shiftSetFirst(this->_timeWindowsStorage, initialGuess); // archive old samples and store initial guess
  if (_numberOfStoredSamples < sizeOfSampleStorage()) {          // together with the initial guess the number of stored samples increases
    _numberOfStoredSamples++;
  }
}

const Eigen::VectorXd Extrapolation::getInitialGuess()
{
  PRECICE_ASSERT(_storageIsInitialized);
  return _timeWindowsStorage.col(0);
}

int Extrapolation::sizeOfSampleStorage()
{
  PRECICE_ASSERT(_storageIsInitialized);
  return _timeWindowsStorage.cols();
}

int Extrapolation::valuesSize()
{
  PRECICE_ASSERT(_storageIsInitialized);
  return _timeWindowsStorage.rows();
}

/**
 * @brief Computes which order may be used for extrapolation.
 *
 * Order of extrapolation is determined by number of stored samples and maximum order defined by the user.
 * Example: If only two samples are available, the maximum order we may use is 1, even if the user demands order 2.
 *
 * @param requestedOrder Order requested by the user.
 * @param numberOfAvailableSamples Samples available for extrapolation.
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

Eigen::VectorXd Extrapolation::extrapolate()
{
  const int usedOrder = computeUsedOrder(_extrapolationOrder, _numberOfStoredSamples);
  PRECICE_ASSERT(_storageIsInitialized);

  if (usedOrder == 0) {
    PRECICE_ASSERT(_numberOfStoredSamples > 0);
    return _timeWindowsStorage.col(0);
  }
  Eigen::VectorXd extrapolatedValue;
  PRECICE_ASSERT(usedOrder == 1);
  PRECICE_DEBUG("Performing first order extrapolation");
  PRECICE_ASSERT(_numberOfStoredSamples > 1);
  extrapolatedValue = _timeWindowsStorage.col(0) * 2.0; // = 2*x^t
  extrapolatedValue -= _timeWindowsStorage.col(1);      // = 2*x^t - x^(t-1)
  return extrapolatedValue;
}

} // namespace precice::cplscheme::impl

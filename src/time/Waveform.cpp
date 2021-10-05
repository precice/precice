#include "time/Waveform.hpp"
#include <algorithm>
#include "logging/LogMacros.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice {
namespace time {

Waveform::Waveform(
    int initializedNumberOfData,
    int extrapolationOrder)
    : _extrapolationOrder(extrapolationOrder)
{
  /**
     * Reserve storage depending on required extrapolation order. Extrapolation happens in-place. Therefore, for zeroth
     * order extrapolation we need one column (to read from and write to), for first order two, for second order three. 
     * Note that extrapolationOrder = 0 is an exception, since we want to always work with at least two samples. One at
     * the beginning and one at the end of the time window. Therefore, we use 2 samples for zeroth and first order
     * extrapolation.
     */
  int initializedNumberOfSamples = std::max({2, _extrapolationOrder + 1});
  _timeWindows                   = Eigen::MatrixXd::Zero(initializedNumberOfData, initializedNumberOfSamples);
  _numberOfValidSamples          = 1; // we assume that upon creation the first sample is always valid.
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

int Waveform::computeUsedOrder(int order)
{
  int usedOrder = -1;
  if (_extrapolationOrder == 0) {
    usedOrder = 0;
  } else if (_extrapolationOrder == 1) {
    if (_numberOfValidSamples < 2) {
      usedOrder = 0;
    } else {
      usedOrder = 1;
    }
  } else if (_extrapolationOrder == 2) {
    if (_numberOfValidSamples < 2) {
      usedOrder = 0;
    } else if (_numberOfValidSamples < 3) {
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
  int usedOrder = computeUsedOrder(_extrapolationOrder);

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

} // namespace time
} // namespace precice

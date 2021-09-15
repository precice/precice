#include "time/Waveform.hpp"
#include <algorithm>
#include "logging/LogMacros.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice {
namespace time {

Waveform::Waveform(
    int initializedNumberOfData,
    int extrapolationOrder)
{
  /**
     * Reserve storage depending on required extrapolation order. Extrapolation happens in-place. Therefore, for zeroth
     * order extrapolation we need one column (to read from and write to), for first order two, for second order three. 
     * Note that extrapolationOrder = 0 is an exception, since we want to always work with at least two samples. One at
     * the beginning and one at the end of the time window. Therefore, we use 2 samples for zeroth and first order
     * extrapolation.
     */
  int initializedNumberOfSamples = std::max({2, extrapolationOrder + 1});
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

void Waveform::moveToNextWindow(int order)
{
  auto initialGuess = extrapolateData(order);
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

Eigen::VectorXd Waveform::extrapolateData(int order)
{
  Eigen::VectorXd extrapolatedValue;
  if ((order == 0) || (_numberOfValidSamples < 2 && order > 0)) {
    PRECICE_ASSERT(this->numberOfSamples() > 0);
    extrapolatedValue = this->_timeWindows.col(0);
  } else if ((order == 1) || (_numberOfValidSamples < 3 && order > 1)) { //timesteps is increased before extrapolate is called
    PRECICE_DEBUG("Performing first order extrapolation");
    PRECICE_ASSERT(this->numberOfSamples() > 1);
    extrapolatedValue = this->_timeWindows.col(0) * 2.0; // = 2*x^t
    extrapolatedValue -= this->_timeWindows.col(1);      // = 2*x^t - x^(t-1)
  } else if (order == 2) {
    // uses formula given in https://doi.org/10.1016/j.compstruc.2008.11.013, p.796, Algorithm line 1
    PRECICE_DEBUG("Performing second order extrapolation");
    PRECICE_ASSERT(this->numberOfSamples() > 2);
    extrapolatedValue = this->_timeWindows.col(0) * 2.5;  // = 2.5*x^t
    extrapolatedValue -= this->_timeWindows.col(1) * 2.0; // = 2.5*x^t - 2*x^(t-1)
    extrapolatedValue += this->_timeWindows.col(2) * 0.5; // = 2.5*x^t - 2*x^(t-1) + 0.5*x^(t-2)
  } else {
    PRECICE_ASSERT(false, "Extrapolation order is invalid.");
  }
  return extrapolatedValue;
}

} // namespace time
} // namespace precice

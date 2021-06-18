#include "time/Waveform.hpp"
#include <algorithm>
#include "logging/LogMacros.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice {
namespace time {

Waveform::Waveform(
    int numberOfData,
    int extrapolationOrder)
{
  /**
     * Reserve storage depending on required extrapolation order. Extrapolation happens in-place. Therefore, for zeroth
     * order extrapolation we need one column (to read from and write to), for first order two, for second order three. 
     * Note that extrapolationOrder = 0 is an exception, since we want to always work with at least two samples. One at
     * the beginning and one at the end of the time window. Therefore, we use 2 samples for zeroth and first order
     * extrapolation.
     */
  int numberOfSamples = std::max(2, extrapolationOrder + 1);
  _timeWindows        = Eigen::MatrixXd::Zero(numberOfData, numberOfSamples);
}

void Waveform::store(Eigen::VectorXd data)
{
  this->_timeWindows.col(0) = data;
}

void Waveform::moveToNextWindow(int timeWindows, int order)
{
  auto initialGuess = extrapolateData(order, timeWindows);
  utils::shiftSetFirst(this->_timeWindows, initialGuess);
}

const Eigen::MatrixXd &Waveform::lastTimeWindows()
{
  return _timeWindows;
}

Eigen::VectorXd Waveform::extrapolateData(int order, int timeWindows)
{
  Eigen::VectorXd extrapolatedValue;
  if ((order == 0) || (timeWindows < 2 && order > 0)) {
    PRECICE_ASSERT(this->_timeWindows.cols() > 0);
    extrapolatedValue = this->_timeWindows.col(0);
  } else if ((order == 1) || (timeWindows < 3 && order > 1)) { //timesteps is increased before extrapolate is called
    PRECICE_DEBUG("Performing first order extrapolation");
    PRECICE_ASSERT(this->_timeWindows.cols() > 1);
    extrapolatedValue = this->_timeWindows.col(0) * 2.0; // = 2*x^t
    extrapolatedValue -= this->_timeWindows.col(1);      // = 2*x^t - x^(t-1)
  } else if (order == 2) {
    PRECICE_DEBUG("Performing second order extrapolation");
    PRECICE_ASSERT(this->_timeWindows.cols() > 2);
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

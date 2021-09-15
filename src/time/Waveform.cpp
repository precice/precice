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

void Waveform::moveToNextWindow(int timeWindows, int order)
{
  auto initialGuess = extrapolateData(order, timeWindows);
  utils::shiftSetFirst(this->_timeWindows, initialGuess);
}

const Eigen::MatrixXd &Waveform::lastTimeWindows()
{
  return _timeWindows;
}

int Waveform::numberOfSamples()
{
  return _timeWindows.cols();
}

int Waveform::numberOfData()
{
  return _timeWindows.rows();
}

Eigen::VectorXd Waveform::extrapolateData(int order, int timeWindows)
{
  Eigen::VectorXd extrapolatedValue;
  if ((order == 0) || (timeWindows < 2 && order > 0)) {
    PRECICE_ASSERT(this->numberOfSamples() > 0);
    extrapolatedValue = this->_timeWindows.col(0);
  } else if ((order == 1) || (timeWindows < 3 && order > 1)) { //timesteps is increased before extrapolate is called
    PRECICE_DEBUG("Performing first order extrapolation");
    PRECICE_ASSERT(this->numberOfSamples() > 1);
    extrapolatedValue = this->_timeWindows.col(0) * 2.0; // = 2*x^t
    extrapolatedValue -= this->_timeWindows.col(1);      // = 2*x^t - x^(t-1)
    // see https://github.com/precice/precice/issues/1089 for derivation.
  } else if (order == 2) {
    PRECICE_DEBUG("Performing second order extrapolation");
    PRECICE_ASSERT(this->numberOfSamples() > 2);
    extrapolatedValue = this->_timeWindows.col(0) * 3;  // = 3*x^t
    extrapolatedValue -= this->_timeWindows.col(1) * 3; // = 3*x^t - 3*x^(t-1)
    extrapolatedValue += this->_timeWindows.col(2);     // = 3*x^t - 3*x^(t-1) + x^(t-2)
    // see https://github.com/precice/precice/issues/1089 for derivation.
  } else {
    PRECICE_ASSERT(false, "Extrapolation order is invalid.");
  }
  return extrapolatedValue;
}

} // namespace time
} // namespace precice

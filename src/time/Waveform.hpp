#pragma once

#include <Eigen/Core>
#include <algorithm>
#include "logging/LogMacros.hpp"
#include "logging/Logger.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice {
namespace time {

class Waveform {
public:
  /**
   * @brief Waveform object which stores data of current and past time windows for performing extrapolation.
   * @param numberOfData defines how many pieces of data one sample in time consists of
   * @param extrapolatioOrder defines the maximum extrapolation order supported by this Waveform and reserves storage correspondingly
   */
  Waveform(int numberOfData,
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

  /**
   * @brief Updates entry in _timeWindows corresponding to this window with given data
   * @param data new sample for this time window
   */
  void updateThisWindow(Eigen::VectorXd data)
  {
    this->_timeWindows.col(0) = data;
  }

  /**
   * @brief Called, when moving to the next time window. All entries in _timeWindows are shifted. The new entry is initialized as the value from the last window (= constant extrapolation)
   * @param timeWindows number of samples that are valid and may be used for extrapolation. Usually number of past time windows.
   */
  void moveToNextWindow(int timeWindows, int order = 0)
  {
    auto initialGuess = extrapolateData(order, timeWindows);
    utils::shiftSetFirst(this->_timeWindows, initialGuess);
  }

  /**
   * @brief getter for Eigen::MatrixXd containing data of current and past time windows. Each column represents a sample in time, with col(0)
   * being the current time window.
   */
  const Eigen::MatrixXd &lastTimeWindows()
  {
    return _timeWindows;
  }

private:
  /// Data values of time windows.
  Eigen::MatrixXd _timeWindows;

  mutable logging::Logger _log{"time::Waveform"};

  /**
   * @brief Extrapolates data _timeWindows using an extrapolation scheme of given order. 
   * 
   * If the order condition cannot be satisfied, since there are not enough samples available, the order is automatically reduced.
   * If order two is required, but only two samples are available, the extrapolation order is automatically reduced to one.
   * 
   * @param order Order of the extrapolation scheme to be used.
   * @param timeWindows number of valid samples.
   */
  Eigen::VectorXd extrapolateData(int order, int timeWindows)
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
};

} // namespace time
} // namespace precice

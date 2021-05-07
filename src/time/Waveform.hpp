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
  Waveform(int numberOfVertices,
           int numberOfSamples)
  {
    numberOfSamples = std::max(2, numberOfSamples); // work with at least two samples (beginning and end of time window)
    _timeWindows    = Eigen::MatrixXd::Zero(numberOfVertices, numberOfSamples);
  }

  /**
   * @brief Updates entry in _timeWindows corresponding to this window with given data
   */
  void updateThisWindow(Eigen::VectorXd data)
  {
    this->_timeWindows.col(0) = data;
  }

  /**
   * @brief Called, when moving to the next time window. All entries in _timeWindows are shifted. The new entry is initialized as the value from the last window (= constant extrapolation)
   */
  void moveToNextWindow(int timeWindows, int order = 0)
  {
    auto initialGuess = extrapolateData(order, timeWindows);
    utils::shiftSetFirst(this->_timeWindows, initialGuess);
  }

  const Eigen::MatrixXd &lastTimeWindows()
  {
    return _timeWindows;
  }

private:
  /// Data values of time windows.
  Eigen::MatrixXd _timeWindows;

  mutable logging::Logger _log{"time::Waveform"};

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

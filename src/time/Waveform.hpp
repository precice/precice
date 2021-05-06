#pragma once

#include <Eigen/Core>
#include <algorithm>
#include "utils/EigenHelperFunctions.hpp"

namespace precice {
namespace time {

class Waveform {
public:
  Waveform(int numberOfVertices,
           int numberOfSamples)
  {
    numberOfSamples  = std::max(2, numberOfSamples); // work with at least two samples (beginning and end of time window)
    _lastTimeWindows = Eigen::MatrixXd::Zero(numberOfVertices, numberOfSamples);
  }

  /**
   * @brief Updates entry in _lastTimeWindows corresponding to this window with given data
   */
  void updateThisWindow(Eigen::VectorXd data)
  {
    this->_lastTimeWindows.col(0) = data;
  }

  /**
   * @brief Called, when moving to the next time window. All entries in _lastTimeWindows are shifted. The new entry is initialized as the value from the last window (= constant extrapolation)
   */
  void moveToNextWindow()
  {
    // @todo try to merge this method with extrapolate data. Signature moveToNextWindow(int computedTimeWindows, int order = 0)
    auto initialGuess = this->_lastTimeWindows.col(0);
    // For extrapolation, treat the initial value as old time windows value
    utils::shiftSetFirst(this->_lastTimeWindows, initialGuess);
    PRECICE_ASSERT(initialGuess == extrapolateData(0, -1), "Initial guess is identical to constant extrapolation");
  }

  Eigen::VectorXd extrapolateData(int order, int timeWindows)
  {
    Eigen::VectorXd extrapolatedValue;
    if (order == 0) {
      PRECICE_ASSERT(this->_lastTimeWindows.cols() > 0);
      extrapolatedValue = this->_lastTimeWindows.col(1);
    } else if ((order == 1) || (timeWindows == 2 && order == 2)) { //timesteps is increased before extrapolate is called
      // PRECICE_INFO("Performing first order extrapolation");
      PRECICE_ASSERT(this->_lastTimeWindows.cols() > 1);
      extrapolatedValue = this->_lastTimeWindows.col(1) * 2.0; // = 2*x^t
      extrapolatedValue -= this->_lastTimeWindows.col(2);      // = 2*x^t - x^(t-1)
    } else if (order == 2) {
      // PRECICE_INFO("Performing second order extrapolation");
      PRECICE_ASSERT(this->_lastTimeWindows.cols() > 2);
      extrapolatedValue = this->_lastTimeWindows.col(1) * 2.5;  // = 2.5*x^t
      extrapolatedValue -= this->_lastTimeWindows.col(2) * 2.0; // = 2.5*x^t - 2*x^(t-1)
      extrapolatedValue += this->_lastTimeWindows.col(3) * 0.5; // = 2.5*x^t - 2*x^(t-1) + 0.5*x^(t-2)
    } else {
      PRECICE_ASSERT(false, "Extrapolation order is invalid.");
    }
    return extrapolatedValue;
  }

  const Eigen::MatrixXd &lastTimeWindows()
  {
    return _lastTimeWindows;
  }

private:
  /// Data values of previous time windows.
  Eigen::MatrixXd _lastTimeWindows;
};

} // namespace time
} // namespace precice
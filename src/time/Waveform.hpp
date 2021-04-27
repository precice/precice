#pragma once

#include <Eigen/Core>
#include "utils/EigenHelperFunctions.hpp"

namespace precice {
namespace time {

class Waveform {
public:
  void addNewWindowData(Eigen::VectorXd data)
  {
    // For extrapolation, treat the initial value as old time windows value
    utils::shiftSetFirst(this->_lastTimeWindows, data);
  }

  Eigen::VectorXd extrapolateData(int order, int timeWindows)
  {
    Eigen::VectorXd extrapolatedValue;
    if ((order == 1) || (timeWindows == 2 && order == 2)) { //timesteps is increased before extrapolate is called
      // PRECICE_INFO("Performing first order extrapolation");
      PRECICE_ASSERT(this->_lastTimeWindows.cols() > 1);
      extrapolatedValue = this->_lastTimeWindows.col(0) * 2.0; // = 2*x^t
      extrapolatedValue -= this->_lastTimeWindows.col(1);      // = 2*x^t - x^(t-1)
    } else if (order == 2) {
      // PRECICE_INFO("Performing second order extrapolation");
      PRECICE_ASSERT(this->_lastTimeWindows.cols() > 2);
      extrapolatedValue = this->_lastTimeWindows.col(0) * 2.5;  // = 2.5*x^t
      extrapolatedValue -= this->_lastTimeWindows.col(1) * 2.0; // = 2.5*x^t - 2*x^(t-1)
      extrapolatedValue += this->_lastTimeWindows.col(2) * 0.5; // = 2.5*x^t - 2*x^(t-1) + 0.5*x^(t-2)
    } else {
      PRECICE_ASSERT(false, "Extrapolation order is invalid.");
    }
    return extrapolatedValue;
  }

  void initialize(int numberOfVertices, int numberOfSamples)
  {
    _lastTimeWindows = Eigen::MatrixXd::Zero(numberOfVertices, numberOfSamples);
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
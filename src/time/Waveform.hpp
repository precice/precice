#pragma once

#include <Eigen/Core>
#include "logging/Logger.hpp"

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
           int extrapolationOrder);

  /**
   * @brief Updates entry in _timeWindows corresponding to this window with given data
   * @param data new sample for this time window
   */
  void store(const Eigen::VectorXd &data);

  /**
   * @brief Called, when moving to the next time window. All entries in _timeWindows are shifted. The new entry is initialized as the value from the last window (= constant extrapolation)
   * @param timeWindows number of samples that are valid and may be used for extrapolation. Usually number of past time windows.
   */
  void moveToNextWindow(int timeWindows, int order = 0);

  /**
   * @brief getter for Eigen::MatrixXd containing data of current and past time windows. Each column represents a sample in time, with col(0)
   * being the current time window.
   */
  const Eigen::MatrixXd &lastTimeWindows();

  /**
   * @brief returns number of samples in time stored by this waveform
   */
  int numberOfSamples();

  /**
   * @brief returns number of data per sample in time stored by this waveform
   */
  int numberOfData();

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
   * For linear extrapolation two equally spaced samples x^t and x^(t-1) are used to compute x^(t+1) using a linear function
   * x^(t+1) = a (t+1) + b. The parameters a and b are computed from the known data.
   *
   * For quadratic extrapolation three equally spaced samples x^t, x^(t-1), x^(t-2) are used to compute x^(t+1) using a quadratic function
   * x^(t+1) = a (t+1)^2 + b (t+1) + c. The parameters a, b and c are computed from the known data.
   * 
   * @param order Order of the extrapolation scheme to be used.
   * @param timeWindows number of valid samples.
   */
  Eigen::VectorXd extrapolateData(int order, int timeWindows);
};

} // namespace time
} // namespace precice

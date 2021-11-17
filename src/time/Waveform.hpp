#pragma once

#include <Eigen/Core>
#include "logging/Logger.hpp"

namespace precice {

namespace testing {
// Forward declaration to friend the boost test struct
class WaveformFixture;
} // namespace testing

namespace time {

class Waveform {
  friend class testing::WaveformFixture; // Make the fixture friend of this class
public:
  /**
   * @brief Waveform object which stores values of current and past time windows for performing extrapolation.
   * @param valuesSize defines how many values one sample in time consists of
   * @param extrapolatioOrder defines the maximum extrapolation order supported by this Waveform and reserves storage correspondingly
   */
  Waveform(const int valuesSize,
           const int extrapolationOrder);

  /**
   * @brief Updates entry in _timeWindows corresponding to this window with given values
   * @param values new sample for this time window
   */
  void store(const Eigen::VectorXd &values);

  /**
   * @brief Called, when moving to the next time window. All entries in _timeWindows are shifted. The new entry is initialized as the value from the last window (= constant extrapolation)
   */
  void moveToNextWindow();

  /**
   * @brief getter for values at the current time window.
   */
  const Eigen::VectorXd getInitialGuess();

private:
  /// Stores values for several time windows.
  Eigen::MatrixXd _timeWindowsStorage;

  /// extrapolation order for this waveform
  const int _extrapolationOrder;

  /// number of stored samples in _timeWindowsStorage
  int _numberOfStoredSamples;

  /**
   * @brief returns number samples in time this waveform can store
   */
  int sizeOfSampleStorage();

  /**
   * @brief returns number of values per sample in time stored by this waveform
   */
  int valuesSize();

  mutable logging::Logger _log{"time::Waveform"};

  /**
   * @brief Extrapolates values from _timeWindowsStorage using an extrapolation scheme of given order. 
   * 
   * If the order condition cannot be satisfied, since there are not enough samples available, the order is automatically reduced.
   * If order two is required, but only two samples are available, the extrapolation order is automatically reduced to one.
   */
  Eigen::VectorXd extrapolate();
};

} // namespace time
} // namespace precice

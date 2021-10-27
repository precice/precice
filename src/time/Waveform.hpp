#pragma once

#include <Eigen/Core>
#include "logging/Logger.hpp"

namespace precice {

namespace testing {
// Forward declaration to friend the boost test struct
struct WaveformFixture;
} // namespace testing

namespace time {

class Waveform {
  friend struct testing::WaveformFixture; // Make the fixture friend of this class
public:
  /**
   * @brief Waveform object which stores data of current and past time windows for performing extrapolation.
   * @param dataCount defines how many pieces of data one sample in time consists of
   * @param extrapolatioOrder defines the maximum extrapolation order supported by this Waveform and reserves storage correspondingly
   */
  Waveform(const int dataCount,
           const int extrapolationOrder);

  /**
   * @brief Updates entry in _timeWindows corresponding to this window with given data
   * @param data new sample for this time window
   */
  void store(const Eigen::VectorXd &data);

  /**
   * @brief Called, when moving to the next time window. All entries in _timeWindows are shifted. The new entry is initialized as the value from the last window (= constant extrapolation)
   */
  void moveToNextWindow();

  /**
   * @brief getter for data at the current time window.
   */
  const Eigen::VectorXd getInitialGuess();

private:
  /// Data values of time windows.
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
   * @brief returns number of data per sample in time stored by this waveform
   */
  int dataCount(); // @todo bad naming, consider renaming. See https://github.com/precice/precice/pull/1094#pullrequestreview-771715472

  mutable logging::Logger _log{"time::Waveform"};

  /**
   * @brief Extrapolates data _timeWindows using an extrapolation scheme of given order. 
   * 
   * If the order condition cannot be satisfied, since there are not enough samples available, the order is automatically reduced.
   * If order two is required, but only two samples are available, the extrapolation order is automatically reduced to one.
   */
  Eigen::VectorXd extrapolateData();
};

} // namespace time
} // namespace precice

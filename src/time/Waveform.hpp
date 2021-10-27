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
  /// To be used, when the interpolation order is not defined for this Waveform.
  static const int UNDEFINED_INTERPOLATION_ORDER;

  /**
   * @brief Waveform object which stores data of current and past time windows for performing extrapolation.
   * @param dataCount defines how many pieces of data one sample in time consists of
   * @param extrapolatioOrder defines the maximum extrapolation order supported by this Waveform and reserves storage correspondingly
   * @param interpolationOrder defines the maximum interpolation order supported by this Waveform and reserves storage correspondingly
   */
  Waveform(const int dataCount,
           const int extrapolationOrder,
           const int interpolationOrder);

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
   * @brief sample Waveform. Uses interpolation with Waveform's interpolation order, if necessary
   * @param normalizedDt time where the sampling inside the window happens. 0 refers to the beginning of the window and 1 to the end.
   */
  Eigen::VectorXd sample(const double normalizedDt);

  /**
   * @brief getter for data at the current time window.
   */
  const Eigen::VectorXd getInitialGuess();

private:
  /// Data values of time windows.
  Eigen::MatrixXd _timeWindowsStorage;

  /// extrapolation order for this waveform
  const int _extrapolationOrder;

  /// interpolation order for this waveform
  const int _interpolationOrder;

  /// number of stored samples in _timeWindowsStorage
  int _numberOfStoredSamples;

  /**
   * @brief returns number samples in time this waveform can store
   */
  int sizeOfSampleStorage();

  /**
   * @brief returns number of data per sample in time stored by this waveform
   */
  int dataSize();

  mutable logging::Logger _log{"time::Waveform"};

  /**
   * @brief Extrapolates data _timeWindows using an extrapolation scheme of given order. 
   * 
   * If the order condition cannot be satisfied, since there are not enough samples available, the order is automatically reduced.
   * If order two is required, but only two samples are available, the extrapolation order is automatically reduced to one.
   */
  Eigen::VectorXd extrapolateData();

  /**
   * @brief Interpolates data inside current time window using an interpolation scheme of the order of this Waveform.
   *
   * @param normalizedDt time where the sampling inside the window happens. 0 refers to the beginning of the window and 1 to the end.
   */
  Eigen::VectorXd interpolateData(const double normalizedDt);
};

} // namespace time
} // namespace precice

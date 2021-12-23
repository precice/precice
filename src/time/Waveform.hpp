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
   * @brief Waveform object which stores values of current and past time windows for performing interpolation. 
   *
   * Storage still needs to be initialized with Waveform::initialize, before the Waveform can be used.
   *
   * @param interpolationOrder defines the maximum interpolation order supported by this Waveform and reserves storage correspondingly
   */
  Waveform(const int interpolationOrder);

  /**
   * @brief Used to initialize _timeWindowsStorage according to required size.
   * @param valuesSize defines how many values one sample in time consists of
   */
  void initialize(const int valuesSize);

  /**
   * @brief Updates entry in _timeWindows corresponding to this window with given values
   * @param values new sample for this time window
   */
  void store(const Eigen::VectorXd &values);

  /**
   * @brief Updates entry in _timeWindowsStorage corresponding to a given column ID with given values
   * @param values new sample for this time window
   * @param columnID ID of column to be updated
   */
  void storeAt(const Eigen::VectorXd values, int columnID);

  /**
   * @brief Called, when moving to the next time window. All entries in _timeWindows are shifted. The new entry is initialized as the value from the last window (= constant extrapolation)
   */
  void moveToNextWindow();

  /**
   * @brief sample Waveform. Uses interpolation with Waveform's interpolation order, if necessary
   * @param normalizedDt time where the sampling inside the window happens. 0 refers to the beginning of the window and 1 to the end.
   */
  Eigen::VectorXd sample(const double normalizedDt);

  // @todo try to make private!
  /**
   * @brief returns number samples in time this waveform can store
   */
  int sizeOfSampleStorage();

  // @todo try to make private!
  /**
   * @brief getter for a sample in _timeWindowsStorage
   */
  Eigen::VectorXd getSample(int sampleID);

private:
  /// Set by initialize. Used for consistency checks.
  bool _storageIsInitialized = false;

  /// Stores values for several time windows.
  Eigen::MatrixXd _timeWindowsStorage;

  /// interpolation order for this waveform
  const int _interpolationOrder;

  /// number of stored samples in _timeWindowsStorage
  int _numberOfStoredSamples;

  mutable logging::Logger _log{"time::Waveform"};

  /**
   * @brief returns number of values per sample in time stored by this waveform
   */
  int valuesSize();

  /**
   * @brief Interpolates values inside current time window using _timeWindowsStorage and an interpolation scheme of the order of this Waveform.
   *
   * @param normalizedDt time where the sampling inside the window happens. 0 refers to the beginning of the window and 1 to the end.
   */
  Eigen::VectorXd interpolate(const double normalizedDt);
};

} // namespace time
} // namespace precice

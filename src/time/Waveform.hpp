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
   * @brief Updates first entry in _timeWindows with given values
   * @param values new sample for this time window
   */
  void storeAtFirstSample(const Eigen::VectorXd &values);

  /**
   * @brief Updates all entries in _timeWindowsStorage with given values
   * @param values new sample for this time window
   */
  void storeAtAllSamples(const Eigen::VectorXd values);

  /**
   * @brief getter for first entry in _timeWindowsStorage
   */
  Eigen::VectorXd readAtFirstSample();

  /**
   * @brief Called, when moving to the next time window. All entries in _timeWindows are shifted. The new entry is initialized as the value from the last window (= constant extrapolation)
   */
  void moveToNextWindow();

  /**
   * @brief sample Waveform. Uses interpolation with Waveform's interpolation order, if necessary.
   * @param normalizedDt time where the sampling inside the window happens. Only allows values between 0 and 1. 0 refers to the beginning of the window and 1 to the end.
   */
  Eigen::VectorXd sample(const double normalizedDt);

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
   * @brief returns maximum number of samples in time this waveform can store
   */
  int maxNumberOfStoredSamples();

  /**
   * @brief Updates entry in _timeWindowsStorage corresponding to a given sampleIndex with given values
   * @param values new sample for this time window
   * @param sampleIndex index of sample to be updated
   */
  void storeAt(const Eigen::VectorXd values, int sampleIndex);

  /**
   * @brief Interpolates values inside current time window using _timeWindowsStorage and an interpolation scheme of the order of this Waveform.
   *
   * @param normalizedDt time where the sampling inside the window happens. 0 refers to the beginning of the window and 1 to the end.
   */
  Eigen::VectorXd interpolate(const double normalizedDt);
};

} // namespace time
} // namespace precice

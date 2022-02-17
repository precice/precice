#pragma once

#include <Eigen/Core>
#include "logging/Logger.hpp"

namespace precice {

namespace testing {
// Forward declaration to friend the boost test struct
class WaveformFixture;
} // namespace testing

namespace time {

/**
 * @brief Stores data samples in time and allows to perform interpolation to obtain new values.
 */
class Waveform {
  friend class testing::WaveformFixture; // Make the fixture friend of this class
public:
  /**
   * @brief Waveform object which stores values of current and past time windows for performing interpolation. 
   *
   * Storage still needs to be initialized with Waveform::initialize, before the Waveform can be used.
   *
   * @param interpolationOrder Defines the maximum interpolation order supported by this Waveform and reserves storage correspondingly
   */
  Waveform(const int interpolationOrder);

  /**
   * @brief Get the _interpolationOrder.
   * 
   * @return int _interpolationOrder
   */
  int getInterpolationOrder() const;

  /**
   * @brief Used to initialize _timeWindowsStorage according to required size.
   * @param valuesSize Defines how many values one sample in time consists of
   */
  void initialize(const int valuesSize);

  /**
   * @brief Updates first entry in _timeWindows with given values.
   * @param values Sample for this time window
   */
  void storeAtFirstSample(const Eigen::VectorXd &values);

  /**
   * @brief Updates all entries in _timeWindowsStorage with given values. Often used to initialize a new Waveform with given data.
   * @param values Sample used for all time window
   */
  void storeAtAllSamples(const Eigen::VectorXd values);

  /**
   * @brief Shifts all entries in _timeWindows. The new entry is initialized as the value from the last window (= constant extrapolation). Usually called, when moving to the next time window.
   */
  void moveToNextWindow();

  /**
   * @brief Evaluate Waveform at specific point in time. Uses interpolation with Waveform's interpolation order, if necessary.
   * @param normalizedDt Time where the sampling inside the window happens. Only allows values between 0 and 1. 0 refers to the beginning of the window and 1 to the end.
   * @return Value of Waveform at time normalizedDt.
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
   * @brief Get number of values per sample in time stored by this waveform.
   * @return Number of values per sample.
   */
  int valuesSize();

  /**
   * @brief Get maximum number of samples in time this waveform can store.
   * @return Maximum number of samples.
   */
  int maxNumberOfStoredSamples();

  /**
   * @brief Updates entry in _timeWindowsStorage corresponding to a given sampleIndex with given values.
   * @param values Input sample.
   * @param sampleIndex Index of sample to be updated to given input sample.
   */
  void storeAt(const Eigen::VectorXd values, int sampleIndex);

  /**
   * @brief Computes which order may be used for interpolation.
   * 
   * Order of interpolation is determined by number of stored samples and maximum order defined by the user.
   * Example: If only two samples are available, the maximum order we may use is 1, even if the user demands order 2.
   *
   * @param requestedOrder Order requested by the user.
   * @param numberOfAvailableSamples Samples available for interpolation.
   * @return Order that may be used.
   */
  int computeUsedOrder(int requestedOrder, int numberOfAvailableSamples);

  /**
   * @brief Interpolates values inside current time window using _timeWindowsStorage and an interpolation scheme of the order of this Waveform.
   * @param normalizedDt time where the sampling inside the window happens. 0 refers to the beginning of the window and 1 to the end.
   * @return Interpolated value at time normalizedDt.
   */
  Eigen::VectorXd interpolate(const double normalizedDt);
};

} // namespace time
} // namespace precice

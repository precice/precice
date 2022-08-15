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
 *
 * When created via Waveform(interpolationOrder) a waveform reserves storage for the samples that are used to create the interpolant.
 * After creation of the waveform it must be initialized with Waveform::initialize(value) to finally reserve the storage.
 * The waveform is initialized with one data value as a constant function.
 * Waveform::store(value) allows the user to update the data sample in the Waveform.
 * Each time the user calls Waveform::moveToNextWindow() the data provided through Waveform::initialize or Waveform::store will be locked.
 * With each call of Waveform::moveToNextWindow() the number of available samples and, therefore, the order of the waveform is increased by one until the interpolation order that was defined during construction is reached.
 * As soon as this interpolation order is reached, the oldest sample will be discarded when Waveform::moveToNextWindow() is called. The order will then stay the same.
 */
class Waveform {
  friend class testing::WaveformFixture; // Make the fixture friend of this class
public:
  /**
   * @brief Waveform object which stores values of current and past time windows for performing interpolation.
   *
   * Storage still needs to be initialized with Waveform::initialize, before the Waveform can be used.
   *
   * @param interpolationOrder Defines the interpolation order supported by this Waveform and reserves storage correspondingly
   */
  Waveform(const int interpolationOrder);

  /**
   * @brief Get the _interpolationOrder.
   *
   * @return int _interpolationOrder
   */
  int getInterpolationOrder() const;

  /**
   * @brief Used to initialize _timeStepsStorage according to required size and initializes Waveform as constant with given values.
   * @param values Defines constant initial value of waveform and its size
   */
  void initialize(const Eigen::VectorXd &values);

  /**
   * @brief Updates first entry in _timeWindows with given value.
   * @param values Sample at end of this time window
   */
  void store(const Eigen::VectorXd &values);

  /**
   * @brief Updates an entry for dt in _timeWindows with given value.
   * @param values Sample at dt in this time window
   * @param normalizedDt normalizedDt associated with this value. Only allows values between 0 and 1. 0 refers to the beginning of the window and 1 to the end.
   */
  void store(const Eigen::VectorXd &values, double normalizedDt);

  /**
   * @brief Shifts all entries in _timeWindows. The new entry is initialized as the value from the last window (= constant extrapolation). Called when moving to the next time window.
   */
  void moveToNextWindow();

  /**
   * @brief Evaluate waveform at specific point in time. Uses interpolation if necessary.
   *
   * Interpolates values inside current time window using _timeStepsStorage and an interpolation scheme of the order of this Waveform.
   *
   * @param normalizedDt Time where the sampling inside the window happens. Only allows values between 0 and 1. 0 refers to the beginning of the window and 1 to the end.
   * @return Value of Waveform at time normalizedDt.
   */
  Eigen::VectorXd sample(const double normalizedDt);

private:
  /// Set by initialize. Used for consistency checks.
  bool _storageIsInitialized = false;

  /// Stores values on the current window.
  std::map<double, Eigen::VectorXd> _timeStepsStorage;

  /// interpolation order for this waveform
  const int _interpolationOrder;

  /// number of stored samples in _timeStepsStorage
  int _numberOfStoredSamples;

  mutable logging::Logger _log{"time::Waveform"};

  /**
   * @brief Get number of values per sample in time stored by this waveform.
   * @return Number of values per sample.
   */
  int valuesSize();

  /**
   * @brief Removes all data stored in _timeStepsStorage
   *
   * Used to clear time steps storage before a new window is entered or a window is repeated
   *
   */
  void clearTimeStepsStorage();

  /**
   * @brief Get maximum number of samples in time this waveform can store.
   * @return Maximum number of samples.
   */
  int maxNumberOfStoredWindows();

  /**
   * @brief Get maximum dt that is stored in this waveform.
   *
   * Used to check whether a user is trying to add a sample associated with a dt that is smaller than the maximum dt. This is forbidden, because the waveform is locked for times that are smaller than the maximum dt.
   *
   * @return the maximum dt from _timeStepsStorage
   */
  double maxStoredDt();

  /**
   * @brief Updates entry in _timeStepsStorage corresponding to a given sampleIndex with given values.
   * @param values Input sample.
   * @param dt dt associated with this value. Only allows values between 0 and 1. 0 refers to the beginning of the window and 1 to the end.
   */
  void storeAt(const Eigen::VectorXd values, double dt);

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
   * @brief Returns point the closest time stored in _timeStepsStorage that is before normalizedDt
   *
   * @param normalizedDt point in time
   * @return double point in time before normalizedDt in _timeStepsStorage
   */
  double findTimeBefore(double normalizedDt);

  /**
   * @brief Returns point the closest time stored in _timeStepsStorage that is after normalizedDt
   *
   * @param normalizedDt point in time
   * @return double point in time after normalizedDt in _timeStepsStorage
   */
  double findTimeAfter(double normalizedDt);
};

} // namespace time
} // namespace precice

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
   * @brief Updates an entry for dt in _timeWindows with given value.
   * @param values Sample at dt in this time window
   * @param normalizedDt normalizedDt associated with this value. Only allows values between 0 and 1. 0 refers to the beginning of the window and 1 to the end.
   */
  void store(const Eigen::VectorXd &values, double normalizedDt = 1.0);

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
  /** @TODO Idea for more efficient data structure and redesign (do this when functionality is working and tested!)
   *   1. use Eigen::MatrixXd instead of map for _timeStepsStorage.
   *   2. create a member std::map<double, int> _timeSteps where (unique) time is mapped to column index of _timeStepsStorage that holds the corresponding sample. (Alternative: Use another Eigen::VectorXd to store times, but this enforces maintaining a consistent order for _timeSteps and _timeStepsStorage. This sounds complicated.)
   */
  /// Stores values on the current window.
  std::map<double, Eigen::VectorXd> _timeStepsStorage;

  /// interpolation order for this waveform
  const int _interpolationOrder;

  mutable logging::Logger _log{"time::Waveform"};

  /**
   * @brief Get maximum dt that is stored in this waveform.
   *
   * Used to check whether a user is trying to add a sample associated with a dt that is smaller than the maximum dt. This is forbidden, because the waveform is locked for times that are smaller than the maximum dt.
   *
   * @return the maximum dt from _timeStepsStorage
   */
  double maxStoredDt();

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
   * @brief Returns point the closest time stored in _timeStepsStorage that is after normalizedDt
   *
   * @param normalizedDt point in time
   * @return double point in time after normalizedDt in _timeStepsStorage
   */
  double findTimeAfter(double normalizedDt);

  /**
   * @brief Get keys of _timeStepsStorage in ascending order. Starting from low to high.
   *
   * @return Eigen::VectorXd
   */
  Eigen::VectorXd getTimesAscending();
};

} // namespace time
} // namespace precice

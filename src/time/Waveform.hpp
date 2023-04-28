#pragma once

#include <Eigen/Core>
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"
#include "time/Storage.hpp"

namespace precice {

namespace testing {
// Forward declaration to friend the boost test struct
class WaveformFixture;
} // namespace testing

namespace time {

/**
 * @brief Stores data samples in time and allows to perform interpolation.
 *
 * The constructor Waveform(interpolationOrder) creates a Storage for the samples that are used to create the interpolant.
 * After creation of the waveform it must be initialized with Waveform::initialize(value) to finally reserve the storage. All newly added values must have the same dimension as the initially provided values.
 * The waveform is initialized with two data values at the beginning and at the end of the window as a constant function. Waveform::store(value) allows the user to provide new data to the Waveform. Interpolation is performed based on these values.
 * The available interpolation order depends on the number of stored samples and can reach the interpolationOrder defined during construction as a maximum. If more samples are available than the maximum order requires, a piecewise interpolation will be used (piecewise constant, piecewise linear and B-Spline interpolation).
 * Interpolation is only performed inside the current time window. If the waveform should be used for the next time window, Waveform::moveToNextWindow() allows the user to use the data at the end of the current window as an initial guess for the next window.
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
   * @param data pointer to data this waveform interpolates
   */
  Waveform(const int interpolationOrder, mesh::PtrData data);

  /**
   * @brief Get the _interpolationOrder.
   *
   * @return int _interpolationOrder
   */
  int getInterpolationOrder() const;

  /**
   * @brief Updates an entry for normalizedDt in _timeWindows with given value.
   * @param values Sample at normalizedDt in this time window
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
   * Interpolates values inside current time window using _storage and an interpolation scheme of the order of this Waveform. The interpolation scheme always uses all available values in _storage and tries to reach _interpolationOrder. If more than the required number of values needed to reach _interpolationOrder are available, a piecewise interpolation strategy will be applied to obtain an interpolation that reaches the requested order and still interpolates all the provided data points.
   *
   * @param normalizedDt Time where the sampling inside the window happens. Only allows values between 0 and 1. 0 refers to the beginning of the window and 1 to the end.
   * @return Value of Waveform at time normalizedDt.
   */
  Eigen::VectorXd sample(const double normalizedDt) const;

private:
  /// Stores values.
  mesh::PtrData _data;

  /// interpolation order for this waveform
  const int _interpolationOrder;

  mutable logging::Logger _log{"time::Waveform"};

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
  int computeUsedOrder(int requestedOrder, int numberOfAvailableSamples) const;
};

} // namespace time
} // namespace precice

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
// @todo Refactor Waveform class. Move sample function inside of Storage::sample
/**
 * @brief Allows to perform interpolation on samples in storage of given data.
 *
 * The constructor Waveform(degree, data) creates a waveform. The samples of the data's storage are used to create the interpolant.
 * The waveform is initialized with two data values at the beginning and at the end of the window as a constant function. Waveform::store(value) allows the user to provide new data to the Waveform. Interpolation is performed based on these values.
 * The maximum allowed polynomial degree depends on the number of stored samples and can reach the interpolationDegree defined during construction as a maximum. If more samples are available than the maximum degree requires, a piecewise interpolation will be used (piecewise constant, piecewise linear and B-Spline interpolation).
 * Interpolation is only performed inside the current time window.
 */
class Waveform {
  friend class testing::WaveformFixture; // Make the fixture friend of this class
public:
  /**
   * @brief Waveform object which stores values of current and past time windows for performing interpolation.
   *
   * @param degree Defines the polynomial degree supported by this Waveform and reserves storage correspondingly
   */
  Waveform(const int degree);

  /**
   * @brief Get the _degree.
   *
   * @return int _degree
   */
  int getDegree() const;

  /// Returns a reference to the _timeStepsStorage.
  time::Storage &timeStepsStorage();

  /// Returns a const reference to the _timeStepsStorage.
  const time::Storage &timeStepsStorage() const;

  /// Returns a the stamples from _timeStepsStorage.
  auto stamples() const
  {
    return _timeStepsStorage.stamples();
  }

  /**
   * @brief Evaluate waveform at specific point in time. Uses interpolation if necessary.
   *
   * Interpolates values inside current time window using _storage and an interpolation scheme of the maximum degree of this Waveform. The interpolation scheme always uses all available values in _storage and tries to reach _degree. If more than the required number of values needed to reach _degree are available, a piecewise interpolation strategy will be applied to obtain an interpolation that reaches the requested polynomial degree and still interpolates all the provided data points.
   *
   * @param time Time where the sampling inside the window happens.
   * @return Value of Waveform at given time.
   */
  Eigen::VectorXd sample(const double time) const;

private:
  /// Stores time steps in the current time window
  time::Storage _timeStepsStorage;

  /// interpolation degree for this waveform
  int _degree;

  mutable logging::Logger _log{"time::Waveform"};

  /**
   * @brief Computes which degree may be used for interpolation.
   *
   * Actual degree of interpolating B-spline is determined by number of stored samples and maximum degree defined by the user.
   * Example: If only two samples are available, the maximum degree we may use is 1, even if the user demands degree 2.
   *
   * @param requestedDegree B-spline degree requested by the user.
   * @param numberOfAvailableSamples Samples available for interpolation.
   * @return B-spline degree that may be used.
   */
  int computeUsedDegree(int requestedDegree, int numberOfAvailableSamples) const;
};

} // namespace time
} // namespace precice

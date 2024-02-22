#pragma once

#include <optional>

#include "math/differences.hpp"
#include "utils/DoubleAggregator.hpp"
#include "utils/assertion.hpp"

namespace precice::cplscheme::impl {

/** Handler for stepping forward in time using time windows
 *
 * The handler respects a maximum time if it was set in the constructor.
 *
 * The goal of this class is to handle all computation around time and time steps
 * with maximizing accuracy and simplifying calling code.
 */
class TimeHandler {
public:
  TimeHandler() = default;

  /// Constructor for setting an optional maximum time
  TimeHandler(std::optional<double> maxTime)
      : _maxTime(maxTime)
  {
  }

  TimeHandler(const TimeHandler &) = default;
  TimeHandler(TimeHandler &&)      = default;

  TimeHandler &operator=(const TimeHandler &) = default;
  TimeHandler &operator=(TimeHandler &&) = default;

  /// @name State altering functions
  /// @{

  /// Resets the handler to the given time
  void resetTo(double timeStart)
  {
    _windowStart = timeStart;
    resetProgress();
  }

  /// Progress the time window by the given amount
  void progressBy(double dt)
  {
    _windowProgress.add(dt);
  }

  /// Resets the progress of the time window back to 0
  void resetProgress()
  {
    _windowProgress = 0.0;
  }

  /** Complete the time window moving the start of the time window
   *
   * The given exact time window size will be used if it is approximately equal to
   * the accumulated progress.
   */
  void completeTimeWindow(double timeWindowSize)
  {
    // Use the exact time window if possible
    if (math::equals(windowProgress(), timeWindowSize)) {
      _windowStart.add(timeWindowSize);
    } else {
      // Handle truncated time windows in case of reaching the end of the simulation
      // This only happens when the final time window is truncated due
      // time window size not being a divider of max-time.
      PRECICE_ASSERT(reachedEnd(), "This can only happened if we reach max-time", timeWindowSize, time(), _maxTime.value_or(-1));
      _windowStart.add(_windowProgress);
    }
    resetProgress();
  }

  /// @}
  /// @name Queries about reaching ends
  /// @{

  /** Has the end of the time window been reached?
   * This respects max time.
   */
  bool reachedEndOfWindow(double timeWindowSize) const
  {
    return math::equals(untilWindowEnd(timeWindowSize), 0.0);
  }

  /** Has the end of the overall time been reached?
   * This is always false if no max time has been defined.
   */
  bool reachedEnd() const
  {
    if (!_maxTime) {
      return false;
    }
    return math::equals(untilTime(*_maxTime), 0.0);
  }

  /// @}
  /// @name Time differences
  /// @{

  /// Returns the time distance to the possibly truncated end of the current time window
  double untilWindowEnd(double timeWindowSize) const
  {
    double maxDt = untilProgress(timeWindowSize);
    if (!_maxTime) {
      return maxDt;
    }
    // Truncate by max Time
    return std::min(maxDt, untilEnd());
  }

  /// Returns the time difference until the overall time reaches the given time
  double untilTime(double t) const
  {
    // rearranged version of: maxtime - windowStart - windowProgress
    return -(_windowStart + _windowProgress.value() - t).value();
  }

  /** Returns the time difference until the end of the overall time.
    * This returns infinity if there is no maximum time defined.
    */
  double untilEnd() const
  {
    if (_maxTime) {
      return untilTime(*_maxTime);
    }
    return std::numeric_limits<double>::max();
  }

  /// Returns the window progress as a double
  double windowProgress() const
  {
    return _windowProgress.value();
  }

  /// @}
  /// @name Time points
  /// @{

  /// Returns the window start as a double
  double windowStart() const
  {
    return _windowStart.value();
  }

  /// Returns the current time as a double
  double time() const
  {
    return (_windowStart + _windowProgress.value()).value();
  }

  /// @}

private:
  /// The optional maximum time
  std::optional<double> _maxTime = std::nullopt;
  /// The aggregator for the start of a time window, handles timestep.
  utils::DoubleAggregator _windowStart;
  /// The aggregator progress inside a time window, handles substeps.
  utils::DoubleAggregator _windowProgress;

  /// Returns the time difference until the progress reaches the given window size
  double untilProgress(double windowSize) const
  {
    // rearranged version of: windowSize - _windowProgress
    return -(_windowProgress - windowSize).value();
  }
};

} // namespace precice::cplscheme::impl

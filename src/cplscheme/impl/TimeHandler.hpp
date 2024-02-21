#pragma once

#include "utils/DoubleAggregator.hpp"

namespace precice::cplscheme::impl {

/// Handler for stepping forward in time using time windows
class TimeHandler {
public:
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

  /// Complete the time window, but use the exact progress instead of the accumulated one
  void completeTimeWindowExact(double exactProgress)
  {
    // @TODO prevent major differences
    _windowStart.add(exactProgress);
    resetProgress();
  }

  /// Complete the time window using the accumulated progress
  void completeTimeWindow()
  {
    _windowStart.add(_windowProgress);
    resetProgress();
  }

  /// Resets the progress of the time window back to 0
  void resetProgress()
  {
    _windowProgress = 0.0;
  }

  /// Returns the time difference until the progress reaches the given window size
  double untilProgress(double windowSize) const
  {
    // rearranged version of: windowSize - _windowProgress
    return -(_windowProgress - windowSize).value();
  }

  /// Returns the time difference until the overall time reaches the given time
  double untilTime(double maxTime) const
  {
    // rearranged version of: maxtime - windowStart - windowProgress
    return -(_windowStart + _windowProgress.value() - maxTime).value();
  }

  /// Returns the window progress as a double
  double windowProgress() const
  {
    return _windowProgress.value();
  }

  /// Returns the window start as a double
  double windowStart() const
  {
    return _windowStart.value();
  }

  /// Returns the current time as a double
  double time() const
  {
    return (_windowStart + _windowProgress).value();
  }

private:
  utils::DoubleAggregator _windowStart;
  utils::DoubleAggregator _windowProgress;
};

} // namespace precice::cplscheme::impl

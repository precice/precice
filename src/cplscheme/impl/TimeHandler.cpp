#include "TimeHandler.hpp"

#include <optional>

#include "math/differences.hpp"
#include "utils/DoubleAggregator.hpp"
#include "utils/assertion.hpp"

namespace precice::cplscheme::impl {

struct TimeHandler::Impl {

  /// The optional maximum time
  std::optional<double> _maxTime = std::nullopt;

  /// The aggregator for the start of a time window, handles timestep.
  utils::DoubleAggregator _windowStart;

  /// The aggregator progress inside a time window, handles substeps.
  utils::DoubleAggregator _windowProgress;
};

// TimeHandler implementation

TimeHandler::TimeHandler()
    : _impl(std::make_unique<TimeHandler::Impl>())
{
}

// Required by unique_ptr
TimeHandler::~TimeHandler() = default;

TimeHandler::TimeHandler(std::optional<double> maxTime)
    : _impl(std::make_unique<TimeHandler::Impl>())
{
  _impl->_maxTime = maxTime;
}

TimeHandler::TimeHandler(const TimeHandler &other)
    : _impl(std::make_unique<TimeHandler::Impl>(*other._impl))
{
}

TimeHandler::TimeHandler(TimeHandler &&other)
    : _impl(std::exchange(other._impl, nullptr))
{
}

TimeHandler &TimeHandler::operator=(const TimeHandler &other)
{
  _impl = std::make_unique<TimeHandler::Impl>(*other._impl);
  return *this;
}

TimeHandler &TimeHandler::operator=(TimeHandler &&other)
{
  _impl = std::exchange(other._impl, nullptr);
  return *this;
}

void TimeHandler::resetTo(double timeStart)
{
  _impl->_windowStart = timeStart;
  resetProgress();
}

void TimeHandler::progressBy(double dt)
{
  _impl->_windowProgress.add(dt);
}

void TimeHandler::resetProgress()
{
  _impl->_windowProgress = 0.0;
}

void TimeHandler::completeTimeWindow(double timeWindowSize)
{
  // Use the exact time window if possible
  if (math::equals(windowProgress(), timeWindowSize)) {
    _impl->_windowStart.add(timeWindowSize);
  } else {
    // Handle truncated time windows in case of reaching the end of the simulation
    // This only happens when the final time window is truncated due
    // time window size not being a divider of max-time.
    PRECICE_ASSERT(reachedEnd(), "This can only happened if we reach max-time", timeWindowSize, time(), _impl->_maxTime.value_or(-1));
    _impl->_windowStart.add(_impl->_windowProgress);
  }
  resetProgress();
}

// Queries about reaching ends

bool TimeHandler::reachedEndOfWindow(double timeWindowSize) const
{
  return math::equals(untilWindowEnd(timeWindowSize), 0.0);
}

bool TimeHandler::reachedEnd() const
{
  if (!_impl->_maxTime) {
    return false; // Directly return preventing lossy computations
  }
  return math::equals(untilTime(*_impl->_maxTime), 0.0);
}

// Time differences

double TimeHandler::untilWindowEnd(double timeWindowSize) const
{
  double maxDt = untilProgress(timeWindowSize);
  if (!_impl->_maxTime) {
    return maxDt; // Return to prevent further lossy computations
  }
  // Truncate by max Time
  return std::min(maxDt, untilEnd());
}

double TimeHandler::untilTime(double t) const
{
  // rearranged version of: maxtime - windowStart - windowProgress
  return -(_impl->_windowStart + _impl->_windowProgress.value() - t).value();
}

double TimeHandler::untilEnd() const
{
  if (!_impl->_maxTime) {
    return std::numeric_limits<double>::max();
  }
  return untilTime(_impl->_maxTime.value());
}

double TimeHandler::windowProgress() const
{
  return _impl->_windowProgress.value();
}

// Time points

double TimeHandler::windowStart() const
{
  return _impl->_windowStart.value();
}

double TimeHandler::time() const
{
  return (_impl->_windowStart + _impl->_windowProgress.value()).value();
}

double TimeHandler::untilProgress(double windowSize) const
{
  // rearranged version of: windowSize - _windowProgress
  return -(_impl->_windowProgress - windowSize).value();
}

} // namespace precice::cplscheme::impl

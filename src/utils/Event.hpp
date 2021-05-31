#pragma once

#include <chrono>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include "logging/Logger.hpp"

namespace precice {
namespace utils {

/// Represents an event that can be started and stopped.
/** Additionally to the duration there is a special property that can be set for a event.
A property is a a key-value pair with a numerical value that can be used to trace certain events,
like MPI calls in an event. It is intended to be set by the user. */
class Event {
public:
  enum class State : int {
    STOPPED = 0,
    STARTED = 1,
    PAUSED  = 2,
  };

  /// Default clock type. All other chrono types are derived from it.
  using Clock = std::chrono::steady_clock;

  using StateChanges = std::vector<std::pair<State, Clock::time_point>>;

  using Data = std::map<std::string, std::vector<int>>;

  /// An Event can't be copied.
  Event(const Event &other) = delete;

  /// Name used to identify the timer. Events of the same name are accumulated to
  std::string name;

  /// Allows to put a non-measured (i.e. with a given duration) Event to the measurements.
  Event(std::string eventName, Clock::duration initialDuration);

  /// Creates a new event and starts it, unless autostart = false, synchronize processes, when barrier == true
  /** Use barrier == true with caution, as it can lead to deadlocks. */
  Event(std::string eventName, bool barrier = false, bool autostart = true);

  /// Stops the event if it's running and report its times to the EventRegistry
  ~Event();

  /// Starts an event. If it's already started it has no effect.
  void start(bool barrier = false);

  /// Stops an event and commit it. If it's already stopped it has no effect.
  void stop(bool barrier = false);

  /// Pauses an event, does not commit. If it's already paused it has no effect.
  void pause(bool barrier = false);

  /// Gets the duration of the event.
  Clock::duration getDuration() const;

  /// Adds named integer data, associated to an event.
  void addData(std::string key, int value);

  Data data;

  StateChanges stateChanges;

private:
  logging::Logger _log{"utils::Events"};

  Clock::time_point starttime;
  Clock::duration   duration = Clock::duration::zero();
  State             state    = State::STOPPED;
  bool              _barrier = false;
};

/// Class that changes the prefix in its scope
class ScopedEventPrefix {
public:
  ScopedEventPrefix(const std::string &name);

  ~ScopedEventPrefix();

private:
  std::string previousName = "";
};

} // namespace utils
} // namespace precice

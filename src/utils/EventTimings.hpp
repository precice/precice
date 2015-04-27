#pragma once

#include <chrono>
#include <map>
#include <string>

namespace precice {
namespace utils {

/// Represents an event that can be started and stopped.
/** Additionally to the duration there is a special property that can be set for a event.
A property is a a key-value pair with a numerical value that can be used to trace certain events,
like MPI calls in an event. It intended to be set by the user. */
class Event
{
public:
  /// Default high precision clock type. All other chrono types are derived from it.
  using Clock = std::chrono::high_resolution_clock;

  /// Name used to identify the timer. Events of the same name are accumulated to
  std::string name;

  Event(std::string eventName, Clock::duration eventDuration);

  /// Creates a new event and starts it, unless autostart = false
  Event(std::string eventName, bool barrier = false, bool autostart = true);

  /// Stops the event if it's running and report its times to the EventRegistry
  ~Event();

  /// Starts an event. If it's already started it has no effect.
  void start(bool barrier = false);

  /// Stops an event. If it's already stopped it has no effect.
  void stop(bool barrier = false);

  /// Gets the duration of the event.
  Clock::duration getDuration();

  /// Map of additional properties that can be set by the user.
  std::map<std::string, double> properties;

  /// Adds the value to the propety.
  void addProp(std::string property, double value);

private:
  Clock::time_point starttime;
  Clock::time_point stoptime;
  Clock::duration duration = Clock::duration::zero();
  bool isStarted = false;
  bool _barrier = false;
};


/// Class that aggregates data (durations and properties) for a specific event.
class EventData
{
public:
  /// Adds an events data.
  void put(Event* event);

  /// Get the average duration of all events so far.
  int getAvg();

  /// Get the maximum duration of all events so far
  int getMax();

  /// Get the minimum duration of all events so far
  int getMin();

  /// Get the total duration of all events so far
  int getTotal();

  /// Get the number of all events so far
  int getCount();

  /// Get the time percentage that the total time of this event took in relation to the globalDuration.
  int getTimePercentage(Event::Clock::duration globalDuration);

  /// Aggregated properties for this event
  using Properties = std::map<std::string, double>;
  Properties properties;

private:
  int count = 0;
  Event::Clock::duration total = Event::Clock::duration::zero();
  Event::Clock::duration max   = Event::Clock::duration::min();
  Event::Clock::duration min   = Event::Clock::duration::max();

};


/// High level object that stores data of all events.
/** Call EventRegistry::intialize at the beginning of your application and
EventRegistry::finalize at the end. Event timings will be usuable without calling this
function at all, but global timings as well as percentages do not work this way.  */
class EventRegistry
{
public:
  /// Sets the global start time
  static void initialize();

  /// Sets the global end time
  static void finalize();

  /// clears the registry. needed for tests
  static void clear();

  /// Finalize the timings and call print. Can be used as a crash handler to still get some timing results.
  static void signal_handler(int signal);

  /// Records the event.
  static void put(Event* event);

  /// Prints the result table to an arbitrary stream.
  /** terse enabled a more machine readable format with one event per line, seperated by whitespace. */
  static void print(std::ostream &out, bool terse = false);

  /// Convenience function: Prints to std::cout
  static void print(bool terse = false);

  /// Convenience function: Prints to filename
  static void print(std::string filename, bool terse = false);

  static void printGlobalDuration();

private:
  static bool initialized;
  static Event::Clock::time_point globalStart;
  static Event::Clock::time_point globalStop;
  static std::map<std::string, EventData> events;
};

/// Convenience function that calls EventRegistry::initalize
void Events_Init();

/// Convenience function that calls EventRegistry::finalize
void Events_Finalize();

/// Convenience function that calls EventRegistry::clear
void Events_Clear();

}} // namespace precice::utils

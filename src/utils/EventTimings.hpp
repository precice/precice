#pragma once

#include <chrono>
#include <list>
#include <map>
#include <vector>
#include <string>
#include "logging/Logger.hpp"

namespace precice {
namespace utils {

/// Represents an event that can be started and stopped.
/** Additionally to the duration there is a special property that can be set for a event.
A property is a a key-value pair with a numerical value that can be used to trace certain events,
like MPI calls in an event. It is intended to be set by the user. */
class Event
{
public:
  
  enum class State {
    STOPPED = 0,
    STARTED = 1,
    PAUSED  = 2,
  };

  /// Default clock type. All other chrono types are derived from it.
  using Clock = std::chrono::steady_clock;

  using StateChanges = std::vector<std::tuple<State, Clock::time_point>>;
    
  /// An Event can't be copied.
  Event(const Event & other) = delete;
  
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
  Clock::duration getDuration();

  std::vector<int> data;

  StateChanges stateChanges;

private:
  
  Clock::time_point starttime;
  // Clock::time_point stoptime;
  Clock::duration duration = Clock::duration::zero();
  State state = State::STOPPED;
  bool _barrier = false;
  static logging::Logger _log;
  
};



/// Class that aggregates durations for a specific event.
class EventData
{
public:
  // Do not add explicit here, it fails on some (older?) compilers
  EventData(std::string _name);
  
  EventData(std::string _name, int _rank, long _count, long _total,
            long _max, long _min, std::vector<int> _data, Event::StateChanges stateChanges);
  
  /// Adds an Events data.
  void put(Event* event);

  std::string getName() const;

  /// Get the average duration of all events so far.
  long getAvg() const;

  /// Get the maximum duration of all events so far
  long getMax() const;

  /// Get the minimum duration of all events so far
  long getMin() const;

  /// Get the total duration of all events so far
  long getTotal() const;

  /// Get the number of all events so far
  long getCount() const;

  /// Get the time percentage that the total time of this event took in relation to the global duration.
  int getTimePercentage() const;

  const std::vector<int> & getData() const;

  void print(std::ostream &out);

  void writeCSV(std::ostream &out);

  void writeEventLog(std::ostream &out);

  Event::Clock::duration max = Event::Clock::duration::min();
  Event::Clock::duration min = Event::Clock::duration::max();
  int rank;
  Event::StateChanges stateChanges;
  
private:
  std::string name;
  long count = 0;
  Event::Clock::duration total = Event::Clock::duration::zero();
  std::vector<int> data;
};

/// Holds data aggregated from all MPI ranks for one event
struct GlobalEventStats
{
  // std::string name;
  int maxRank, minRank;
  Event::Clock::duration max   = Event::Clock::duration::min();
  Event::Clock::duration min   = Event::Clock::duration::max();
};


using GlobalEvents = std::multimap<std::string, EventData>;

/// High level object that stores data of all events.
/** Call EventRegistry::intialize at the beginning of your application and
EventRegistry::finalize at the end. Event timings will be usuable without calling this
function at all, but global timings as well as percentages do not work this way.  */
class EventRegistry
{
public:
  /// Deleted copy operator for singleton pattern
  EventRegistry(EventRegistry const &) = delete;
  
  /// Deleted assigment operator for singleton pattern
  void operator=(EventRegistry const &) = delete;

  static EventRegistry & instance();
  
  /// Sets the global start time
  /**
   * @param[in] applicationName A name that is added to the logfile to distinguish different participants
   */
  void initialize(std::string appName = "");

  /// Sets the global end time
  void finalize();

  /// Clears the registry. needed for tests
  void clear();

  /// Finalizes the timings and calls print. Can be used as a crash handler to still get some timing results.
  void signal_handler(int signal);

  /// Records the event.
  void put(Event* event);

  /// Make this returning a reference or smart ptr?
  Event & getStoredEvent(std::string name);

  /// Returns the timestamp of the run, i.e. when the run finished
  std::chrono::system_clock::time_point getTimestamp();
  
  /// Returns the duration of the run in ms, either currently running, or fixed when run is stopped.
  Event::Clock::duration getDuration();

  /// Prints a verbose report to stdout and a terse one to EventTimings-AppName.log
  void printAll();

  /// Prints the result table to an arbitrary stream.
  /** terse enables a more machine readable format with one event per line, seperated by whitespace. */
  void print(std::ostream &out);

  /// Convenience function: Prints to std::cout
  void print();

  void writeCSV(std::string filename);

  void writeEventLogs(std::string filename);
  
  void printGlobalStats();

private:
  /// Private, empty constructor for singleton pattern
  EventRegistry()
    : globalEvent("_GLOBAL", true, false) // Unstarted, it's started in initialize
  {}
  
  /// Gather EventData from all ranks on rank 0.
  void collect();
  
  /// Returns length of longest name
  size_t getMaxNameWidth();

  /// Event for measuring global time, also acts as a barrier
  Event globalEvent;
  
  bool initialized = false;
  Event::Clock::time_point starttime;
  Event::Clock::duration duration;
  
  /// Timestamp when the run finished
  std::chrono::system_clock::time_point timestamp;

  /// Map of name -> events for this rank only
  std::map<std::string, EventData> events;

  std::map<std::string, Event> storedEvents;

  /// Multimap of name -> EventData of events for all ranks
  GlobalEvents globalEvents;

  /// A name that is added to the logfile to distinguish different participants
  std::string applicationName;
};
}} // namespace precice::utils

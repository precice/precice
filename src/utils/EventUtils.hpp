#pragma once

#include <chrono>
#include <iosfwd>
#include <map>
#include <stddef.h>
#include <string>
#include <utility>
#include <vector>
#include "Event.hpp"

#ifndef PRECICE_NO_MPI
#include <mpi.h>
#else
#include "utils/MPI_Mock.hpp"
#endif

namespace precice {
namespace utils {

/// Class that aggregates durations for a specific event.
class EventData {
public:
  explicit EventData(std::string _name);

  EventData(std::string _name, long _count, long _total, long _max, long _min,
            Event::Data data, Event::StateChanges stateChanges);

  /// Adds an Events data.
  void put(Event const &event);

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

  Event::Data const &getData() const;

  Event::Clock::duration max   = Event::Clock::duration::min();
  Event::Clock::duration min   = Event::Clock::duration::max();
  Event::Clock::duration total = Event::Clock::duration::zero();

  Event::StateChanges stateChanges;

private:
  std::string                             name;
  long                                    count = 0;
  std::map<std::string, std::vector<int>> data;
};

/// Holds all EventData of one particular rank
class RankData {
public:
  /// Records the initialized timestamp
  void initialize();

  /// Records the finalized timestamp
  void finalize();

  /// Adds a new event
  void put(Event const &event);

  /// Adds aggregated data for a specific event
  void addEventData(EventData ed);

  /// Normalizes all Events to zero time of t0
  void normalizeTo(std::chrono::system_clock::time_point t0);

  /// Clears all Event data
  void clear();

  /// Map of EventName -> EventData, should be private later
  std::map<std::string, EventData> evData;

  std::chrono::system_clock::duration getDuration() const;

  std::chrono::system_clock::time_point initializedAt;
  std::chrono::system_clock::time_point finalizedAt;

private:
  std::chrono::steady_clock::time_point initializedAtTicks;
  std::chrono::steady_clock::time_point finalizedAtTicks;

  bool isFinalized = true;
};

/// Holds data aggregated from all MPI ranks for one event
struct GlobalEventStats {
  int                    maxRank, minRank;
  Event::Clock::duration max = Event::Clock::duration::min();
  Event::Clock::duration min = Event::Clock::duration::max();
};

/// High level object that stores data of all events.
/** Call EventRegistry::intialize at the beginning of your application and
EventRegistry::finalize at the end. Event timings will be usuable without calling this
function at all, but global timings as well as percentages do not work this way.  */
class EventRegistry {
public:
  /// Deleted copy operator for singleton pattern
  EventRegistry(EventRegistry const &) = delete;

  /// Deleted assigment operator for singleton pattern
  void operator=(EventRegistry const &) = delete;

  /// Returns the only instance (singleton) of the EventRegistry class
  static EventRegistry &instance();

  /// Sets the global start time
  /**
   * @param[in] applicationName A name that is added to the logfile to distinguish different participants
   * @param[in] runName A name of the run, will be printed as a separate column with each Event.
   * @param[in] comm MPI communicator which is used for barriers and collecting information from ranks.
   */
  void initialize(std::string applicationName = "", std::string runName = "", MPI_Comm comm = MPI_COMM_WORLD);

  /// Sets the global end time
  void finalize();

  /// Clears the registry. needed for tests
  void clear();

  /// Finalizes the timings and calls print. Can be used as a crash handler to still get some timing results.
  void signal_handler(int signal);

  /// Records the event.
  void put(Event const &event);

  /// Returns or creates a stored event, i.e., an event with life beyond the current scope
  Event &getStoredEvent(std::string const &name);

  /// Prints a pretty report to stdout and a JSON report to appName-events.json
  void printAll() const;

  /// Prints the result table to an arbitrary stream, only prints at rank 0.
  void writeSummary(std::ostream &out) const;

  /// Writes the aggregated timings and state changes at JSON, only at rank 0.
  void writeJSON(std::ostream &out) const;

  MPI_Comm const &getMPIComm() const;

  /// Currently active prefix. Changing that applies only to newly created events.
  std::string prefix;

  /// A name that is added to the logfile to identify a run
  std::string runName;

private:
  /// Private, empty constructor for singleton pattern
  EventRegistry()
      : globalEvent("_GLOBAL", true, false) // Unstarted, it's started in initialize
  {
  }

  RankData localRankData;

  /// Holds RankData from all ranks, only populated at rank 0
  std::vector<RankData> globalRankData;

  /// Gather EventData from all ranks on rank 0.
  void collect();

  /// Normalize times among all ranks
  void normalize();

  /// Collects first initialize and last finalize time at rank 0.
  std::pair<std::chrono::system_clock::time_point, std::chrono::system_clock::time_point> collectInitAndFinalize();

  /// Returns length of longest name
  size_t getMaxNameWidth() const;

  /// Finds the first initialized time and last finalized time in globalRankData
  std::pair<std::chrono::system_clock::time_point, std::chrono::system_clock::time_point> findFirstAndLastTime() const;

  /// Event for measuring global time
  Event globalEvent;

  bool initialized = false;

  bool finalized = false;

  std::map<std::string, Event> storedEvents;

  /// A name that is added to the logfile to distinguish different participants
  std::string applicationName;

  /// MPI Communicator
  MPI_Comm comm;
};

} // namespace utils
} // namespace precice

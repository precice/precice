#pragma once

#include <chrono>
#include <cstddef>
#include <fstream>
#include <iosfwd>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "logging/Logger.hpp"
#include "profiling/Event.hpp"
#include "utils/assertion.hpp"

namespace precice::profiling {

/// The Mode of the Event utility
enum struct Mode {
  All,
  Fundamental,
  Off
};

enum struct EventClass : bool {
  Normal      = false,
  Fundamental = true
};

inline EventClass toEventClass(bool isFundamental)
{
  return static_cast<EventClass>(isFundamental);
}

/// An event that has been recorded and it waiting to be written to file
struct TimedEntry {
  TimedEntry(int eid, Event::Clock::time_point c)
      : eid(eid), clock(c) {}

  int                      eid;
  Event::Clock::time_point clock;
};

struct StartEntry : TimedEntry {
  using TimedEntry::TimedEntry;
  static constexpr char type = 'b';
};

struct StopEntry : TimedEntry {
  using TimedEntry::TimedEntry;
  static constexpr char type = 'e';
};

struct DataEntry : TimedEntry {
  DataEntry(int eid, Event::Clock::time_point c, int did, int dv)
      : TimedEntry(eid, c), did(did), dvalue(dv) {}

  static constexpr char type = 'd';
  int                   did;
  int                   dvalue;
};

struct NameEntry {
  static constexpr char type = 'n';
  std::string           name;
  int                   id;
};

using PendingEntry = std::variant<StartEntry, StopEntry, DataEntry, NameEntry>;

/** High level object that stores data of all events.
 *
 * Call EventRegistry::initialize at the beginning of your application and
 * EventRegistry::finalize at the end.
 *
 * Use \ref setWriteQueueMax() to adjust buffering behaviour.
 */
class EventRegistry {
public:
  ~EventRegistry();

  /// Deleted copy and move SMFs for singleton pattern
  EventRegistry(EventRegistry const &) = delete;
  EventRegistry(EventRegistry &&)      = delete;
  EventRegistry &operator=(EventRegistry const &) = delete;
  EventRegistry &operator=(EventRegistry &&) = delete;

  /// Returns the only instance (singleton) of the EventRegistry class
  static EventRegistry &instance();

  /// Sets the global start time
  /**
   * @param[in] applicationName A name that is added to the logfile to distinguish different participants
   * @param[in] rank the current number of the parallel instance
   * @param[in] size the total number of a parallel instances
   */
  void initialize(std::string_view applicationName, int rank = 0, int size = 1);

  /// Sets the maximum size of the writequeue before calling flush(). Use 0 to flush on destruction.
  void setWriteQueueMax(std::size_t size);

  /// Sets the directory where to write the event files to
  void setDirectory(std::string_view directory);

  /// Sets the operational mode of the registry.
  void setMode(Mode mode);

  /// Create the file and starts the filestream if profiling is turned on
  void startBackend();

  /// Sets the global end time and flushes buffers
  void finalize();

  /// Clears the registry.
  void clear();

  /// Records an event
  void put(PendingEntry pe);

  /// Writes all recorded events to file and flushes the buffer.
  void flush();

  /// Should an event of this class be forwarded to the registry?
  inline bool accepting(EventClass ec) const
  {
    return _mode == Mode::All || (ec == EventClass::Fundamental && _mode == Mode::Fundamental);
  }

  int nameToID(std::string_view name);

  /// Currently active prefix. Changing that applies only to newly created events.
  std::string prefix;

private:
  /// The name of the current participant
  std::string _applicationName;

  std::string _directory;

  /// The operational mode of the registry
  Mode _mode = Mode::Fundamental;

  /// The rank/number of parallel instance of the current program
  int _rank = 0;

  /// The amount of parallel instances of the current program
  int _size = 1;

  /// Indicator for the first record to be written
  bool _firstwrite = true;

  /// The id of the global event
  std::optional<int> _globalId;

  /// Private, empty constructor for singleton pattern
  EventRegistry() = default;

  std::map<std::string, int, std::less<>> _nameDict;

  std::vector<PendingEntry> _writeQueue;
  std::size_t               _writeQueueMax = 0;

  std::ofstream _output;

  bool _initialized = false;

  bool _finalized = false;

  bool _isBackendRunning = false;

  /// The initial time clock, used to take runtime measurements.
  Event::Clock::time_point _initClock;

  /// The initial time, used to describe when the run started.
  std::chrono::system_clock::time_point _initTime;

  /// Stops the global event, flushes the buffers and closes the filestream
  void stopBackend();

  logging::Logger _log{"Events"};
};

} // namespace precice::profiling

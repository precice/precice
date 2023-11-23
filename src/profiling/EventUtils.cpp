#include <algorithm>
#include <array>
#include <boost/filesystem/operations.hpp>
#include <cassert>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <memory>
#include <optional>
#include <ratio>
#include <string>
#include <string_view>
#include <sys/types.h>
#include <tuple>
#include <utility>
#include <variant>

#include "logging/LogMacros.hpp"
#include "profiling/Event.hpp"
#include "profiling/EventUtils.hpp"
#include "utils/assertion.hpp"
#include "utils/fmt.hpp"

namespace precice::profiling {

using sys_clk  = std::chrono::system_clock;
using stdy_clk = std::chrono::steady_clock;

/// Converts the time_point into a string like "2019-01-10T18:30:46.834"
std::string timepoint_to_string(sys_clk::time_point c)
{
  using namespace std::chrono;
  std::time_t ts = sys_clk::to_time_t(c);
  auto        ms = duration_cast<microseconds>(c.time_since_epoch()) % 1000;

  std::stringstream ss;
  ss << std::put_time(std::localtime(&ts), "%FT%T") << "." << std::setw(3) << std::setfill('0') << ms.count();
  return ss.str();
}

// -----------------------------------------------------------------------

EventRegistry::~EventRegistry()
{
  finalize();
}

EventRegistry &EventRegistry::instance()
{
  static EventRegistry instance;
  return instance;
}

void EventRegistry::initialize(std::string_view applicationName, int rank, int size)
{
  auto initClock = Event::Clock::now();
  auto initTime  = std::chrono::system_clock::now();

  this->_applicationName = std::move(applicationName);
  this->_rank            = rank;
  this->_size            = size;
  this->_initTime        = initTime;
  this->_initClock       = initClock;

  _writeQueue.clear();
  _firstwrite = true;
  _globalId   = std::nullopt;

  _initialized = true;
  _finalized   = false;
}

void EventRegistry::setWriteQueueMax(std::size_t size)
{
  _writeQueueMax = size;
}

void EventRegistry::setDirectory(std::string_view directory)
{
  _directory = directory;
}

void EventRegistry::setMode(Mode mode)
{
  _mode = mode;
}

namespace {
std::string toString(Mode m)
{
  switch (m) {
  case (Mode::Off):
    return "off";
  case (Mode::Fundamental):
    return "fundamental";
  case (Mode::All):
    return "all";
  }
  PRECICE_UNREACHABLE("Unknown mode");
}
} // namespace

void EventRegistry::startBackend()
{
  if (_mode == Mode::Off) {
    PRECICE_DEBUG("Profiling is turned off. Backend will not start.");
    return;
  }
  // Create the directory if necessary
  bool isLocal = _directory.empty() || _directory == ".";
  if (!isLocal) {
    auto exists = boost::filesystem::exists(_directory);
    PRECICE_CHECK(
        !(exists && !boost::filesystem::is_directory(_directory)),
        "The destination folder \"{}\" exists but isn't a directory. Please remove the directory \"precice-run\" and try again.",
        _directory);
    if (!exists) {
      boost::filesystem::create_directories(_directory);
    }
  }
  auto filename = fmt::format("{}/{}-{}-{}.json", _directory, _applicationName, _rank, _size);
  PRECICE_DEBUG("Starting backend with events-file: \"{}\"", filename);
  _output.open(filename);
  PRECICE_CHECK(_output, "Unable to open the events-file: \"{}\"", filename);
  _globalId = nameToID("_GLOBAL");
  PRECICE_ASSERT(_globalId.has_value());
  _writeQueue.emplace_back(StartEntry{_globalId.value(), _initClock});

  // write header
  fmt::print(_output,
             R"({{
  "meta":{{
  "name": "{}",
  "rank": "{}",
  "size": "{}",
  "unix_us": "{}",
  "tinit": "{}",
  "mode": "{}"
  }},
  "events":[
  )",
             _applicationName,
             _rank,
             _size,
             std::chrono::duration_cast<std::chrono::microseconds>(_initTime.time_since_epoch()).count(),
             timepoint_to_string(_initTime),
             toString(_mode));
  _output.flush();
}

void EventRegistry::stopBackend()
{
  if (_mode == Mode::Off) {
    return;
  }
  // create end of global event
  auto now = Event::Clock::now();
  if (_globalId) {
    put(StopEntry{_globalId.value(), now});
  }
  // flush the queue
  flush();
  _output << "]}";
  _output.close();
  _nameDict.clear();
}

void EventRegistry::finalize()
{
  if (_finalized)
    return;

  stopBackend();

  _initialized = false;
  _finalized   = true;
}

void EventRegistry::clear()
{
  _writeQueue.clear();
}

void EventRegistry::put(PendingEntry pe)
{
  PRECICE_ASSERT(_mode != Mode::Off, "The profiling is off.")
  _writeQueue.emplace_back(std::move(pe));
  if (_writeQueueMax > 0 && _writeQueue.size() > _writeQueueMax) {
    flush();
  }
}

namespace {
struct EventWriter {
  std::ostream &           out;
  Event::Clock::time_point initClock;
  std::string              prefix;

  auto sinceInit(Event::Clock::time_point tp)
  {
    return std::chrono::duration_cast<std::chrono::microseconds>(tp - initClock).count();
  }

  void operator()(const StartEntry &se)
  {
    fmt::print(out,
               R"({}{{"et":"{}","eid":{},"ts":{}}})",
               prefix, se.type, se.eid, sinceInit(se.clock));
  }

  void operator()(const StopEntry &se)
  {
    fmt::print(out,
               R"({}{{"et":"{}","eid":{},"ts":{}}})",
               prefix, se.type, se.eid, sinceInit(se.clock));
  }

  void operator()(const DataEntry &de)
  {
    fmt::print(out,
               R"({}{{"et":"{}","eid":{},"ts":{},"dn":{},"dv":"{}"}})",
               prefix, de.type, de.eid, sinceInit(de.clock), de.did, de.dvalue);
  }

  void operator()(const NameEntry &ne)
  {
    fmt::print(out,
               R"({}{{"et":"n","en":"{}","eid":{}}})",
               prefix, ne.name, ne.id);
  }
};
} // namespace

void EventRegistry::flush()
try {
  if (_mode == Mode::Off || _writeQueue.empty()) {
    return;
  }
  PRECICE_ASSERT(_output, "Filestream doesn't exist.");

  auto first = _writeQueue.begin();
  // Don't prefix the first write with a comma
  if (_firstwrite) {
    PRECICE_ASSERT(!_writeQueue.empty() && !_writeQueue.front().valueless_by_exception());
    std::visit(EventWriter{_output, _initClock, ""}, _writeQueue.front());
    ++first;
    _firstwrite = false;
  }

  EventWriter ew{_output, _initClock, ","};
  std::for_each(first, _writeQueue.end(), [&ew](const auto &pe) { std::visit(ew, pe); });

  _output.flush();
  _writeQueue.clear();
} catch (const std::bad_variant_access &e) {
  PRECICE_UNREACHABLE(e.what());
}

int EventRegistry::nameToID(std::string_view name)
{
  if (auto iter = _nameDict.find(name);
      iter == _nameDict.end()) {
    int id = _nameDict.size();
    _nameDict.insert(iter, {std::string(name), id});
    _writeQueue.emplace_back(NameEntry{std::string(name), id});
    return id;
  } else {
    return iter->second;
  }
}

} // namespace precice::profiling

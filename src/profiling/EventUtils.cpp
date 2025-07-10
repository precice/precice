#include <algorithm>
#include <array>
#include <cassert>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <lzma.h>
#include <memory>
#include <optional>
#include <ratio>
#include <sstream>
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

/// The version of the Events file. Increase on changes
constexpr int file_version{1};

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
  _nameDict.clear();

  _globalId = nameToID("_GLOBAL");
  _writeQueue.emplace_back(StartEntry{_globalId.value(), _initClock});

  _initialized = true;
  _finalized   = false;
  if (_isBackendRunning) {
    stopBackend();
  }
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
  case (Mode::API):
    return "api";
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

  PRECICE_ASSERT(!_isBackendRunning);

  // Create the directory if necessary
  bool isLocal = _directory.empty() || _directory == ".";
  if (!isLocal) {
    auto exists = std::filesystem::exists(_directory);
    PRECICE_CHECK(
        !(exists && !std::filesystem::is_directory(_directory)),
        "The destination folder \"{}\" exists but isn't a directory. Please remove the directory \"precice-run\" and try again.",
        _directory);
    if (!exists) {
      std::filesystem::create_directories(_directory);
    }
  }
  auto filename = fmt::format("{}/{}-{}-{}.txt", _directory, _applicationName, _rank, _size);
  PRECICE_DEBUG("Starting backend with events-file: \"{}\"", filename);
  _output.open(filename);
  PRECICE_CHECK(_output, "Unable to open the events-file: \"{}\"", filename);

  // write header
  fmt::println(_output,
               R"({{"name":"{}","rank":{},"size":{},"unix_us":"{}","tinit":"{}","mode":"{}","file_version":{}}})",
               _applicationName,
               _rank,
               _size,
               std::chrono::duration_cast<std::chrono::microseconds>(_initTime.time_since_epoch()).count(),
               timepoint_to_string(_initTime),
               toString(_mode),
               ::precice::profiling::file_version);

  _strm = LZMA_STREAM_INIT;
  lzma_options_lzma options;
  lzma_lzma_preset(&options, 0);
  options.mode = LZMA_MODE_FAST;
  if (lzma_alone_encoder(&_strm, &options) != LZMA_OK) {
    throw std::runtime_error("Failed to init LZMA encoder");
  }

  _output.flush();
  _isBackendRunning = true;
}

void EventRegistry::stopBackend()
{
  if (_mode == Mode::Off || !_isBackendRunning) {
    return;
  }
  // create end of global event
  auto now = Event::Clock::now();
  put(StopEntry{*_globalId, now});
  // flush the queue
  flush();

  lzma_ret ret;
  do {
    _strm.next_out  = reinterpret_cast<uint8_t *>(_buf.data());
    _strm.avail_out = _buf.size();

    ret = lzma_code(&_strm, LZMA_FINISH);
    if (ret != LZMA_OK && ret != LZMA_STREAM_END) {
      throw std::runtime_error("Finalization failed");
    }
    _output.write(_buf.data(), _buf.size() - _strm.avail_out);
  } while (ret != LZMA_STREAM_END);

  lzma_end(&_strm);

  _output.close();
  _nameDict.clear();

  _isBackendRunning = false;
}

void EventRegistry::finalize()
{
  if (_finalized || !_initialized) {
    return;
  }

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
  PRECICE_ASSERT(_mode != Mode::Off, "The profiling is off.");

  // avoid flushing the queue when we start measuring but only if we don't explicitly want to write every entry
  auto skipFlush = _writeQueueMax != 1 && std::holds_alternative<StartEntry>(pe);

  _writeQueue.emplace_back(std::move(pe));
  if (!skipFlush && _writeQueueMax > 0 && _writeQueue.size() > _writeQueueMax) {
    flush();
  }
}

void EventRegistry::putCritical(PendingEntry pe)
{
  PRECICE_ASSERT(_mode != Mode::Off, "The profiling is off.");
  _writeQueue.emplace_back(std::move(pe));
}

namespace {
struct EventWriter {
  std::ostream            &out;
  Event::Clock::time_point initClock;

  auto sinceInit(Event::Clock::time_point tp)
  {
    return std::chrono::duration_cast<std::chrono::microseconds>(tp - initClock).count();
  }

  void operator()(const StartEntry &se)
  {
    fmt::print(out,
               "B{}:{}\n",
               se.eid, sinceInit(se.clock));
  }

  void operator()(const StopEntry &se)
  {
    fmt::print(out,
               "E{}:{}\n",
               se.eid, sinceInit(se.clock));
  }

  void operator()(const DataEntry &de)
  {
    fmt::print(out,
               "D{}:{}:{}:{}\n",
               de.eid, sinceInit(de.clock), de.did, de.dvalue);
  }

  void operator()(const NameEntry &ne)
  {
    fmt::print(out,
               "N{}:{}\n",
               ne.id, ne.name);
  }
};
} // namespace

void EventRegistry::flush()
try {
  if (_mode == Mode::Off || _writeQueue.empty()) {
    return;
  }
  PRECICE_ASSERT(_output, "Filestream doesn't exist.");

  std::ostringstream oss{};
  EventWriter        ew{oss, _initClock};
  std::for_each(_writeQueue.begin(), _writeQueue.end(), [&ew](const auto &pe) { std::visit(ew, pe); });
  auto str = oss.str();

  _strm.next_in  = reinterpret_cast<uint8_t *>(str.data());
  _strm.avail_in = str.size();

  while (_strm.avail_in > 0) {
    _strm.next_out  = reinterpret_cast<uint8_t *>(_buf.data());
    _strm.avail_out = _buf.size();
    lzma_ret ret    = lzma_code(&_strm, LZMA_RUN);
    PRECICE_ASSERT(ret == LZMA_OK, "Compression failed");
    _output.write(_buf.data(), _buf.size() - _strm.avail_out);
  }

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

#include "Event.hpp"
#include "EventUtils.hpp"
#include "logging/LogMacros.hpp"

namespace precice {
namespace utils {

Event::Event(std::string eventName, Clock::duration initialDuration)
    : name(EventRegistry::instance().prefix + eventName),
      duration(initialDuration)
{
  EventRegistry::instance().put(*this);
}

Event::Event(std::string eventName, bool barrier, bool autostart)
    : name(eventName),
      _barrier(barrier)
{
  // Set prefix here: workaround to omit data lock between instance() and Event ctor
  if (eventName != "_GLOBAL")
    name = EventRegistry::instance().prefix + eventName;
  if (autostart) {
    start(_barrier);
  }
}

Event::~Event()
{
  stop(_barrier);
}

void Event::start(bool barrier)
{
  if (barrier)
    MPI_Barrier(EventRegistry::instance().getMPIComm());

  state = State::STARTED;
  stateChanges.push_back(std::make_pair(State::STARTED, Clock::now()));
  starttime = Clock::now();
  PRECICE_DEBUG("Started event " << name);
}

void Event::stop(bool barrier)
{
  if (state == State::STARTED or state == State::PAUSED) {
    if (barrier)
      MPI_Barrier(EventRegistry::instance().getMPIComm());

    if (state == State::STARTED) {
      auto stoptime = Clock::now();
      duration += Clock::duration(stoptime - starttime);
    }
    stateChanges.push_back(std::make_pair(State::STOPPED, Clock::now()));
    state = State::STOPPED;
    EventRegistry::instance().put(*this);
    data.clear();
    stateChanges.clear();
    duration = Clock::duration::zero();
    PRECICE_DEBUG("Stopped event " << name);
  }
}

void Event::pause(bool barrier)
{
  if (state == State::STARTED) {
    if (barrier)
      MPI_Barrier(EventRegistry::instance().getMPIComm());

    auto stoptime = Clock::now();
    stateChanges.emplace_back(State::PAUSED, Clock::now());
    state = State::PAUSED;
    duration += Clock::duration(stoptime - starttime);
    PRECICE_DEBUG("Paused event " << name);
  }
}

Event::Clock::duration Event::getDuration() const
{
  return duration;
}

void Event::addData(std::string key, int value)
{
  data[key].push_back(value);
}

// -----------------------------------------------------------------------

ScopedEventPrefix::ScopedEventPrefix(std::string const &name)
{
  previousName = EventRegistry::instance().prefix;
  EventRegistry::instance().prefix += name;
}

ScopedEventPrefix::~ScopedEventPrefix()
{
  EventRegistry::instance().prefix = previousName;
}

} // namespace utils
} // namespace precice

#include "profiling/Event.hpp"
#include "profiling/EventUtils.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

namespace precice::profiling {

Event::Event(std::string_view eventName, Options options)
    : _fundamental(options.fundamental), _synchronize(options.synchronized)
{
  auto &er = EventRegistry::instance();
  _eid     = er.nameToID(std::string(EventRegistry::instance().prefix).append(eventName));
  start();
}

Event::~Event()
{
  if (_state == State::RUNNING) {
    stop();
  }
}

void Event::start()
{
  if (_synchronize) {
    ::precice::utils::IntraComm::synchronize();
  }
  auto timestamp = Clock::now();
  PRECICE_ASSERT(_state == State::STOPPED, _eid);
  _state = State::RUNNING;

  if (EventRegistry::instance().accepting(toEventClass(_fundamental))) {
    EventRegistry::instance().put(StartEntry{_eid, timestamp});
  }
}

void Event::stop()
{
  auto timestamp = Clock::now();
  PRECICE_ASSERT(_state == State::RUNNING, _eid);
  _state = State::STOPPED;

  if (EventRegistry::instance().accepting(toEventClass(_fundamental))) {
    EventRegistry::instance().put(StopEntry{_eid, timestamp});
  }
}

void Event::addData(std::string_view key, int value)
{
  auto timestamp = Clock::now();
  PRECICE_ASSERT(_state == State::RUNNING, _eid);

  auto &er = EventRegistry::instance();
  if (er.accepting(toEventClass(_fundamental))) {
    auto did = er.nameToID(key);
    er.put(DataEntry{_eid, timestamp, did, value});
  }
}

// -----------------------------------------------------------------------

ScopedEventPrefix::ScopedEventPrefix(std::string_view name)
{
  previousName = EventRegistry::instance().prefix;
  EventRegistry::instance().prefix += name;
}

ScopedEventPrefix::~ScopedEventPrefix()
{
  pop();
}

void ScopedEventPrefix::pop()
{
  EventRegistry::instance().prefix = previousName;
}

} // namespace precice::profiling

#include "profiling/Event.hpp"
#include "profiling/EventUtils.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

namespace precice::profiling {

Event::Event(std::string_view eventName, Options options)
    : _fundamental(options.fundamental), _synchronize(options.synchronized)
{
  auto &er   = EventRegistry::instance();
  auto  name = std::string(eventName);
  PRECICE_ASSERT(eventName.find('/') == std::string_view::npos);
  _eid = er.nameToID(name);
  if (_synchronize && er.parallel()) {
    _sid = er.nameToID(name + ".sync");
  }
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
  PRECICE_ASSERT(_state == State::STOPPED, _eid);
  _state         = State::RUNNING;
  auto &registry = EventRegistry::instance();
  if (!registry.accepting(toEventClass(_fundamental))) {
    return;
  }

  if (_synchronize && registry.parallel() && ::precice::utils::IntraComm::willSynchronize()) {
    // We need to synchronize, so we record a sync event
    PRECICE_ASSERT(_sid != -1);
    registry.putCritical(StartEntry{_sid, Clock::now()});
    ::precice::utils::IntraComm::synchronize();
    // end of sync
    auto timestamp = Clock::now();
    registry.putCritical(StopEntry{_sid, timestamp});
    registry.put(StartEntry{_eid, timestamp});
  } else {
    registry.put(StartEntry{_eid, Clock::now()});
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

} // namespace precice::profiling

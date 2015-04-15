#include "EventTimings.hpp"

#include "MasterSlave.hpp"
#include "Parallel.hpp"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <ctime>
#include <vector>
#include <map>
#include <chrono>

namespace precice {
namespace utils {

Event::Event(std::string eventName, Clock::duration eventDuration)
    : name(eventName)
    , duration(eventDuration)
    , isStarted(false)
    , _barrier(false) {
  EventRegistry::put(this);
}

Event::Event(std::string eventName, bool barrier, bool autostart)
    : _barrier(barrier) {
  name = eventName;

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
    Parallel::synchronizeProcesses();

  isStarted = true;
  starttime = Clock::now();
}

void Event::stop(bool barrier)
{
  if (isStarted) {
    if (barrier)
      Parallel::synchronizeProcesses();

    stoptime = Clock::now();
    isStarted = false;
    duration = Clock::duration(stoptime - starttime);
    EventRegistry::put(this);
  }
}

void Event::addProp(std::string property, double value)
{
  properties[property] += value;
}


Event::Clock::duration Event::getDuration()
{
  return duration;
}

// -----------------------------------------------------------------------


void EventData::put(Event* event)
{
  if (not precice::utils::MasterSlave::_slaveMode) {
    count++;
    Event::Clock::duration duration = event->getDuration();
    total += duration;
    min = std::min(duration, min);
    max = std::max(duration, max);

    for (auto p : event->properties) {
      properties[p.first] += p.second;
    }
  }
}

int EventData::getAvg()
{
  return (std::chrono::duration_cast<std::chrono::milliseconds>(total) / count).count();

}

int EventData::getMax()
{
  return std::chrono::duration_cast<std::chrono::milliseconds>(max).count();
}

int EventData::getMin()
{
  return std::chrono::duration_cast<std::chrono::milliseconds>(min).count();
}

int EventData::getTotal()
{
  return std::chrono::duration_cast<std::chrono::milliseconds>(total).count();
}

int EventData::getCount()
{
  return count;
}

int EventData::getTimePercentage(Event::Clock::duration globalDuration)
{
  return ((double) total.count() / globalDuration.count()) * 100;
}

// -----------------------------------------------------------------------

// Static members need to be initalized like that
std::map<std::string, EventData> EventRegistry::events;
Event::Clock::time_point EventRegistry::globalStart;
Event::Clock::time_point EventRegistry::globalStop;
bool EventRegistry::initialized = false;

void EventRegistry::initialize()
{
  globalStart = Event::Clock::now();
  initialized = true;
}

void EventRegistry::finalize()
{
  globalStop = Event::Clock::now();
  initialized = false;
}

void EventRegistry::clear()
{
  events.clear();
}

void EventRegistry::signal_handler(int signal)
{
  if (initialized) {
    finalize();
    print();
  }
}

void EventRegistry::put(Event* event)
{
  EventData data = events[event->name];
  data.put(event);
  events[event->name] = data;
}

void EventRegistry::print(std::ostream &out, bool terse)
{
  if (not precice::utils::MasterSlave::_slaveMode) {
    using std::endl;
    using std::setw; using std::setprecision;
    using std::left; using std::right;
    EventData::Properties allProps;
    Event::Clock::duration globalDuration = globalStop - globalStart;

    std::time_t currentTime = std::time(NULL);
    out << "Run finished at " << std::asctime(std::localtime(&currentTime));

    if (not terse) {
      out << "Global runtime = "
          << std::chrono::duration_cast<std::chrono::milliseconds>(globalDuration).count() << "ms / "
          << std::chrono::duration_cast<std::chrono::seconds>(globalDuration).count() << "s"
          << "\n" << "\n";

      out << "Event                Count    Total[ms]     Max[ms]     Min[ms]     Avg[ms]   T%" << "\n";
      out << "--------------------------------------------------------------------------------" << "\n";

      for (auto e : events) {
        out << setw(14) << left << e.first << right
            << setw(12) << e.second.getCount()
            << setw(12) << e.second.getTotal()
            << setw(12) << e.second.getMax()
            << setw(12) << e.second.getMin()
            << setw(12) << e.second.getAvg()
            << setw(6)  << e.second.getTimePercentage(globalDuration)
            << "\n";
        for (auto p : e.second.properties) {
          allProps[p.first] += p.second;

          out << "  " << setw(12) << left << p.first
              << setw(12) << right << std::fixed << std::setprecision(5) << p.second
              << "\n";
        }
        out << "\n";
      }

      out << "Properties from all Events, accumulated" << "\n";
      out << "---------------------------------------" << "\n";
      for (auto a : allProps) {
        out << setw(14) << left << a.first << right
            << setw(12) << std::fixed << std::setprecision(5) << a.second << "\n";
      }
    }
    else // terse output
    {
      auto global = std::chrono::duration_cast<std::chrono::milliseconds>(globalDuration).count();

      out << "# Eventname Count Total Max Min Avg T%" << "\n";
      out << "\"GLOBAL\" "  << 1 << " "        // Eventname Count
          << global << " "  << global << " "   // Total Max
          <<  global << " "  << global << " "  // Min Avg
          << 100 << "\n";                      // T%
      for (auto e : events) {
        out << "\"" << e.first << "\" "
            << e.second.getCount() << " " << e.second.getTotal() << " "
            << e.second.getMax()   << " " << e.second.getMin()   << " "
            << e.second.getAvg()   << " " << e.second.getTimePercentage(globalDuration) << "\n";
      }
    }

    out << endl;
  }
}

void EventRegistry::print(bool terse)
{
  EventRegistry::print(std::cout, terse);
}

void EventRegistry::print(std::string filename, bool terse)
{
  std::ofstream outfile;
  outfile.open(filename, std::ios::out | std::ios::app);
  EventRegistry::print(outfile, terse);
  outfile.close();
}

void EventRegistry::printGlobalDuration()
{
  Event::Clock::duration globalDuration = Event::Clock::now() - globalStart;

  std::cout << "Global Duration = "
            << std::chrono::duration_cast<std::chrono::milliseconds>(
                   globalDuration).count() << "ms" << std::endl;
}

void Events_Init()
{
  EventRegistry::initialize();
}

void Events_Finalize()
{
  EventRegistry::finalize();
}

void Events_Clear()
{
  EventRegistry::clear();
}

}} // namespace precice::utils

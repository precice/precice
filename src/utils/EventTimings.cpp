#include "EventTimings.hpp"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <ctime>
#include <vector>
#include <map>
#include <chrono>

#include "utils/MasterSlave.hpp"

namespace precice {
namespace utils {

Event::Event(std::string eventName, bool autostart)
{
  name = eventName;
  if (autostart) {
    isStarted = true;
    starttime = Clock::now();
  }
}

Event::~Event()
{
  stop();
}

void Event::start()
{
  isStarted = true;
  starttime = Clock::now();
}

void Event::stop()
{
  if (isStarted) {
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

void EventRegistry::print(std::ostream &out)
{
  if (not precice::utils::MasterSlave::_slaveMode) {
    using std::endl;
    using std::setw; using std::setprecision;
    using std::left; using std::right;
    EventData::Properties allProps;
    Event::Clock::duration globalDuration = globalStop - globalStart;

    std::time_t currentTime = std::time(NULL);
    out << "Run finished at " << std::asctime(std::localtime(&currentTime));

    out << "Global runtime = "
        << std::chrono::duration_cast<std::chrono::milliseconds>(globalDuration).count() << "ms / "
        << std::chrono::duration_cast<std::chrono::seconds>(globalDuration).count() << "s"
        << endl << endl;

    out << "Event                Count    Total[ms]     Max[ms]     Min[ms]     Avg[ms]   T%" << endl;
    out << "--------------------------------------------------------------------------------" << endl;
        
    for (auto e : events) {
      out << setw(14) << left << e.first << right 
          << setw(12) << e.second.getCount()
          << setw(12) << e.second.getTotal()
          << setw(12) << e.second.getMax() 
          << setw(12) << e.second.getMin()
          << setw(12) << e.second.getAvg()
          << setw(6)  << e.second.getTimePercentage(globalDuration)
          << endl;
      for (auto p : e.second.properties) {
        allProps[p.first] += p.second;
      
        out << "  " << setw(12) << left << p.first
            << setw(12) << right << std::fixed << std::setprecision(5) << p.second
            << endl;
      }
      out << endl;
    }

    out << "All Events, accumulated" << endl;
    out << "--------------------------" << endl;
    for (auto a : allProps) {
      out << setw(14) << left << a.first << right
          << setw(12) << std::fixed << std::setprecision(5) << a.second << endl;
    }
  }
}

void EventRegistry::print()
{
  EventRegistry::print(std::cout);
}

void EventRegistry::print(std::string filename)
{
  std::ofstream outfile;
  outfile.open(filename, std::ios::out | std::ios::app);
  EventRegistry::print(outfile);
  outfile.close();
}

void Events_Init()
{
  EventRegistry::initialize();
}

void Events_Finalize()
{
  EventRegistry::finalize();
}

}} // namespace precice::utils

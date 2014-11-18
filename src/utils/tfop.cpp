#include "tfop.hpp"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <chrono>

using namespace std;


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
  if (isStarted) {
    std::cerr << "Event " << name << " started which was already started. Restarting." << std::endl;
  }
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

void Event::addProp(std::string property, int value)
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
  count++;
  Event::Clock::duration duration = event->getDuration();
  total += duration;
  if (min > duration) {
    min = duration;
  }
  if (max < duration) {
    max = duration;
  }
  for (auto p : event->properties) {
    properties[p.first] += p.second;
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

void EventRegistry::print()
{
  using std::cout; using std::endl;
  using std::setw; using std::setprecision;
  using std::left; using std::right;
  EventData::Properties allProps;
  Event::Clock::duration globalDuration = globalStop - globalStart;

  cout << "Global runtime = "
       << std::chrono::duration_cast<std::chrono::milliseconds>(globalDuration).count() << "ms / "
       << std::chrono::duration_cast<std::chrono::seconds>(globalDuration).count() << "s"
       << endl << endl;

  cout << "Event                Count    Total[ms]     Max[ms]     Min[ms]     Avg[ms]   T%" << endl;
  cout << "--------------------------------------------------------------------------------" << endl;
        
  for (auto e : events) {
    cout << setw(14) << left << e.first << right 
         << setw(12) << e.second.getCount()
         << setw(12) << e.second.getTotal()
         << setw(12) << e.second.getMax() 
         << setw(12) << e.second.getMin()
         << setw(12) << e.second.getAvg()
         << setw(6) << e.second.getTimePercentage(globalDuration)
         << endl;
    for (auto p : e.second.properties) {
      allProps[p.first] += p.second;
      
      cout << "  " << setw(12) << left << p.first
           << setw(12) << right << p.second
           << endl;
    }
    cout << endl;
  }

  cout << "All Events, accumulated" << endl;
  cout << "--------------------------" << endl;
  for (auto a : allProps) {
    cout << setw(14) << left << a.first << right
         << setw(12) << a.second << endl;
  }
  
}



void TFOP_Init()
{
  std::cout << "Initialize TFOP" << std::endl;
  EventRegistry::initialize();
  
}

void TFOP_Finalize()
{
  std::cout << "Finalize TLOP" << std::endl;
  EventRegistry::finalize();
}

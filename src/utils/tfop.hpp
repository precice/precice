#pragma once

#include <chrono>
#include <map>
#include <string>


class Event
{
public:
  using Clock = std::chrono::high_resolution_clock;
  std::string name;  
  Event(std::string eventName, bool autostart = true);
  ~Event();
  void start();
  void stop();
  Clock::duration getDuration();

  std::map<std::string, double> properties;
  // void setProp(std::string property, double value);
  // void setProp(std::string property, int value);

  // void addProp(std::string property, double value);
  void addProp(std::string property, int value);

  // void incProp(std::string property);
  // void incProp(std::string property);

  

private:  
  Clock::time_point starttime;
  Clock::time_point stoptime;
  Clock::duration duration;
  bool isStarted = false;

};

class EventData
{
public:
  void put(Event* event);
  int getAvg();
  int getMax();
  int getMin();
  int getTotal();
  int getCount();
  int getTimePercentage(Event::Clock::duration globalDuration);

  using Properties = std::map<std::string, double>;
  Properties properties; // evtl. friend machen
  
private:
  int count = 0;
  Event::Clock::duration total = Event::Clock::duration::zero();
  Event::Clock::duration max   = Event::Clock::duration::min();
  Event::Clock::duration min   = Event::Clock::duration::max();

};



class EventRegistry
{
public:
  static std::map<std::string, EventData> events;
  static void put(Event* event);
  static Event::Clock::time_point globalStart;
  static Event::Clock::time_point globalStop;
  static void print(int signal);
};


void TFOP_Init();
void TFOP_Finalize();

#include "EventTimings.hpp"

#include "MasterSlave.hpp"
#include "Parallel.hpp"
#include <cassert>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstring>
#include <ctime>
#include <vector>
#include <map>
#include <chrono>
#include <utility>
#include <limits>
#ifndef PRECICE_NO_MPI
#include <mpi.h>
#endif
#include "prettyprint.hpp"
#include "TableWriter.hpp"

namespace precice {
namespace utils {

logging::Logger Event::_log("utils::Events");

struct MPI_EventData
{
  char name[255] = {'\0'};
  int rank, count = 0;
  long total = 0, max = 0, min = 0;
  int dataSize = 0, stateChangesSize = 0;
};

Event::Event(std::string eventName, Clock::duration initialDuration)
  : name(eventName),
    duration(initialDuration)
{
  EventRegistry::instance().put(this);
}

Event::Event(std::string eventName, bool barrier, bool autostart)
  : name(eventName),
    _barrier(barrier)
{
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
    
  state = State::STARTED;
  stateChanges.push_back(std::make_tuple(State::STARTED, Clock::now()));
  starttime = Clock::now();
  DEBUG("Started event " << name);
}

void Event::stop(bool barrier)
{
  if (state == State::STARTED or state == State::PAUSED) {
    if (barrier)
      Parallel::synchronizeProcesses();

    if (state == State::STARTED) {
      auto stoptime = Clock::now();
      duration += Clock::duration(stoptime - starttime);
    }
    stateChanges.push_back(std::make_tuple(State::STOPPED, Clock::now()));
    state = State::STOPPED;
    EventRegistry::instance().put(this);
    data.clear();
    duration = Clock::duration::zero();
    DEBUG("Stopped event " << name);
  }
}

void Event::pause(bool barrier)
{
  if (state == State::STARTED) {
    if (barrier)
      Parallel::synchronizeProcesses();

    auto stoptime = Clock::now();
    stateChanges.push_back(std::make_tuple(State::PAUSED, Clock::now()));
    state = State::PAUSED;
    duration += Clock::duration(stoptime - starttime);
    DEBUG("Paused event " << name);
  }
}

Event::Clock::duration Event::getDuration()
{
  return duration;
}

// -----------------------------------------------------------------------

EventData::EventData(std::string _name) :
  name(_name)
{
  rank =  Parallel::getProcessRank();
}

EventData::EventData(std::string _name, int _rank, long _count, long _total,
                     long _max, long _min, std::vector<int> _data, Event::StateChanges _stateChanges)
  :  max(std::chrono::milliseconds(_max)),
     min(std::chrono::milliseconds(_min)),
     rank(_rank),
     stateChanges(_stateChanges),
     name(_name),
     count(_count),
     total(std::chrono::milliseconds(_total)),
     data(_data)
{}

void EventData::put(Event* event)
{
  count++;
  Event::Clock::duration duration = event->getDuration();
  total += duration;
  min = std::min(duration, min);
  max = std::max(duration, max);
  data.insert(std::end(data), std::begin(event->data), std::end(event->data));
  stateChanges.insert(std::end(stateChanges), std::begin(event->stateChanges), std::end(event->stateChanges));
}

std::string EventData::getName() const
{
  return name;
}

long EventData::getAvg() const
{
  return (std::chrono::duration_cast<std::chrono::milliseconds>(total) / count).count();
}

long EventData::getMax() const
{
  return std::chrono::duration_cast<std::chrono::milliseconds>(max).count();
}

long EventData::getMin() const
{
  return std::chrono::duration_cast<std::chrono::milliseconds>(min).count();
}

long EventData::getTotal() const
{
  return std::chrono::duration_cast<std::chrono::milliseconds>(total).count();
}

long EventData::getCount() const
{
  return count;
}

int EventData::getTimePercentage() const
{
  return (static_cast<double>(total.count()) / EventRegistry::instance().getDuration().count()) * 100;
}

const std::vector<int> & EventData::getData() const
{
  return data;
}


void EventData::print(std::ostream &out)
{
  using std::setw;
  out << setw(30) << std::left << name << std::right
      << setw(12) << getCount()
      << setw(12) << getTotal()
      << setw(12) << getMax()
      << setw(12) << getMin()
      << setw(12) << getAvg()
      << setw(10) << getTimePercentage()
      << "\n";
}

void EventData::writeCSV(std::ostream &out)
{
  std::time_t ts = std::chrono::system_clock::to_time_t(EventRegistry::instance().getTimestamp());
      
  out << std::put_time(std::localtime(&ts), "%F %T") << ","
      << rank << ","
      << getName() << ","
      << getCount() << ","
      << getTotal() << ","
      << getMin() << "," << getMax() << "," << getAvg() << ","
      << getTimePercentage();

  bool first = true;
  out << ",\"[";
  for (auto & d : data) {
    if (not first)
      out << ",";
    out << d;
    first = false;
  }
  out << "]\"" << std::endl;
}

void EventData::writeEventLog(std::ostream &out)
{
  using namespace std::chrono;
  std::time_t its = system_clock::to_time_t(EventRegistry::instance().getTimestamp());
  for (auto & sc : stateChanges) {
    out << std::put_time(std::localtime(&its), "%F %T") << ","
        << name << ","
        << rank << ","
        << duration_cast<milliseconds>(std::get<1>(sc).time_since_epoch()).count() << ","
        << static_cast<int>(std::get<0>(sc)) << std::endl;
  }
}


// -----------------------------------------------------------------------

EventRegistry & EventRegistry::instance()
{
  static EventRegistry instance;
  return instance;
}

void EventRegistry::initialize(std::string appName)
{
  starttime = Event::Clock::now();
  applicationName = appName;

  globalEvent.start(true);
  initialized = true;
}

void EventRegistry::finalize()
{
  globalEvent.stop(true); // acts as a barrier
  duration = Event::Clock::now() - starttime;
  timestamp = std::chrono::system_clock::now();
  initialized = false;
  for (auto & e : storedEvents)
    e.second.stop();
  collect();
}

void EventRegistry::clear()
{
  events.clear();
}

void EventRegistry::signal_handler(int signal)
{
  if (initialized) {
    finalize();
    printAll();
  }
}

void EventRegistry::put(Event* event)
{
  auto data = std::get<0>(events.emplace(event->name, event->name));
  data->second.put(event);
}

Event & EventRegistry::getStoredEvent(std::string name)
{
  auto insertion = storedEvents.emplace(std::piecewise_construct,
                                        std::forward_as_tuple(name),
                                        std::forward_as_tuple(name, false, false));

  return std::get<0>(insertion)->second;
}


std::chrono::system_clock::time_point EventRegistry::getTimestamp()
{
  return timestamp;
}

Event::Clock::duration EventRegistry::getDuration()
{
  if (duration == Event::Clock::duration::zero())
    return Event::Clock::now() - starttime;
  else
      return duration;
}

void EventRegistry::printAll()
{
  print();

  std::string csvFile, logFile;
  if (applicationName.empty()) {
    csvFile = "EventTimings.log";
    logFile = "Events.log";
  }
  else {
    csvFile = "EventTimings-" + applicationName + ".log";
    logFile = "Events-" + applicationName + ".log";
  }
  writeCSV(csvFile);
  writeEventLogs(logFile);
  
}


void EventRegistry::print(std::ostream &out)
{
  int rank =  Parallel::getProcessRank();
  int size = Parallel::getCommunicatorSize();
  
  if (rank == 0) {
    using std::endl;
    using std::setw; using std::setprecision;
    using std::left; using std::right;
    
    std::time_t ts = std::chrono::system_clock::to_time_t(timestamp);
    auto globalDuration = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
    
    // out << "Run finished at " << std::put_time(std::localtime(&currentTime), "%c") << std::endl;
    out << "Run finished at " << std::asctime(std::localtime(&ts));

    out << "Global runtime       = "
        << globalDuration << "ms / "
        << std::chrono::duration_cast<std::chrono::seconds>(duration).count() << "s" << std::endl
        << "Number of processors = " << size << std::endl
        << "# Rank: " << rank << std::endl << std::endl;

    Table table( {
        {getMaxNameWidth(), "Event"},
        {10, "Count"},
        {10, "Total[ms]"},
        {10, "Max[ms]"},
        {10, "Min[ms]"},
        {10, "Avg[ms]"},
        {10, "T[%]"}
      });
    table.printHeader();
      
    for (auto & e : events) {
      auto & ev = e.second;
      table.printLine(ev.getName(), ev.getCount(), ev.getTotal(), ev.getMax(),
                      ev.getMin(), ev.getAvg(), ev.getTimePercentage());
    }        

    out << endl;
    printGlobalStats();
    out << endl << std::flush;
  }
}

void EventRegistry::print()
{
  EventRegistry::print(std::cout);
}

void EventRegistry::writeCSV(std::string filename)
{
  int rank =  Parallel::getProcessRank();
  int size = Parallel::getCommunicatorSize();
  
  if (rank != 0)
    return;
  
  bool fileExists = std::ifstream(filename).is_open();
   
  std::ofstream outfile;
  outfile.open(filename, std::ios::out | std::ios::app);
  if (not fileExists)
    outfile << "Timestamp,Rank,Name,Count,Total,Min,Max,Avg,T%,Data" << "\n";
   
  std::time_t ts = std::chrono::system_clock::to_time_t(timestamp);
  // auto globalDuration = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
  std::tm tm = *std::localtime(&ts);

  outfile << "# Run finished at: " << std::put_time(&tm, "%F %T") << std::endl
          << "# Number of processors: " << size << std::endl;
  
  // table.printLine("GLOBAL", 0, 1, globalDuration, globalDuration, globalDuration, globalDuration, 100);
      
  for (auto e : globalEvents) {
    auto & ev = e.second;
    ev.writeCSV(outfile);
  }
  
  outfile.close();
}

void EventRegistry::writeEventLogs(std::string filename)
{
  int rank =  Parallel::getProcessRank();
    
  if (rank != 0)
    return;
  
  bool fileExists = std::ifstream(filename).is_open();
   
  std::ofstream outfile;
  outfile.open(filename, std::ios::out | std::ios::app);
  if (not fileExists)
    outfile << "RunTimestamp,Name,Rank,Timestamp,State" << "\n";
   
  for (auto e : globalEvents) {
    auto & ev = e.second;
    ev.writeEventLog(outfile);
  }
  
  outfile.close();
}


 
std::map<std::string, GlobalEventStats> getGlobalStats(GlobalEvents events)
{
  std::map<std::string, GlobalEventStats> globalStats;
  for (auto & e : events) {
    auto ev = e.second;
    GlobalEventStats & stats = globalStats[e.first];
    if (ev.max > stats.max) {
      stats.max = ev.max;
      stats.maxRank = ev.rank;
    }
    if (ev.min < stats.min) {
      stats.min = ev.min;
      stats.minRank = ev.rank;
    }    
  }
  return globalStats;
}

void EventRegistry::printGlobalStats()
{
  Table t({ {getMaxNameWidth(), "Name"},
      {10, "Max"}, {10, "MaxOnRank"}, {10, "Min"}, {10, "MinOnRank"}, {10, "Min/Max"} });
  t.printHeader();
  
  auto stats = getGlobalStats(globalEvents);
  for (auto & e : stats) {
    auto & ev = e.second;
    double rel = static_cast<double>(ev.min.count()) / static_cast<double>(ev.max.count());
    t.printLine(e.first, ev.max, ev.maxRank, ev.min, ev.minRank, rel);
  }
}


void EventRegistry::collect()
{
  #ifndef PRECICE_NO_MPI
  // Register MPI datatype
  MPI_Datatype MPI_EVENTDATA;
  int blocklengths[] = {255, 2, 3, 2};
  MPI_Aint displacements[] = {offsetof(MPI_EventData, name), offsetof(MPI_EventData, rank),
                              offsetof(MPI_EventData, total), offsetof(MPI_EventData, dataSize)};
  MPI_Datatype types[] = {MPI_CHAR, MPI_INT, MPI_LONG, MPI_INT};
  MPI_Type_create_struct(4, blocklengths, displacements, types, &MPI_EVENTDATA);
  MPI_Type_commit(&MPI_EVENTDATA);
 
  int rank, MPIsize;
  MPI_Comm_rank(Parallel::getGlobalCommunicator(), &rank);
  MPI_Comm_size(Parallel::getGlobalCommunicator(), &MPIsize);
  
  std::vector<MPI_Request> requests;
  std::vector<int> eventsPerRank(MPIsize);
  size_t eventsSize = events.size();
  MPI_Gather(&eventsSize, 1, MPI_INT, eventsPerRank.data(), 1, MPI_INT, 0, Parallel::getGlobalCommunicator());

  std::vector<MPI_EventData> eventSendBuf(events.size());
  std::vector<std::unique_ptr<char[]>> packSendBuf(events.size());
  int i = 0;
  for (const auto & ev : events) {
    MPI_EventData eventdata;
    MPI_Request req;
    
    assert(ev.first.size() <= 255);
    ev.first.copy(eventSendBuf[i].name, 255);
    eventSendBuf[i].rank = rank;
    eventSendBuf[i].count = ev.second.getCount();
    eventSendBuf[i].total = ev.second.getTotal();
    eventSendBuf[i].max = ev.second.getMax();
    eventSendBuf[i].min = ev.second.getMin();
    eventSendBuf[i].dataSize = ev.second.getData().size();
    eventSendBuf[i].stateChangesSize = ev.second.stateChanges.size();
    
    int packSize = 0, pSize = 0;
    // int packSize = sizeof(int) * ev.second.getData().size() +
      // sizeof(Event::StateChanges::value_type) * ev.second.stateChanges.size();
    MPI_Pack_size(ev.second.getData().size(), MPI_INT, Parallel::getGlobalCommunicator(), &pSize);
    packSize += pSize;
    MPI_Pack_size(ev.second.stateChanges.size() * sizeof(Event::StateChanges::value_type),
                  MPI_BYTE, Parallel::getGlobalCommunicator(), &pSize);
    packSize += pSize;
    
    packSendBuf[i] = std::unique_ptr<char[]>(new char[packSize]);
    int position = 0;
    MPI_Pack(const_cast<int*>(ev.second.getData().data()), ev.second.getData().size(),
             MPI_INT, packSendBuf[i].get(), packSize, &position, Parallel::getGlobalCommunicator());
    MPI_Pack(const_cast<Event::StateChanges::pointer>(ev.second.stateChanges.data()),
             ev.second.stateChanges.size() * sizeof(Event::StateChanges::value_type),
             MPI_BYTE, packSendBuf[i].get(), packSize, &position, Parallel::getGlobalCommunicator());

    MPI_Isend(&eventSendBuf[i], 1, MPI_EVENTDATA, 0, 0, Parallel::getGlobalCommunicator(), &req);
    requests.push_back(req);
    MPI_Isend(packSendBuf[i].get(), position, MPI_PACKED, 0, 0, Parallel::getGlobalCommunicator(), &req);
    requests.push_back(req);
    ++i;
  }
  
  if (rank == 0) {
    for (int i = 0; i < MPIsize; ++i) {
      for (int j = 0; j < eventsPerRank[i]; ++j) {
        MPI_EventData ev;
        MPI_Recv(&ev, 1, MPI_EVENTDATA, i, MPI_ANY_TAG, Parallel::getGlobalCommunicator(), MPI_STATUS_IGNORE);
        std::vector<int> recvData(ev.dataSize);
        Event::StateChanges recvStateChanges(ev.stateChangesSize);
        MPI_Status status;
        int packSize = 0, position = 0;
        MPI_Probe(i, MPI_ANY_TAG, Parallel::getGlobalCommunicator(), &status);
        MPI_Get_count(&status, MPI_PACKED, &packSize);
        char packBuffer[packSize];
        MPI_Recv(packBuffer, packSize, MPI_PACKED, i, MPI_ANY_TAG,
                 Parallel::getGlobalCommunicator(), MPI_STATUS_IGNORE);
        MPI_Unpack(packBuffer, packSize, &position, recvData.data(), ev.dataSize, MPI_INT,
                   Parallel::getGlobalCommunicator());
        MPI_Unpack(packBuffer, packSize, &position, recvStateChanges.data(),
                   ev.stateChangesSize * sizeof(Event::StateChanges::value_type), MPI_BYTE,
                   Parallel::getGlobalCommunicator());
          
        globalEvents.emplace(std::piecewise_construct, std::forward_as_tuple(ev.name),
                             std::forward_as_tuple(ev.name, ev.rank, ev.count, ev.total, ev.max, ev.min,
                                                   recvData, recvStateChanges));

      }
    }
  }    
  MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
  MPI_Type_free(&MPI_EVENTDATA);
  #endif
}


size_t EventRegistry::getMaxNameWidth()
{
  size_t maxEventWidth = 0;
  for (auto & ev : events)
    if (ev.second.getName().size() > maxEventWidth)
      maxEventWidth = ev.second.getName().size();
  
  return maxEventWidth;
}

}} // namespace precice::utils

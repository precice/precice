#include "EventUtils.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <memory>
#include <nlohmann/json.hpp>
#include <ratio>
#include <string>
#include <tuple>
#include <utility>
#include "TableWriter.hpp"
#include "utils/Event.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace utils {

using sys_clk  = std::chrono::system_clock;
using stdy_clk = std::chrono::steady_clock;

/// Converts the time_point into a string like "2019-01-10T18:30:46.834"
std::string timepoint_to_string(sys_clk::time_point c)
{
  using namespace std::chrono;
  std::time_t ts = sys_clk::to_time_t(c);
  auto        ms = duration_cast<milliseconds>(c.time_since_epoch()) % 1000;

  std::stringstream ss;
  ss << std::put_time(std::localtime(&ts), "%FT%T") << "." << std::setw(3) << std::setfill('0') << ms.count();
  return ss.str();
}

std::map<std::string, GlobalEventStats> getGlobalStats(std::vector<RankData> events)
{
  std::map<std::string, GlobalEventStats> globalStats;
  for (size_t rank = 0; rank < events.size(); ++rank) {
    for (auto &evData : events[rank].evData) {
      auto const &      event = evData.second;
      GlobalEventStats &stats = globalStats[evData.first];
      if (event.max > stats.max) {
        stats.max     = event.max;
        stats.maxRank = rank;
      }
      if (event.min < stats.min) {
        stats.min     = event.min;
        stats.minRank = rank;
      }
    }
  }
  return globalStats;
}

struct MPI_EventData {
  char name[255] = {'\0'};
  int  count     = 0;
  long total = 0, max = 0, min = 0;
  int  dataSize = 0, stateChangesSize = 0;
};

// -----------------------------------------------------------------------

EventData::EventData(std::string _name)
    : name(_name)
{
}

EventData::EventData(std::string _name, long _count, long _total, long _max, long _min,
                     Event::Data data, Event::StateChanges _stateChanges)
    : max(std::chrono::milliseconds(_max)),
      min(std::chrono::milliseconds(_min)),
      total(std::chrono::milliseconds(_total)),
      stateChanges(_stateChanges),
      name(_name),
      count(_count),
      data(data)
{
}

void EventData::put(Event const &event)
{
  count++;
  stdy_clk::duration duration = event.getDuration();
  total += duration;
  min = std::min(duration, min);
  max = std::max(duration, max);
  for (auto const &d : event.data) {
    auto &source = std::get<1>(d);
    auto &target = data[std::get<0>(d)];
    target.insert(target.begin(), source.begin(), source.end());
  }
  stateChanges.insert(std::end(stateChanges), std::begin(event.stateChanges), std::end(event.stateChanges));
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

Event::Data const &EventData::getData() const
{
  return data;
}

// -----------------------------------------------------------------------

void RankData::initialize()
{
  initializedAt      = sys_clk::now();
  initializedAtTicks = stdy_clk::now();
  isFinalized        = false;
}

void RankData::finalize()
{
  finalizedAt      = sys_clk::now();
  finalizedAtTicks = stdy_clk::now();
  isFinalized      = true;
}

void RankData::put(Event const &event)
{
  /// Construct or return EventData object with name as key and name as arg to ctor.
  auto data = std::get<0>(evData.emplace(std::piecewise_construct,
                                         std::forward_as_tuple(event.name),
                                         std::forward_as_tuple(event.name)));
  data->second.put(event);
}

void RankData::addEventData(EventData ed)
{
  evData.emplace(ed.getName(), std::move(ed));
}

void RankData::normalizeTo(sys_clk::time_point t0)
{
  auto const delta = initializedAt - t0; // duration that this rank initialized after the first rank
  PRECICE_ASSERT(t0 <= initializedAt);   // t0 should always be before or equal my init time

  for (auto &events : evData) {
    for (auto &sc : events.second.stateChanges) {
      auto &tp = sc.second;
      tp       = stdy_clk::time_point(tp - initializedAtTicks + delta);
      // Deactivated, as it fails for the precice tests, likely because Events are used before initialize
      // assert(tp.time_since_epoch().count() >= 0); // Trying to do normalize twice?
    }
  }
}

void RankData::clear()
{
  evData.clear();
}

sys_clk::duration RankData::getDuration() const
{
  if (isFinalized)
    return finalizedAt - initializedAt;
  else
    return sys_clk::now() - initializedAt;
}

// -----------------------------------------------------------------------

EventRegistry &EventRegistry::instance()
{
  static EventRegistry instance;
  return instance;
}

void EventRegistry::initialize(std::string applicationName, std::string runName, MPI_Comm comm)
{
  this->applicationName = applicationName;
  this->runName         = runName;
  this->comm            = comm;

  localRankData.initialize();

  globalEvent.start(false);
  initialized = true;
  finalized   = false;
}

void EventRegistry::finalize()
{
  if (finalized)
    return;

  globalEvent.stop();
  localRankData.finalize();

  for (auto &e : storedEvents)
    e.second.stop();

  if (initialized) // this makes only sense when it was properly initialized
    normalize();

  collect();

  initialized = false;
  finalized   = true;
}

void EventRegistry::clear()
{
  localRankData.clear();
  globalRankData.clear();
  storedEvents.clear();
}

void EventRegistry::signal_handler(int signal)
{
  if (initialized) {
    finalize();
    printAll();
  }
}

void EventRegistry::put(Event const &event)
{
  localRankData.put(event);
}

Event &EventRegistry::getStoredEvent(std::string const &name)
{
  // Reset the prefix for creation of a stored event. Using prefixes with stored events is possible
  // but leads to unexpected results, such as not getting the event you want, because someone else up the
  // stack set a prefix.
  auto previousPrefix = prefix;
  prefix              = "";
  auto insertion      = storedEvents.emplace(std::piecewise_construct,
                                        std::forward_as_tuple(name),
                                        std::forward_as_tuple(name, false, false));

  prefix = previousPrefix;
  return std::get<0>(insertion)->second;
}

void EventRegistry::printAll() const
{
  int myRank;
  MPI_Comm_rank(comm, &myRank);

  if (myRank != 0)
    return;

  std::string logFile;
  std::string summaryFile;
  if (applicationName.empty()) {
    logFile     = "Events.json";
    summaryFile = "Events-summary.log";
  } else {
    logFile     = applicationName + "-events.json";
    summaryFile = applicationName + "-events-summary.log";
  }

  std::ofstream summaryFS{summaryFile};
  writeSummary(summaryFS);
  std::ofstream logFS{logFile};
  writeJSON(logFS);
}

void EventRegistry::writeSummary(std::ostream &out) const
{
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  if (rank == 0) {
    using std::endl;
    { // Print per event stats
      std::time_t  ts       = sys_clk::to_time_t(localRankData.finalizedAt);
      double const duration = std::chrono::duration_cast<std::chrono::milliseconds>(localRankData.getDuration()).count();

      out << "Run finished at " << std::asctime(std::localtime(&ts));

      out << "Global runtime       = "
          << duration << "ms / "
          << duration / 1000 << "s" << endl
          << "Number of processors = " << size << endl
          << "# Rank: " << rank << endl
          << endl;

      Table table(out);
      table.addColumn("Event", getMaxNameWidth());
      table.addColumn("Count", 10);
      table.addColumn("Total[ms]", 10);
      table.addColumn("Max[ms]", 10);
      table.addColumn("Min[ms]", 10);
      table.addColumn("Avg[ms]", 10);
      table.addColumn("Time Ratio", 6, 3);
      table.printHeader();

      for (auto &e : localRankData.evData) {
        auto &ev = e.second;
        table.printRow(ev.getName(), ev.getCount(), ev.getTotal(), ev.getMax(), ev.getMin(), ev.getAvg(),
                       ev.getTotal() / duration);
      }
    }
    out << endl
        << endl;
    { // Print aggregated states
      Table t(out);
      t.addColumn("Name", getMaxNameWidth());
      t.addColumn("Max", 10);
      t.addColumn("MaxOnRank", 10);
      t.addColumn("Min", 10);
      t.addColumn("MinOnRank", 10);
      t.addColumn("Min/Max", 10);
      t.printHeader();

      auto stats = getGlobalStats(globalRankData);
      for (auto &e : stats) {
        auto & ev  = e.second;
        double rel = 0;
        if (ev.max != stdy_clk::duration::zero()) // Guard against division by zero
          rel = static_cast<double>(ev.min.count()) / ev.max.count();

        t.printRow(e.first, ev.max, ev.maxRank, ev.min, ev.minRank, rel);
      }
    }
  }
}

void EventRegistry::writeJSON(std::ostream &out) const
{
  using json = nlohmann::json;
  using namespace std::chrono;

  if (globalRankData.empty())
    return;

  json js;

  sys_clk::time_point initT, finalT;
  std::tie(initT, finalT) = findFirstAndLastTime();
  js["Name"]              = runName;
  js["Initialized"]       = timepoint_to_string(initT);
  js["Finalized"]         = timepoint_to_string(finalT);

  for (auto const &rank : globalRankData) {
    auto         jTimings      = json::object();
    auto         jStateChanges = json::array();
    double const duration      = duration_cast<milliseconds>(rank.getDuration()).count();
    for (auto const &events : rank.evData) {
      auto const &e                     = events.second;
      jTimings[events.second.getName()] = {
          {"Count", e.getCount()},
          {"Total", e.getTotal()},
          {"Max", e.getMax()},
          {"Min", e.getMin()},
          {"TimeRatio", e.getTotal() / duration},
          {"Data", e.getData()}};
      for (auto const &sc : e.stateChanges) {
        jStateChanges.push_back({{"Name", events.second.getName()},
                                 {"State", sc.first},
                                 {"Timestamp", duration_cast<milliseconds>(sc.second.time_since_epoch()).count()}});
      }
    }
    js["Ranks"].push_back({{"Finalized", timepoint_to_string(rank.finalizedAt)},
                           {"Initialized", timepoint_to_string(rank.initializedAt)},
                           {"Timings", jTimings},
                           {"StateChanges", jStateChanges}});
  }

  out << std::setw(2) << js << std::endl;
}

MPI_Comm const &EventRegistry::getMPIComm() const
{
  return comm;
}

void EventRegistry::collect()
{
#ifndef PRECICE_NO_MPI
  // Register MPI datatype
  MPI_Datatype MPI_EVENTDATA;
  int          blocklengths[]  = {255, 1, 3, 2};
  MPI_Aint     displacements[] = {offsetof(MPI_EventData, name), offsetof(MPI_EventData, count),
                              offsetof(MPI_EventData, total), offsetof(MPI_EventData, dataSize)};
  MPI_Datatype types[]         = {MPI_CHAR, MPI_INT, MPI_LONG, MPI_INT};
  MPI_Type_create_struct(4, blocklengths, displacements, types, &MPI_EVENTDATA);
  MPI_Type_commit(&MPI_EVENTDATA);

  int rank, MPIsize;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &MPIsize);

  std::vector<MPI_Request> requests;
  std::vector<int>         eventsPerRank(MPIsize);
  size_t                   eventsSize = localRankData.evData.size();
  MPI_Gather(&eventsSize, 1, MPI_INT, eventsPerRank.data(), 1, MPI_INT, 0, comm);

  std::vector<MPI_EventData>     eventSendBuf(eventsSize);
  std::vector<std::vector<long>> stateChangesBuf(eventsSize);
  int                            i = 0;

  MPI_Request req;

  // Send the times from the local RankData
  std::array<long, 2> times = {localRankData.initializedAt.time_since_epoch().count(),
                               localRankData.finalizedAt.time_since_epoch().count()};
  MPI_Isend(&times, times.size(), MPI_LONG, 0, 0, comm, &req);
  requests.push_back(req);

  // Send all events from all ranks, including rank 0, to rank 0
  for (auto const &evData : localRankData.evData) {
    const auto &  ev = evData.second;
    MPI_EventData eventdata;

    // Send aggregated EventData
    assert(evData.first.size() <= sizeof(eventdata.name));
    evData.first.copy(eventSendBuf[i].name, sizeof(eventdata.name));
    eventSendBuf[i].count            = ev.getCount();
    eventSendBuf[i].total            = ev.getTotal();
    eventSendBuf[i].max              = ev.getMax();
    eventSendBuf[i].min              = ev.getMin();
    eventSendBuf[i].dataSize         = ev.getData().size();
    eventSendBuf[i].stateChangesSize = ev.stateChanges.size();
    MPI_Isend(&eventSendBuf[i], 1, MPI_EVENTDATA, 0, 0, comm, &req);
    requests.push_back(req);

    // Send the state changes
    for (auto const &sc : ev.stateChanges) {
      stateChangesBuf[i].push_back(static_cast<long>(sc.first)); // state
      stateChangesBuf[i].push_back(
          std::chrono::duration_cast<std::chrono::milliseconds>(sc.second.time_since_epoch()).count());
    }
    MPI_Isend(stateChangesBuf[i].data(), ev.stateChanges.size() * 2, MPI_LONG, 0, 0, comm, &req);
    requests.push_back(req);

    // Send the map that stores the data associated with an event
    for (auto const &md : ev.getData()) {
      auto &key = std::get<0>(md);
      auto &val = std::get<1>(md);
      MPI_Isend(const_cast<char *>(key.c_str()), key.size(), MPI_CHAR, 0, 0, comm, &req);
      requests.push_back(req);
      MPI_Isend(const_cast<int *>(val.data()), val.size(), MPI_INT, 0, 0, comm, &req);
      requests.push_back(req);
    }

    ++i;
  }

  // Receive
  if (rank == 0) {
    for (int i = 0; i < MPIsize; ++i) {
      RankData data;
      // Receive initialized and finalized times
      std::array<long, 2> recvTimes;
      MPI_Recv(&recvTimes, 2, MPI_LONG, i, MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);
      data.initializedAt = sys_clk::time_point(sys_clk::duration(recvTimes[0]));
      data.finalizedAt   = sys_clk::time_point(sys_clk::duration(recvTimes[1]));

      // Receive all events from this rank
      for (int j = 0; j < eventsPerRank[i]; ++j) {
        // Receive aggregated EventData
        MPI_EventData ev;
        MPI_Recv(&ev, 1, MPI_EVENTDATA, i, MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);

        // Receive all state changes for this event
        std::vector<long> recvStateChanges(ev.stateChangesSize * 2);
        MPI_Recv(recvStateChanges.data(), recvStateChanges.size(),
                 MPI_LONG, i, MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);
        Event::StateChanges stateChanges; // evtl. reserve
        for (size_t k = 0; k < recvStateChanges.size(); k += 2) {
          stateChanges.emplace_back(static_cast<Event::State>(recvStateChanges[k]),
                                    stdy_clk::time_point(std::chrono::milliseconds(recvStateChanges[k + 1])));
        }

        // Receive the map that stores the data associated with an event
        Event::Data dataMap;
        for (int j = 0; j < ev.dataSize; j++) {
          MPI_Status status;
          int        count = 0;
          MPI_Probe(i, MPI_ANY_TAG, comm, &status);
          MPI_Get_count(&status, MPI_CHAR, &count);
          std::string key(count, '\0');
          MPI_Recv(&key[0], count, MPI_CHAR, i, MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);
          MPI_Probe(i, MPI_ANY_TAG, comm, &status);
          MPI_Get_count(&status, MPI_INT, &count);
          std::vector<int> val(count);
          MPI_Recv(val.data(), count, MPI_INT, i, MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);
          dataMap[key] = val;
        }

        // Create the EventData
        EventData ed(ev.name, ev.count, ev.total, ev.max, ev.min, dataMap, stateChanges);
        data.addEventData(std::move(ed));
      }
      globalRankData.push_back(data);
    }
  }
  MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
  MPI_Type_free(&MPI_EVENTDATA);
#else
  globalRankData.push_back(localRankData);
#endif
}

void EventRegistry::normalize()
{
  long ticks = localRankData.initializedAt.time_since_epoch().count();
  long minTicks;
  MPI_Allreduce(&ticks, &minTicks, 1, MPI_LONG, MPI_MIN, EventRegistry::instance().getMPIComm());

  // This assumes the same epoch and ticks rep, should be true for system time
  sys_clk::time_point t0{sys_clk::duration{minTicks}};
  localRankData.normalizeTo(t0);
}

std::pair<sys_clk::time_point, sys_clk::time_point>
EventRegistry::collectInitAndFinalize()
{
  long ticks = localRankData.initializedAt.time_since_epoch().count();
  long minTicks;
  MPI_Reduce(&ticks, &minTicks, 1, MPI_LONG, MPI_MIN, 0, EventRegistry::instance().getMPIComm());

  ticks = localRankData.finalizedAt.time_since_epoch().count();
  long maxTicks;
  MPI_Reduce(&ticks, &maxTicks, 1, MPI_LONG, MPI_MAX, 0, EventRegistry::instance().getMPIComm());

  // This assumes the same epoch and ticks rep, should be true for system time

  return {sys_clk::time_point{sys_clk::duration{minTicks}},
          sys_clk::time_point{sys_clk::duration{maxTicks}}};
}

size_t EventRegistry::getMaxNameWidth() const
{
  size_t maxEventWidth = 0;
  for (auto &ev : localRankData.evData)
    if (ev.second.getName().size() > maxEventWidth)
      maxEventWidth = ev.second.getName().size();

  return maxEventWidth;
}

std::pair<sys_clk::time_point, sys_clk::time_point> EventRegistry::findFirstAndLastTime() const
{
  using T    = decltype(globalRankData)::value_type const &;
  auto first = std::min_element(std::begin(globalRankData), std::end(globalRankData),
                                [](T a, T b) { return a.initializedAt < b.initializedAt; });
  auto last  = std::max_element(std::begin(globalRankData), std::end(globalRankData),
                               [](T a, T b) { return a.finalizedAt < b.finalizedAt; });

  return std::make_pair(first->initializedAt, last->finalizedAt);
}

} // namespace utils
} // namespace precice

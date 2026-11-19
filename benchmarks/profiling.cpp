#include <benchmark/benchmark.h>

#include "profiling/Event.hpp"
#include "profiling/EventUtils.hpp"

using namespace precice::profiling;

// this file contains benchmarks to profile the event creation overhead of the profiling library

// profiles the time it takes to emit an known event
static void knownEvent(benchmark::State &state)
{
  auto &er = EventRegistry::instance();
  er.initialize("bench-k", 0, 1);
  er.setWriteQueueMax(0); // keep in memory
  er.setMode(Mode::All);

  std::string_view name = "event";
  // create name record
  er.nameToID(name);

  for (auto _ : state) {
    // profile the event only
    Event e(name);
  }
}

// profiles the time it takes to emit an known event
static void knownFundamentalEvent(benchmark::State &state)
{
  auto &er = EventRegistry::instance();
  er.initialize("bench-kf", 0, 1);
  er.setWriteQueueMax(0); // keep in memory
  er.setMode(Mode::All);

  std::string_view name = "event";
  // create name record
  er.nameToID(name);

  for (auto _ : state) {
    // profile the event only
    Event e(name, Fundamental);
  }
}

// profiles the time it takes to emit an unknown event into a freshly initialized registry
static void unknownEvent(benchmark::State &state)
{
  auto  &er = EventRegistry::instance();
  size_t n  = 0;
  for (auto _ : state) {
    state.PauseTiming();
    er.initialize("bench-u", 0, 1);
    er.setWriteQueueMax(0); // keep in memory
    er.setMode(Mode::All);

    std::string      sname = fmt::format("Event{}", n++);
    std::string_view name  = sname;
    state.ResumeTiming();

    Event e(name);
  }
}

// profiles the time it takes to emit an unknown fundamental event into a freshly initialized registry
static void unknownFundamentalEvent(benchmark::State &state)
{
  auto  &er = EventRegistry::instance();
  size_t n  = 0;
  for (auto _ : state) {
    state.PauseTiming();
    er.initialize("bench-uf", 0, 1);
    er.setWriteQueueMax(0); // keep in memory
    er.setMode(Mode::All);

    std::string      sname = fmt::format("Event{}", n);
    std::string_view name  = sname;
    state.ResumeTiming();

    Event e(name, Fundamental);
  }
}

// profiles the time it takes to reject an event
static void eventOff(benchmark::State &state)
{
  auto &er = EventRegistry::instance();
  er.initialize("bench-o", 0, 1);
  er.setWriteQueueMax(0); // keep in memory
  er.setMode(Mode::Off);

  std::string_view name = "event";
  for (auto _ : state) {
    Event e(name);
  }
}

// profiles the time it takes to reject a fundamental event
static void fundamentalEventOff(benchmark::State &state)
{
  auto &er = EventRegistry::instance();
  er.initialize("bench-of", 0, 1);
  er.setWriteQueueMax(0); // keep in memory
  er.setMode(Mode::Off);

  std::string_view name = "event";
  for (auto _ : state) {
    Event e(name, Fundamental);
  }
}

BENCHMARK(knownEvent)->Name("profiling known event");
BENCHMARK(knownFundamentalEvent)->Name("profiling known fundamental event");

BENCHMARK(unknownEvent)->Name("profiling new event");
BENCHMARK(unknownFundamentalEvent)->Name("profiling new fundamental event");

BENCHMARK(eventOff)->Name("profiling off");
BENCHMARK(fundamentalEventOff)->Name("profiling fundamental off");

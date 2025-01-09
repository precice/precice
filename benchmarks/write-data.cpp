#include <benchmark/benchmark.h>
#include <precice/precice.hpp>
#include <vector>

#include "helper.hpp"

/// Benchmarks writing vectorial data to n verticles of a participant before calling initialize
static void writeDataVector(benchmark::State &state)
{
  const int            nv = state.range(0);
  precice::Participant p("A", benchConfig("write-data.xml"), 0, 1);

  std::vector<int> vids(nv);
  auto             pos = generate3DHalton(nv);
  p.setMeshVertices("MeshA", pos, vids);
  std::vector<double> data(nv * 3, 3.14159);

  for (auto _ : state) {
    p.writeData("MeshA", "Vector", vids, data);
  }
  state.SetComplexityN(state.range(0));
}

BENCHMARK(writeDataVector)->Name("Write vector data to participant")->RangeMultiplier(2)->Range(1 << 10, 1 << 20)->Complexity();

/// Benchmarks writing scalar data to n verticles of a participant before calling initialize
static void writeDataScalar(benchmark::State &state)
{
  const int            nv = state.range(0);
  precice::Participant p("A", benchConfig("write-data.xml"), 0, 1);

  std::vector<int> vids(nv);
  auto             pos = generate3DHalton(nv);
  p.setMeshVertices("MeshA", pos, vids);
  std::vector<double> data(nv, 3.14159);

  for (auto _ : state) {
    p.writeData("MeshA", "Scalar", vids, data);
  }
  state.SetComplexityN(state.range(0));
}

BENCHMARK(writeDataScalar)->Name("Write scalar data to participant")->RangeMultiplier(2)->Range(1 << 10, 1 << 20)->Complexity();

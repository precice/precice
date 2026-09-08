#include <benchmark/benchmark.h>

#include "helper.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/NearestNeighborMapping.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"

using namespace precice;

namespace {

/// Two 3d Halton meshes of n vertices, offset so no vertices coincide
std::pair<mesh::PtrMesh, mesh::PtrMesh> makeMeshPair(int n)
{
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", 3, 0));
  auto          inPos = generate3DHalton(n);
  for (auto row : inPos.rowwise()) {
    inMesh->createVertex(row);
  }

  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", 3, 1));
  auto          outPos = generate3DHalton(n);
  outPos.array() += 1e-3;
  for (auto row : outPos.rowwise()) {
    outMesh->createVertex(row);
  }
  return {inMesh, outMesh};
}

} // namespace

/// Benchmarks computeMapping of a consistent nearest-neighbour mapping between two 3d meshes
static void computeNearestNeighborMapping(benchmark::State &state)
{
  int n          = state.range(0);
  auto [in, out] = makeMeshPair(n);

  for (auto _ : state) {
    state.PauseTiming();
    mapping::NearestNeighborMapping mapping(mapping::Mapping::CONSISTENT, 3);
    mapping.setMeshes(in, out);
    in->index().clear();
    state.ResumeTiming();

    mapping.computeMapping();
    benchmark::DoNotOptimize(mapping);
  }
  state.SetComplexityN(state.range(0));
}

BENCHMARK(computeNearestNeighborMapping)->Name("NN computeMapping")->RangeMultiplier(4)->Range(1 << 10, 1 << 18)->Complexity();

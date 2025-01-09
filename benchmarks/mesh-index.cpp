#include <benchmark/benchmark.h>

#include "helper.hpp"
#include "mesh/Mesh.hpp"
#include "query/Index.hpp"

using namespace precice;

static void indexVertices(benchmark::State &state)
{
  int        n = state.range(0);
  mesh::Mesh m{"A", 3, 0};
  auto       pos = generate3DHalton(n);
  for (auto row : pos.rowwise()) {
    m.createVertex(row);
  }
  for (auto _ : state) {
    query::Index index(m);
    auto         bounds = index.getRtreeBounds();
    benchmark::DoNotOptimize(bounds);
  }
  state.SetComplexityN(state.range(0));
}

BENCHMARK(indexVertices)->Name("Index vertices")->RangeMultiplier(2)->Range(1 << 10, 1 << 20)->Complexity();

static void queryNearestToCenter(benchmark::State &state)
{
  int        n = state.range(0);
  mesh::Mesh m{"A", 3, 0};
  auto       pos = generate3DHalton(n);
  for (auto row : pos.rowwise()) {
    m.createVertex(row);
  }
  auto &index = m.index();
  index.getRtreeBounds();

  Eigen::Vector3d center{0.5, 0.5, 0.5};

  for (auto _ : state) {
    auto match = index.getClosestVertex(center);
    benchmark::DoNotOptimize(match);
  }
  state.SetComplexityN(state.range(0));
}

BENCHMARK(queryNearestToCenter)->Name("Closest vertex to center")->RangeMultiplier(2)->Range(1 << 10, 1 << 20)->Complexity();

static void query10NearestToCenter(benchmark::State &state)
{
  int        n = state.range(0);
  mesh::Mesh m{"A", 3, 0};
  auto       pos = generate3DHalton(n);
  for (auto row : pos.rowwise()) {
    m.createVertex(row);
  }
  auto &index = m.index();
  index.getRtreeBounds();

  Eigen::Vector3d center{0.5, 0.5, 0.5};

  for (auto _ : state) {
    auto match = index.getClosestVertices(center, 10);
    benchmark::DoNotOptimize(match);
  }
  state.SetComplexityN(state.range(0));
}

BENCHMARK(query10NearestToCenter)->Name("10 closest vertices to center")->RangeMultiplier(2)->Range(1 << 10, 1 << 20)->Complexity();

static void queryInnerBB(benchmark::State &state)
{
  int        n = state.range(0);
  mesh::Mesh m{"A", 3, 0};
  auto       pos = generate3DHalton(n);
  for (auto row : pos.rowwise()) {
    m.createVertex(row);
  }
  auto &index = m.index();
  index.getRtreeBounds();

  Eigen::Vector3d   bbMin{0.25, 0.25, 0.25};
  Eigen::Vector3d   bbMax{0.75, 0.75, 0.75};
  mesh::BoundingBox bb{bbMin, bbMax};

  for (auto _ : state) {
    auto match = index.getVerticesInsideBox(bb);
    benchmark::DoNotOptimize(match);
  }
}

BENCHMARK(queryInnerBB)->Name("vertices inside box")->Arg(100000)->Arg(1000000);

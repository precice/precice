#include <benchmark/benchmark.h>

#include "helper.hpp"
#include "mesh/Mesh.hpp"
#include "query/Index.hpp"

using namespace precice;

static void tagInBBIndex(benchmark::State &state)
{
  int        n = state.range(0);
  mesh::Mesh m{"A", 3, 0};
  auto       pos = generate3DHalton(n);
  for (auto row : pos.rowwise()) {
    m.createVertex(row);
  }

  Eigen::Vector3d   bbMin{0.25, 0.25, 0.25};
  Eigen::Vector3d   bbMax{0.75, 0.75, 0.75};
  mesh::BoundingBox bb{bbMin, bbMax};

  for (auto _ : state) {
    query::Index index(m);
    auto         matches = index.getVerticesInsideBox(bb);
    for (auto vid : matches) {
      m.vertex(vid).tag();
    }
    benchmark::DoNotOptimize(m);
  }
}

BENCHMARK(tagInBBIndex)->Name("Tag inside box with index")->Arg(100000)->Arg(1000000);

static void tagInBBManual(benchmark::State &state)
{
  int        n = state.range(0);
  mesh::Mesh m{"A", 3, 0};
  auto       pos = generate3DHalton(n);
  for (auto row : pos.rowwise()) {
    m.createVertex(row);
  }

  Eigen::Vector3d   bbMin{0.25, 0.25, 0.25};
  Eigen::Vector3d   bbMax{0.75, 0.75, 0.75};
  mesh::BoundingBox bb{bbMin, bbMax};

  for (auto _ : state) {
    for (auto &v : m.vertices()) {
      if (bb.contains(v)) {
        v.tag();
      }
    }
    benchmark::DoNotOptimize(m);
  }
}

BENCHMARK(tagInBBManual)->Name("Tag inside box with BB")->Arg(100000)->Arg(1000000);

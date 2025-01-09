#include <benchmark/benchmark.h>

#include <Eigen/Core>
#include <mesh/BoundingBox.hpp>
#include <mesh/Vertex.hpp>

using namespace precice;

/// Benchmarks 2D BoundingBox contains with a match
static void BBcontains2D(benchmark::State &state)
{
  Eigen::Vector2d   min{-1, -3.4};
  Eigen::Vector2d   max{1, 3.4};
  mesh::BoundingBox bb{min, max};

  Eigen::Vector2d point{0, 0};
  mesh::Vertex    v(point, 0);

  for (auto _ : state) {
    bool inside = bb.contains(v);
    benchmark::DoNotOptimize(inside);
  }
}

BENCHMARK(BBcontains2D)->Name("2D BoundingBox contains");

/// Benchmarks 3D BoundingBox contains with a match
static void BBcontains3D(benchmark::State &state)
{
  Eigen::Vector3d   min{-1, -3.4, -1.0};
  Eigen::Vector3d   max{1, 3.4, 1.0};
  mesh::BoundingBox bb{min, max};

  Eigen::Vector3d point{0, 0, 0};
  mesh::Vertex    v(point, 0);

  for (auto _ : state) {
    bool inside = bb.contains(v);
    benchmark::DoNotOptimize(inside);
  }
}

BENCHMARK(BBcontains3D)->Name("3D BoundingBox contains");

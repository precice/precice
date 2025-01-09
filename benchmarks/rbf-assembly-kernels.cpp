#include <benchmark/benchmark.h>
#include <boost/range/irange.hpp>

#include "mapping/RadialBasisFctSolver.hpp"
#include "mapping/impl/BasisFunctions.hpp"

using namespace precice;

template <typename RBF>
void assembleA(benchmark::State &state, RBF f, mapping::Polynomial polynomial)
{
  std::array<bool, 3> activeAxis{true, true, true};
  const int           dimension = 3;

  mesh::Mesh inMesh("inMesh", dimension, -1);

  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 10; ++j)
      inMesh.createVertex(Eigen::Vector3d(0, i * 3.1, j * 1.4315));

  mesh::Mesh outMesh("outMesh", dimension, -1);

  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 10; ++j)
      outMesh.createVertex(Eigen::Vector3d(1, i * 1.2111, j * 5.4315));

  for (auto _ : state) {
    Eigen::MatrixXd mat = mapping::buildMatrixA(f, inMesh, boost::irange<Eigen::Index>(0, inMesh.nVertices()), outMesh, boost::irange<Eigen::Index>(0, outMesh.nVertices()), activeAxis, polynomial);
    benchmark::DoNotOptimize(mat);
  }
}

BENCHMARK_CAPTURE(assembleA, 1, mapping::CompactPolynomialC6{1}, mapping::Polynomial::ON)->Name("Assemble RBF matrix A for C6 with polynomials");
BENCHMARK_CAPTURE(assembleA, 2, mapping::CompactPolynomialC6{1}, mapping::Polynomial::OFF)->Name("Assemble RBF matrix A for C6");
BENCHMARK_CAPTURE(assembleA, 3, mapping::ThinPlateSplines{}, mapping::Polynomial::ON)->Name("Assemble RBF matrix A for TPS with polynomials");
BENCHMARK_CAPTURE(assembleA, 4, mapping::ThinPlateSplines{}, mapping::Polynomial::OFF)->Name("Assemble RBF matrix A for TPS");

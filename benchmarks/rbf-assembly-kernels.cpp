#include <benchmark/benchmark.h>
#include <boost/range/irange.hpp>

#include "mapping/RadialBasisFctSolver.hpp"
#include "mapping/impl/BasisFunctions.hpp"

using namespace precice;

// Benchmark of our assemble routines

class SetupMeshes : public benchmark::Fixture {
public:
  void SetUp(const ::benchmark::State &state)
  {
    std::cout << "Setup" << std::endl;
    polynomial    = mapping::Polynomial::OFF;
    activeAxis    = {{true, true, true}};
    int dimension = 3;

    inMesh = std::make_unique<mesh::Mesh>("inMesh", dimension, -1);

    for (int i = 0; i < 10; ++i)
      for (int j = 0; j < 10; ++j)
        inMesh->createVertex(Eigen::Vector3d(0, i * 3.1, j * 1.4315));

    outMesh = std::make_unique<mesh::Mesh>("outMesh", dimension, -1);

    for (int i = 0; i < 10; ++i)
      for (int j = 0; j < 10; ++j)
        outMesh->createVertex(Eigen::Vector3d(1, i * 1.2111, j * 5.4315));

    std::cout << "end" << std::endl;
  }

  void TearDown(const ::benchmark::State &state)
  {
    inMesh.reset();
    outMesh.reset();
  }

  std::unique_ptr<mesh::Mesh> inMesh;
  std::unique_ptr<mesh::Mesh> outMesh;
  std::array<bool, 3>         activeAxis;
  mapping::Polynomial         polynomial{};
};

BENCHMARK_F(SetupMeshes, assembleA)
(benchmark::State &state)
{
  mapping::CompactPolynomialC6 f(1);
  Eigen::MatrixXd              mat;
  for (auto _ : state) {
    mat = mapping::buildMatrixA(f, *inMesh.get(), boost::irange<Eigen::Index>(0, inMesh->vertices().size()), *outMesh.get(), boost::irange<Eigen::Index>(0, outMesh->vertices().size()), activeAxis, polynomial);
  }
  benchmark::DoNotOptimize(mat);
}

#include <benchmark/benchmark.h>
#include <precice/SolverInterface.hpp>
#include <random>
#include <vector>

// A very simple write block vector data test. Works only if preCICE was built without MPI

class MyFixture : public benchmark::Fixture {
public:
  void SetUp(const ::benchmark::State &state)
  {
    p = std::make_unique<precice::SolverInterface>("A", "write-block-vector-data.xml", 0, 1);

    vids.resize(nv);
    std::vector<double> pos(nv * 3);
    std::generate(pos.begin(), pos.end(), []() -> double {
      static std::random_device               rd;
      static std::mt19937                     gen(rd());
      static std::uniform_real_distribution<> dis(0.0, 100.0);
      return dis(gen);
    });

    p->setMeshVertices(p->getMeshID("MeshA"), nv, pos.data(), vids.data());
    data.resize(nv * 3, 3.14159);
  }

  void TearDown(const ::benchmark::State &state)
  {
    p.reset();
  }

  std::unique_ptr<precice::SolverInterface> p;
  const int                                 nv = 1000000;
  std::vector<int>                          vids;
  std::vector<double>                       data;
};

BENCHMARK_DEFINE_F(MyFixture, writeblockScalarData)
(benchmark::State &state)
{
  const auto mid = p->getMeshID("MeshA");
  const auto did = p->getDataID("Vector", mid);

  for (auto _ : state) {
    p->writeBlockVectorData(did, nv, vids.data(), data.data());
  }
}

BENCHMARK_REGISTER_F(MyFixture, writeblockScalarData);

#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_CASE(PrimaryRankSockets)
{
  PRECICE_TEST("ParallelSolver"_on(3_ranks), "SerialSolver"_on(1_rank));

  std::string myMeshName;
  if (context.isNamed("ParallelSolver")) {
    myMeshName = "ParallelMesh";
  } else {
    myMeshName = "SerialMesh";
  }
  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);
  auto                     meshID      = myMeshName;
  double                   position[2] = {0, 0};
  interface.setMeshVertex(meshID, position);
  interface.initialize();
  interface.advance(1.0);
  interface.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel

#endif // PRECICE_NO_MPI

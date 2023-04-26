#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_CASE(LocalRBFPartitioning)
{
  PRECICE_TEST("SolverOne"_on(3_ranks), "SolverTwo"_on(1_rank));

  if (context.name == "SolverOne") {
    precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);
    auto                     meshName = "MeshOne";
    auto                     dataName = "Data2";

    int    vertexIDs[2];
    double xCoord       = context.rank * 0.4;
    double positions[4] = {xCoord, 0.0, xCoord + 0.2, 0.0};
    interface.setMeshVertices(meshName, 2, positions, vertexIDs);
    interface.initialize();
    double values[2];
    interface.advance(1.0);
    double preciceDt = interface.getMaxTimeStepSize();
    interface.readBlockScalarData(meshName, dataName, 2, vertexIDs, preciceDt, values);
    interface.finalize();
  } else {
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);
    auto                     meshName = "MeshTwo";
    int                      vertexIDs[6];
    double                   positions[12] = {0.0, 0.0, 0.2, 0.0, 0.4, 0.0, 0.6, 0.0, 0.8, 0.0, 1.0, 0.0};
    interface.setMeshVertices(meshName, 6, positions, vertexIDs);
    interface.initialize();
    auto   dataName  = "Data2";
    double values[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    interface.writeBlockScalarData(meshName, dataName, 6, vertexIDs, values);
    interface.advance(1.0);
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel

#endif // PRECICE_NO_MPI

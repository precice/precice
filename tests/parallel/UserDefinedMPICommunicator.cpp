#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

// Tests SolverInterface() with a user-defined MPI communicator.
BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_CASE(UserDefinedMPICommunicator)
{
  PRECICE_TEST("SolverOne"_on(3_ranks), "SolverTwo"_on(1_rank));

  /// @todo simplify once #1191 is merged
  if (context.isNamed("SolverOne")) {
    MPI_Comm                 myComm = precice::utils::Parallel::current()->comm;
    precice::SolverInterface interface(context.name, context.config(), context.rank, context.size, &myComm);
    auto                     meshName = "MeshOne";

    int    vertexIDs[2];
    double xCoord       = context.rank * 0.4;
    double positions[4] = {xCoord, 0.0, xCoord + 0.2, 0.0};
    interface.setMeshVertices(meshName, 2, positions, vertexIDs);
    interface.initialize();
    interface.finalize();
  } else {
    MPI_Comm                 myComm = precice::utils::Parallel::current()->comm;
    precice::SolverInterface interface(context.name, context.config(), context.rank, context.size, &myComm);
    auto                     meshName = "MeshTwo";
    int                      vertexIDs[6];
    double                   positions[12] = {0.0, 0.0, 0.2, 0.0, 0.4, 0.0, 0.6, 0.0, 0.8, 0.0, 1.0, 0.0};
    interface.setMeshVertices(meshName, 6, positions, vertexIDs);
    interface.initialize();
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel

#endif // PRECICE_NO_MPI

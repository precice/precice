#ifndef PRECICE_NO_MPI
#ifndef PRECICE_NO_PETSC

#include "testing/Testing.hpp"

#include <precice/Participant.hpp>
#include <vector>

// Tests Participant() with a user-defined MPI communicator.
// Since PETSc also uses MPI, we use petrbf mapping here.
BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_CASE(UserDefinedMPICommunicatorPetRBF)
{
  PRECICE_TEST("SolverOne"_on(3_ranks), "SolverTwo"_on(1_rank));

  if (context.isNamed("SolverOne")) {

    MPI_Comm             myComm = precice::utils::Parallel::current()->comm;
    precice::Participant interface(context.name, context.config(), context.rank, context.size, &myComm);

    auto   meshName = "MeshOne";
    int    vertexIDs[2];
    double xCoord       = context.rank * 0.4;
    double positions[4] = {xCoord, 0.0, xCoord + 0.2, 0.0};
    interface.setMeshVertices(meshName, positions, vertexIDs);
    interface.initialize();
    interface.finalize();
  } else {

    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    precice::Participant interface(context.name, context.config(), context.rank, context.size);

    auto   meshName = "MeshTwo";
    int    vertexIDs[6];
    double positions[12] = {0.0, 0.0, 0.2, 0.0, 0.4, 0.0, 0.6, 0.0, 0.8, 0.0, 1.0, 0.0};
    interface.setMeshVertices(meshName, positions, vertexIDs);
    interface.initialize();
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel

#endif // PRECICE_NO_PETSC
#endif // PRECICE_NO_MPI

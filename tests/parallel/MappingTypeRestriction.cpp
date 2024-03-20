#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_CASE(MappingTypeRestriction)

/**
 * @brief Reproduces bug that disallowed read-conservative mapping on serial participant
 *
 * SolverOne runs in parallel and uses a read-consistent mapping.
 * SolverTwo runs in serial and uses a read-conservative mapping.
 * A previous bug disallowed the read-conservative mapping during the configuration on the parallel participant.
 */
{
  PRECICE_TEST("SolverOne"_on(3_ranks), "SolverTwo"_on(1_rank));

  if (context.name == "SolverOne") {
    precice::Participant interface(context.name, context.config(), context.rank, context.size);
    auto                 meshName = "MeshOne";
    int                  vertexIDs[1];
    double               positions[2] = {0.0, 0.0};

    interface.setMeshVertices(meshName, positions, vertexIDs);
    interface.initialize();
    interface.finalize();
  } else {
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    precice::Participant interface(context.name, context.config(), context.rank, context.size);
    auto                 meshName = "MeshTwo";
    int                  vertexIDs[1];
    double               positions[2] = {0.0, 0.0};

    interface.setMeshVertices(meshName, positions, vertexIDs);
    interface.initialize();
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel

#endif // PRECICE_NO_MPI

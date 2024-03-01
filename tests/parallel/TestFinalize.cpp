#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/impl/ParticipantImpl.hpp>
#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_CASE(TestFinalize)
{
  PRECICE_TEST("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks));

  if (context.isNamed("SolverOne")) {
    precice::Participant participant(context.name, context.config(), context.rank, context.size);
    auto                 meshName = "MeshOne";
    double               xCoord   = 0.0 + context.rank;
    double               v[]      = {xCoord, 0, 0};
    participant.setMeshVertex(meshName, v);
    participant.initialize();
    BOOST_TEST(precice::testing::WhiteboxAccessor::impl(participant).mesh("MeshOne").nVertices() == 1);
    BOOST_TEST(precice::testing::WhiteboxAccessor::impl(participant).mesh("MeshTwo").nVertices() == 1);
    participant.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    precice::Participant participant(context.name, context.config(), context.rank, context.size);
    auto                 meshName = "MeshTwo";
    double               xCoord   = 0.0 + context.rank;
    double               v[]      = {xCoord, 0, 0};
    participant.setMeshVertex(meshName, v);
    participant.initialize();
    BOOST_TEST(precice::testing::WhiteboxAccessor::impl(participant).mesh("MeshTwo").nVertices() == 1);
    participant.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel

#endif // PRECICE_NO_MPI

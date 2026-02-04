#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Explicit)
BOOST_AUTO_TEST_SUITE(SerialCoupling)

/**
 * @brief ParticipantFirst leads to rounding errors in time management super quickly
 * Triggers an assertion in tw=3
 */
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(FirstParticipantTimeBug)
{
  PRECICE_TEST();

  Participant precice(context.name, context.config(), 0, 1);

  std::string meshName;
  if (context.isNamed("SolverOne")) {
    meshName = "SolverOne-Mesh";
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    meshName = "SolverTwo-Mesh";
  }

  double   v0[]     = {0, 0, 0};
  VertexID vertexID = precice.setMeshVertex(meshName, v0);
  precice.initialize();
  while (precice.isCouplingOngoing()) {
    if (context.isNamed("SolverOne")) {
      precice.advance(101.1);
    } else {
      // Hits assertion in tw=3
      precice.advance(precice.getMaxTimeStepSize());
    }
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Explicit
BOOST_AUTO_TEST_SUITE_END() // SerialCoupling

#endif // PRECICE_NO_MPI

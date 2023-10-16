#include <string>
#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Explicit)
BOOST_AUTO_TEST_SUITE(SerialCoupling)

/**
 * @brief Test to run a simple serial coupling.
 *
 * Performs many, small time steps to test for floating point inaccuracies.
 */
BOOST_AUTO_TEST_CASE(DoNothingWithSmallSteps)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  Participant precice(context.name, context.config(), 0, 1);

  double value = 1; // actual value does not matter, but we need to call readData to test for assertions in the bspline class.

  std::string meshName, writeDataName, readDataName;
  if (context.isNamed("SolverOne")) {
    meshName     = "MeshOne";
    readDataName = "DataTwo";
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    meshName     = "MeshTwo";
    readDataName = "DataOne";
  }

  double   v0[]     = {0, 0, 0};
  VertexID vertexID = precice.setMeshVertex(meshName, v0);
  precice.initialize();
  double dt = precice.getMaxTimeStepSize();

  while (precice.isCouplingOngoing()) {
    double dt = precice.getMaxTimeStepSize();
    precice.readData(meshName, readDataName, {&vertexID, 1}, dt, {&value, 1});
    precice.advance(dt);
  }

  BOOST_TEST(not precice.isCouplingOngoing());
  precice.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Explicit
BOOST_AUTO_TEST_SUITE_END() // SerialCoupling

#endif // PRECICE_NO_MPI

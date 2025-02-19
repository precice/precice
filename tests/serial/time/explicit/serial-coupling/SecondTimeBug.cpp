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
 * @brief Test based on #2160 reported by Claudio and Leonardo
 * Triggers an assertion in tw=7, t=60.06
 * Original issue used Dt=0.013 and crashed at t=64.012 after 4925 time windows
 */
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(SecondTimeBug)
{
  PRECICE_TEST();

  Participant precice(context.name, context.config(), 0, 1);

  std::array<double, 3> value{1, 1, 1}; // value doesn't matter

  std::string meshName, writeDataName, readDataName;
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
    double dt = precice.getMaxTimeStepSize();
    if (context.isNamed("SolverTwo")) {
      // Bug #2160, read fails in tw=7 t=60.06
      precice.readData(meshName, "Data-One", {&vertexID, 1}, dt, value);
    } else {
      precice.writeData(meshName, "Data-One", {&vertexID, 1}, value);
    }
    precice.advance(dt);
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Explicit
BOOST_AUTO_TEST_SUITE_END() // SerialCoupling

#endif // PRECICE_NO_MPI

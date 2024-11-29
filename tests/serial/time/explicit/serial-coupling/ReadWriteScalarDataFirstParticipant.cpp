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
 * @brief Test to run a simple serial coupling where the first participant prescribes the time window size.
 *
 * Ensures that time window sizes are passed correctly and that reading and writing is possible.
 */
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataFirstParticipant)
{
  PRECICE_TEST();

  Participant precice(context.name, context.config(), 0, 1);

  // SolverOne prescribes these, thus SolverTwo expect these (we use "first-participant" as dt method)
  std::vector<double> timestepSizes{1.0, 2.0, 3.0};

  // some dummy values, to check the actual values is not the point of this test,
  // but more to test whether reading / writing is possible at all
  double expectedDataValue = 2.5;
  double actualDataValue   = -1.0;

  std::string meshName, writeDataName, readDataName;
  if (context.isNamed("SolverOne")) {
    meshName      = "MeshOne";
    writeDataName = "DataOne";
    readDataName  = "DataTwo";
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    meshName      = "MeshTwo";
    writeDataName = "DataTwo";
    readDataName  = "DataOne";
  }

  double   v0[]     = {0, 0, 0};
  VertexID vertexID = precice.setMeshVertex(meshName, v0);
  BOOST_TEST(!precice.requiresInitialData());
  precice.initialize();

  int window = 0;

  for (auto prescribed_dt : timestepSizes) {
    BOOST_TEST(precice.isCouplingOngoing());
    double dt = precice.getMaxTimeStepSize();

    precice.readData(meshName, readDataName, {&vertexID, 1}, dt, {&actualDataValue, 1});
    if (window == 0 && context.isNamed("SolverOne")) {
      BOOST_TEST(actualDataValue == 0);
    } else {
      BOOST_TEST(actualDataValue == expectedDataValue);
    }
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&expectedDataValue, 1});

    if (context.isNamed("SolverOne")) {
      precice.advance(prescribed_dt);
    } else if (context.isNamed("SolverTwo")) {
      BOOST_TEST(dt == prescribed_dt);
      precice.advance(dt);
    }
    window++;
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

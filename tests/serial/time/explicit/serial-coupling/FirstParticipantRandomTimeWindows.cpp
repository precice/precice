#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <random>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Explicit)
BOOST_AUTO_TEST_SUITE(SerialCoupling)

PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(FirstParticipantRandomTimeWindows)
{
  PRECICE_TEST();

  Participant precice(context.name, context.config(), 0, 1);

  std::mt19937                        gen(101); // seed for reproducibility
  std::lognormal_distribution<double> dist(-2.0, 0.5);

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
    // SolverOne takes random steps
    double dt = context.isNamed("SolverOne") ? dist(gen) : precice.getMaxTimeStepSize();

    if (context.isNamed("SolverTwo")) {
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

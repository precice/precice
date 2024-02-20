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
 * @brief Test to ensure that correct max-time is reached, if time-window-size does not fit. See https://github.com/precice/precice/issues/1922 for context.
 *
 */
BOOST_AUTO_TEST_CASE(DoNonfittingWindows)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  Participant precice(context.name, context.config(), 0, 1);

  std::string meshName, writeDataName, readDataName;

  if (context.isNamed("SolverOne")) {
    meshName = "MeshOne";
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    meshName = "MeshTwo";
  }

  double v0[] = {0, 0, 0};
  precice.setMeshVertex(meshName, v0);
  BOOST_TEST(precice.requiresInitialData());
  precice.initialize();

  BOOST_TEST(precice.getMaxTimeStepSize() == 0.75);
  precice.advance(precice.getMaxTimeStepSize());
  BOOST_TEST(precice.isTimeWindowComplete());
  BOOST_TEST(precice.isCouplingOngoing());

  BOOST_TEST(precice.getMaxTimeStepSize() == 0.25);
  precice.advance(precice.getMaxTimeStepSize());
  BOOST_TEST(precice.isTimeWindowComplete());
  BOOST_TEST(!precice.isCouplingOngoing());
  BOOST_TEST(precice.getMaxTimeStepSize() == 0.0);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Explicit
BOOST_AUTO_TEST_SUITE_END() // SerialCoupling

#endif // PRECICE_NO_MPI

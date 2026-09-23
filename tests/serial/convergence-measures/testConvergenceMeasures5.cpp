#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(ConvergenceMeasures)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(testConvergenceMeasures5)
{
  PRECICE_TEST();
  using Eigen::Vector2d;

  std::string meshName = context.isNamed("SolverOne") ? "MeshOne" : "MeshTwo";

  precice::Participant p(context.name, context.config(), 0, 1);

  Vector2d vertex{0.0, 0.0};

  VertexID vertexID = p.setMeshVertex(meshName, vertex);

  p.initialize();
  int numberOfAdvanceCalls = 0;
  int numberOfIterations   = -1;
  int timestep             = 0;

  // Data3 converges on iteration 2 for time window 1, and iteration 3 for time window 2.
  // Data3 is written by SolverTwo on its provided MeshTwo but NOT exchanged.
  // This tests issue #2442: convergence measures using non-exchanged data.
  std::vector<double> writeValuesData3 = {1.0, 1.01, 2.0, 2.5, 2.51};

  // Expected iterations per time window based on the convergence of Data3
  std::vector<int> expectedIterations = {2, 3};

  while (p.isCouplingOngoing()) {
    BOOST_REQUIRE(numberOfAdvanceCalls < static_cast<int>(writeValuesData3.size()));
    if (p.requiresWritingCheckpoint()) {
      numberOfIterations = 0;
    }

    if (context.isNamed("SolverTwo")) {
      // Write Data3 (non-exchanged, used for convergence measure)
      p.writeData(meshName, "Data3", {&vertexID, 1}, {&writeValuesData3.at(numberOfAdvanceCalls), 1});
    }

    p.advance(p.getMaxTimeStepSize());
    ++numberOfAdvanceCalls;
    ++numberOfIterations;

    if (p.requiresReadingCheckpoint()) {
    } else { // converged
      BOOST_REQUIRE(timestep < static_cast<int>(expectedIterations.size()));
      BOOST_TEST(numberOfIterations == expectedIterations.at(timestep));
      ++timestep;
    }
  }
  p.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // ConvergenceMeasures

#endif // PRECICE_NO_MPI

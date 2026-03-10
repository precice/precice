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

  precice::Participant interface(context.name, context.config(), 0, 1);

  Vector2d vertex{0.0, 0.0};

  VertexID vertexID = interface.setMeshVertex(meshName, vertex);

  interface.initialize();
  int numberOfAdvanceCalls = 0;
  int numberOfIterations   = -1;
  int timestep             = 0;

  // Data3 converges on iteration 2 for time window 1, and iteration 3 for time window 2.
  // These values are written by SolverTwo on its provided MeshTwo.
  // Data3 is NOT exchanged, but used in the convergence measure.
  std::vector<double> writeValuesData3 = {1.0, 1.01, 2.0, 2.5, 2.51};

  // Expected iterations per time window based on the convergence of Data3
  std::vector<int> expectedIterations = {2, 3};

  while (interface.isCouplingOngoing()) {
    if (interface.requiresWritingCheckpoint()) {
      numberOfIterations = 0;
    }

    if (context.isNamed("SolverTwo")) {
      // Write Data3 (non-exchanged, used for convergence measure)
      interface.writeData(meshName, "Data3", {&vertexID, 1}, {&writeValuesData3.at(numberOfAdvanceCalls), 1});
    }

    interface.advance(1.0);
    ++numberOfAdvanceCalls;
    ++numberOfIterations;

    if (interface.requiresReadingCheckpoint()) {
    } else { // converged
      BOOST_TEST(numberOfIterations == expectedIterations.at(timestep));
      ++timestep;
    }
  }
  interface.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // ConvergenceMeasures

#endif // PRECICE_NO_MPI

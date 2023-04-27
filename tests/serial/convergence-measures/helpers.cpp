#ifndef PRECICE_NO_MPI

#include "helpers.hpp"
#include "testing/Testing.hpp"

#include "precice/SolverInterface.hpp"

void testConvergenceMeasures(const std::string configFile, TestContext const &context, std::vector<int> &expectedIterations)
{
  using Eigen::Vector2d;

  std::string meshName = context.isNamed("SolverOne") ? "MeshOne" : "MeshTwo";

  precice::SolverInterface interface(context.name, configFile, 0, 1);

  Vector2d vertex{0.0, 0.0};

  std::vector<double> writeValues = {1.0, 1.01, 2.0, 2.5, 2.8, 2.81};

  VertexID vertexID = interface.setMeshVertex(meshName, vertex.data());

  interface.initialize();
  int numberOfAdvanceCalls = 0;
  int numberOfIterations   = -1;
  int timestep             = 0;

  while (interface.isCouplingOngoing()) {
    if (interface.requiresWritingCheckpoint()) {
      numberOfIterations = 0;
    }

    if (context.isNamed("SolverTwo")) {
      auto dataName = "Data2";
      interface.writeData(meshName, dataName, {&vertexID, 1}, {&writeValues.at(numberOfAdvanceCalls), 1});
    }

    interface.advance(1.0);
    ++numberOfAdvanceCalls;
    ++numberOfIterations;

    if (interface.requiresReadingCheckpoint()) {
    } else { //converged
      BOOST_TEST(numberOfIterations == expectedIterations.at(timestep));
      ++timestep;
    }
  }
  interface.finalize();
}

#endif

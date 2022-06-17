#ifndef PRECICE_NO_MPI

#include "helpers.hpp"
#include "testing/Testing.hpp"

#include "precice/SolverInterface.hpp"

void testConvergenceMeasures(const std::string configFile, TestContext const &context, std::vector<int> &expectedIterations)
{
  using Eigen::Vector2d;
  using namespace precice::constants;

  std::string meshName = context.isNamed("SolverOne") ? "MeshOne" : "MeshTwo";

  precice::SolverInterface interface(context.name, configFile, 0, 1);
  const int                meshID = interface.getMeshID(meshName);

  Vector2d vertex{0.0, 0.0};

  std::vector<double> writeValues = {1.0, 1.01, 2.0, 2.5, 2.8, 2.81};

  VertexID vertexID = interface.setMeshVertex(meshID, vertex.data());

  interface.initialize();
  interface.initializeData();
  int numberOfAdvanceCalls = 0;
  int numberOfIterations   = -1;
  int timestep             = 0;

  while (interface.isCouplingOngoing()) {
    if (interface.isActionRequired(actionWriteIterationCheckpoint())) {
      numberOfIterations = 0;
      interface.markActionFulfilled(actionWriteIterationCheckpoint());
    }

    if (context.isNamed("SolverTwo")) {
      precice::DataID dataID = interface.getDataID("Data2", meshID);
      interface.writeScalarData(dataID, vertexID, writeValues.at(numberOfAdvanceCalls));
    }

    interface.advance(1.0);
    ++numberOfAdvanceCalls;
    ++numberOfIterations;

    if (interface.isActionRequired(actionReadIterationCheckpoint())) {
      interface.markActionFulfilled(actionReadIterationCheckpoint());
    } else { //converged
      BOOST_TEST(numberOfIterations == expectedIterations.at(timestep));
      ++timestep;
    }
  }
  interface.finalize();
}

#endif
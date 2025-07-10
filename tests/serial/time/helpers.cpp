#ifndef PRECICE_NO_MPI

#include "helpers.hpp"
#include "testing/Testing.hpp"

#include <fstream>
#include "precice/precice.hpp"

void doManySteps(TestContext const &context)
{
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

  while (precice.isCouplingOngoing()) {
    double dt = precice.getMaxTimeStepSize();
    precice.requiresWritingCheckpoint(); // to avoid error for implicit coupling schemes
    precice.readData(meshName, readDataName, {&vertexID, 1}, 0, {&value, 1});
    precice.readData(meshName, readDataName, {&vertexID, 1}, dt, {&value, 1});
    precice.advance(dt);
    precice.requiresReadingCheckpoint(); // to avoid error for implicit coupling schemes
  }

  BOOST_TEST(not precice.isCouplingOngoing());
  precice.finalize();
}

#endif

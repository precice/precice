#include <string>
#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Implicit)
BOOST_AUTO_TEST_SUITE(SerialCoupling)

PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(DirtyTimestep)
{
  PRECICE_TEST();

  std::string meshName, dataWriteName, dataReadName;
  if (context.isNamed("SolverOne")) {
    meshName      = "MeshOne";
    dataWriteName = "DataOne";
    dataReadName  = "DataTwo";
  } else {
    meshName      = "MeshTwo";
    dataWriteName = "DataTwo";
    dataReadName  = "DataOne";
  }

  Participant participant(context.name, context.config(), 0, 1);

  const int dimensions       = participant.getMeshDimensions(meshName);
  const int numberOfVertices = 3;

  std::vector<double> vertices(numberOfVertices * dimensions);
  std::vector<int>    vertexIDs(numberOfVertices);
  for (int i = 0; i < numberOfVertices; i++)
    for (int j = 0; j < dimensions; j++)
      vertices.at(j + dimensions * i) = i;
  participant.setMeshVertices(meshName, vertices, vertexIDs);

  // Data is scalar in the config -> one component per vertex.
  const int           dataSize = numberOfVertices;
  std::vector<double> readData(dataSize, 0.0);
  std::vector<double> writeData(dataSize, 0.0);

  participant.initialize();

  double t           = 0.0;
  double tCheckpoint = -1;

  while (participant.isCouplingOngoing()) {
    if (participant.requiresWritingCheckpoint()) {
      tCheckpoint = t;
    }

    const double precice_dt = participant.getMaxTimeStepSize();

    // naive way of the solver integrating its own clock and reporting the lossy difference
    const double tOld = t;
    t += precice_dt;
    const double dt = t - tOld;

    // read both 0.0 and dt for testing purposes
    participant.readData(meshName, dataReadName, vertexIDs, 0.0, readData);
    participant.readData(meshName, dataReadName, vertexIDs, dt, readData);

    for (int i = 0; i < dataSize; i++) {
      writeData.at(i) = readData.at(i) + 1.0;
    }
    participant.writeData(meshName, dataWriteName, vertexIDs, writeData);

    participant.advance(dt);

    if (participant.requiresReadingCheckpoint()) {
      t = tCheckpoint;
    }
  }

  participant.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Implicit
BOOST_AUTO_TEST_SUITE_END() // SerialCoupling

#endif // PRECICE_NO_MPI

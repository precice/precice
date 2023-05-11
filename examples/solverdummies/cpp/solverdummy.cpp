#include <iostream>
#include <sstream>
#include <vector>
#include "precice/Participant.hpp"

int main(int argc, char **argv)
{
  int commRank = 0;
  int commSize = 1;

  using namespace precice;

  if (argc != 3) {
    std::cout << "The solverdummy was called with an incorrect number of arguments. Usage: ./solverdummy configFile solverName\n\n";
    std::cout << "Parameter description\n";
    std::cout << "  configurationFile: Path and filename of preCICE configuration\n";
    std::cout << "  solverName:        SolverDummy participant name in preCICE configuration\n";
    return EXIT_FAILURE;
  }

  std::string configFileName(argv[1]);
  std::string solverName(argv[2]);
  std::string meshName;
  std::string dataWriteName;
  std::string dataReadName;

  std::cout << "DUMMY: Running solver dummy with preCICE config file \"" << configFileName << "\" and participant name \"" << solverName << "\".\n";

  Participant interface(solverName, configFileName, commRank, commSize);

  if (solverName == "SolverOne") {
    dataWriteName = "Data-One";
    dataReadName  = "Data-Two";
    meshName      = "SolverOne-Mesh";
  }
  if (solverName == "SolverTwo") {
    dataReadName  = "Data-One";
    dataWriteName = "Data-Two";
    meshName      = "SolverTwo-Mesh";
  }

  int dimensions       = interface.getMeshDimensions(meshName);
  int numberOfVertices = 3;

  std::vector<double> readData(numberOfVertices * dimensions);
  std::vector<double> writeData(numberOfVertices * dimensions);
  std::vector<double> vertices(numberOfVertices * dimensions);
  std::vector<int>    vertexIDs(numberOfVertices);

  for (int i = 0; i < numberOfVertices; i++) {
    for (int j = 0; j < dimensions; j++) {
      vertices.at(j + dimensions * i)  = i;
      readData.at(j + dimensions * i)  = i;
      writeData.at(j + dimensions * i) = i;
    }
  }

  interface.setMeshVertices(meshName, vertices, vertexIDs);

  if (interface.requiresInitialData()) {
    std::cout << "DUMMY: Writing initial data\n";
  }

  interface.initialize();

  while (interface.isCouplingOngoing()) {

    if (interface.requiresWritingCheckpoint()) {
      std::cout << "DUMMY: Writing iteration checkpoint\n";
    }

    double dt = interface.getMaxTimeStepSize();
    interface.readData(meshName, dataReadName, vertexIDs, dt, readData);

    for (int i = 0; i < numberOfVertices * dimensions; i++) {
      writeData.at(i) = readData.at(i) + 1;
    }

    interface.writeData(meshName, dataWriteName, vertexIDs, writeData);

    interface.advance(dt);

    if (interface.requiresReadingCheckpoint()) {
      std::cout << "DUMMY: Reading iteration checkpoint\n";
    } else {
      std::cout << "DUMMY: Advancing in time\n";
    }
  }

  interface.finalize();
  std::cout << "DUMMY: Closing C++ solver dummy...\n";

  return 0;
}

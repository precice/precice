// To compile use:
// mpic++ -I$PRECICE_ROOT/src main.cpp -lprecice -o solverdummy

#include <iostream>
#include <sstream>
#include "precice/SolverInterface.hpp"

int main(int argc, char **argv)
{
  int commRank = 0;
  int commSize = 1;

  using namespace precice;
  using namespace precice::constants;

  if (argc != 4) {
    std::cout << "Usage: ./solverdummy configFile solverName meshName\n\n";
    std::cout << "Parameter description\n";
    std::cout << "  configurationFile: Path and filename of preCICE configuration\n";
    std::cout << "  solverName:        SolverDummy participant name in preCICE configuration\n";
    std::cout << "  meshName:          Mesh in preCICE configuration that carries read and write data\n";
    return 1;
  }

  std::string configFileName(argv[1]);
  std::string solverName(argv[2]);
  std::string meshName(argv[3]);

  std::cout << "DUMMY: Running solver dummy with preCICE config file \"" << configFileName << "\", participant name \"" << solverName << "\", and mesh name \"" << meshName << "\".\n";

  SolverInterface interface(solverName, configFileName, commRank, commSize);

  int meshID     = interface.getMeshID(meshName);
  int dimensions = interface.getDimensions();
  std::string dataWriteName; 
  std::string dataReadName;
  int N = 3;            // Number of vertices

  if (solverName == "SolverOne"){
    dataWriteName="Forces";
    dataReadName="Velocities";
  }
  if (solverName == "SolverTwo"){
    dataReadName="Forces";
    dataWriteName="Velocities";
  }
  const int readDataID = interface.getDataID(dataReadName,meshID);
  const int writeDataID = interface.getDataID(dataWriteName,meshID);

  double vertexPos[N * dimensions];
  double readData[N * dimensions];
  double writeData[N * dimensions];
  int vertexIDs[N];

  for (int i = 0; i < N; i++){
    for (int j = 0; j < dimensions; j++) {
      vertexPos[j + N*i]=i;
      readData[j + N*i] = i;
      writeData[j + N*i] = i;
    }
    vertexIDs[i] = i;
    std::cout << "DUMMY: Vertex: " << vertexPos[i*N + 0] << " ; " << vertexPos[i*N + 1] << " ; " << vertexPos[i*N + 2] << "\n";
    std::cout << "DUMMY: Read Data: " << readData[i*N + 0] << " ; " << readData[i*N + 1] << " ; " << readData[i*N + 2] << "\n";
    std::cout << "DUMMY: Write Data: " << writeData[i*N + 0] << " ; " << writeData[i*N + 1] << " ; " << writeData[i*N + 2] << "\n";
  }
  interface.setMeshVertices(meshID,N,vertexPos,vertexIDs);

  double dt = interface.initialize();

  while (interface.isCouplingOngoing()) {

    if (interface.isActionRequired(actionWriteIterationCheckpoint())) {
      interface.writeBlockVectorData(writeDataID,N,vertexIDs,writeData);
      for (int i = 0; i < N; i++){
        std::cout << "DUMMY: Write Data: " << writeData[i*N + 0] << " ; " << writeData[i*N + 1] << " ; " << writeData[i*N + 2] << "\n";
      }
      std::cout << "DUMMY: Writing iteration checkpoint\n";
      interface.markActionFulfilled(actionWriteIterationCheckpoint());
    }


    dt = interface.advance(dt);

    if (interface.isActionRequired(actionReadIterationCheckpoint())) {
      interface.readBlockVectorData(readDataID,N,vertexIDs,readData);
      for (int i = 0; i < N; i++){
        std::cout << "DUMMY: Read Data: " << readData[i*N + 0] << " ; " << readData[i*N + 1] << " ; " << readData[i*N + 2] << "\n";
      }
      std::cout << "DUMMY: Reading iteration checkpoint\n";
      interface.markActionFulfilled(actionReadIterationCheckpoint());
    } else {
      std::cout << "DUMMY: Advancing in time\n";
    }

    for (int i = 0; i < N*dimensions; i++){
      writeData[i] = readData[i] + 1;
    }

  }

  interface.finalize();
  std::cout << "DUMMY: Closing C++ solver dummy...\n";

  return 0;
}

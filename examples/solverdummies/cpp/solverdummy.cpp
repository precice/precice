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
  int N = 6;            // Number of vertices

  if (solverName == "SolverOne"){
    dataWriteName="dataOne";
    dataReadName="dataTwo";
  }
  if (solverName == "SolverTwo"){
    dataReadName="dataOne";
    dataWriteName="dataTwo";
  }
  const int readDataID = interface.getDataID(dataReadName,meshID);
  const int writeDataID = interface.getDataID(dataWriteName,meshID);

  std::vector<double> readData(N*dimensions);
  std::vector<double> writeData(N*dimensions);
  std::vector<double> vertex(dimensions);
  std::vector<int> vertexIDs(N);

  for (int i = 0; i < N; i++){
    if (i == 1){
      vertex[0] = 1.0;
      vertex[1] = 0.0;
      vertex[2] = 0.0;
    } else if (i == 2) {
      vertex[0] = 0.0;
      vertex[1] = 0.0;
      vertex[2] = 0.0;
    } else if (i == 3) {
      vertex[0] = 2.0;
      vertex[1] = 2.0;
      vertex[2] = 0.0;
    } else if (i == 4) {
      vertex[0] = 0.0;
      vertex[1] = 1.0;
      vertex[2] = 0.0;
    }
    else if (i == 0) {
      vertex[0] = 0.0;
      vertex[1] = 1.0;
      vertex[2] = 1.0;
    }
    else if (i == 5) {
      vertex[0] = 0.0;
      vertex[1] = 1.0;
      vertex[2] = 2.0;
    }
    /*
    if (i == 0){
      vertex[0] = 0.5;
      vertex[1] = 0.4;
      vertex[2] = 0.0;
    } else if (i == 1) {
      vertex[0] = 1.0;
      vertex[1] = 1.0;
      vertex[2] = 0.0;
    } else if (i == 2) {
      vertex[0] = 1.0;
      vertex[1] = 0.0;
      vertex[2] = 0.0;
    } else {
      vertex[0] = 0.0;
      vertex[1] = 1.0;
      vertex[2] = 0.0;
    }
    */
    for (int j = 0; j < dimensions; j++) {
      //vertex[j] = i;
      readData[j + i*dimensions] = i;
      writeData[j + i*dimensions] = i;
      
    }
    vertexIDs[i] = interface.setMeshVertex(meshID, vertex.data());
  }

  interface.setMeshEdge(meshID, 1, 2);
  interface.setMeshEdge(meshID, 4, 3);
  interface.setMeshEdge(meshID, 4, 2);
  interface.setMeshEdge(meshID, 1, 3);

  interface.setMeshQuad(meshID, 0, 1, 2, 3);
  //interface.setMeshQuadWithEdges(meshID, 0, 1, 2, 3);

  double dt = interface.initialize();

  while (interface.isCouplingOngoing()) {

    interface.readBlockVectorData(readDataID,N,vertexIDs.data(),readData.data());
    std::cout << "DUMMY: Reading iteration checkpoint\n";
    interface.markActionFulfilled(actionReadIterationCheckpoint());

    for (int i = 0; i < N*dimensions; i++){
      writeData[i] = readData[i] + 1;
    }

    interface.writeBlockVectorData(writeDataID,N,vertexIDs.data(),writeData.data());
    std::cout << "DUMMY: Writing iteration checkpoint\n";
    interface.markActionFulfilled(actionWriteIterationCheckpoint());

    dt = interface.advance(dt);

    std::cout << "DUMMY: Advancing in time\n";

  }

  interface.finalize();
  std::cout << "DUMMY: Closing C++ solver dummy...\n";

  return 0;
}
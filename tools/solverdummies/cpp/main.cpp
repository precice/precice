// To compile use:
// mpic++ -I$PRECICE_ROOT/src main.cpp -lprecice -o solverdummy

#include <iostream>
#include <sstream>
#include "precice/SolverInterface.hpp"
#include "precice/Constants.hpp"
#include <mpi.h>

int main (int argc, char **argv)
{
  std::cout << "Starting SolverDummy..." << std::endl;
  MPI_Init(&argc, &argv);
  int commRank = 0;
  int commSize = 1;

  using namespace precice;
  using namespace precice::constants;

  if (argc != 4){
    std::cout << "Usage: ./solverdummy configFile solverName meshName" << std::endl;
    std::cout << std::endl;
    std::cout << "Parameter description" << std::endl;
    std::cout << "  configurationFile: Path and filename of preCICE configuration" << std::endl;
    std::cout << "  solverName:        SolverDummy participant name in preCICE configuration" << std::endl;
    std::cout << "  meshName:          Mesh in preCICE configuration that carries read and write data" << std::endl;
    return 1;
}

  std::string configFileName(argv[1]);
  std::string solverName(argv[2]);
  std::string meshName(argv[3]);
  int N = 1;

  SolverInterface interface(solverName, commRank, commSize);
  interface.configure(configFileName);

  int meshID = interface.getMeshID(meshName);
  int dimensions = interface.getDimensions();
  std::vector<std::vector<double>> vertices(N, std::vector<double>(dimensions, 0)); 
  std::vector<int> dataIndices(N,0);

  for(int i=0; i<N; i++){
    dataIndices[i] = interface.setMeshVertex(meshID, vertices[i].data());
  }

  double dt = interface.initialize();

  while (interface.isCouplingOngoing()){

    if (interface.isActionRequired(actionWriteIterationCheckpoint())){
      std::cout << "DUMMY: Writing iteration checkpoint" << std::endl;
      interface.fulfilledAction(actionWriteIterationCheckpoint());
    }

    dt = interface.advance(dt);

    if (interface.isActionRequired(actionReadIterationCheckpoint())){
      std::cout << "DUMMY: Writing iteration checkpoint" << std::endl;
      interface.fulfilledAction(actionReadIterationCheckpoint());
    }
    else{
      std::cout << "DUMMY: Advancing in time" << std::endl;
    }
  }

  interface.finalize();
  std::cout << "DUMMY: Closing python solver dummy..." << std::endl;

  return 0;
}

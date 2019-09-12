/*
 * To compile, run the following command in MATLAB:
 * mex solverdummy.cpp -lprecice -output solverdummy
 * The result will be a function callable in MATLAB.
 *
 * To run, run the following command in MATLAB:
 * solverdummy configFile solverName meshName
 * OR, equivalently,
 * solverdummy('configFile','solverName','meshName')
 * since MATLAB syntax works like this.
*/

#include "precice/SolverInterface.hpp"
#include "precice/Constants.hpp"
#include <matrix.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mexPrintf("Starting SolverDummy...\n");
  int commRank = 0;
  int commSize = 1;

  using namespace precice;
  using namespace precice::constants;

  if (nrhs != 3){
    
    mexPrintf("Usage: ./solverdummy configFile solverName meshName\n\n");
    mexPrintf("Parameter description:\n");
    mexPrintf(" configFile: Path and filename of preCICE configuration\n");
    mexPrintf(" solverName: SolverDummy participant name in preCICE configuration\n");
    mexPrintf(" meshName:   Mesh in preCICE configuration that carries read and write data\n");
    return;
  }
  
  std::string configFileName(mxArrayToString(prhs[0]));
  std::string solverName(mxArrayToString(prhs[1]));
  std::string meshName(mxArrayToString(prhs[2]));
  int N = 1;
  SolverInterface interface(solverName, commRank, commSize);
  interface.configure(configFileName);

  int meshID = interface.getMeshID(meshName);
  int dimensions = interface.getDimensions();
  mexPrintf(meshName.c_str());
  std::vector<std::vector<double>> vertices(N, std::vector<double>(dimensions, 0)); 
  std::vector<int> dataIndices(N,0);

  for(int i=0; i<N; i++){
    dataIndices[i] = interface.setMeshVertex(meshID, vertices[i].data());
  }
  double dt = interface.initialize();

  while (interface.isCouplingOngoing()){

    if (interface.isActionRequired(actionWriteIterationCheckpoint())){
      mexPrintf("DUMMY: Writing iteration checkpoint\n");
      interface.fulfilledAction(actionWriteIterationCheckpoint());
    }

    dt = interface.advance(dt);

    if (interface.isActionRequired(actionReadIterationCheckpoint())){
      mexPrintf("DUMMY: Writing iteration checkpoint\n");
      interface.fulfilledAction(actionReadIterationCheckpoint());
    }
    else{
      mexPrintf("DUMMY: Advancing in time\n");
    }
  }
  interface.finalize();
  mexPrintf("DUMMY: Closing C++ solver dummy...\n");
}

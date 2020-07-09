#include "precice/SolverInterfaceC.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
  double  dt                 = 0.0;
  int     solverProcessIndex = 0;
  int     solverProcessSize  = 1;
  int     dimensions         = -1;
  double *vertex;
  double *readData;
  double *writeData;
  int     meshID   = -1;
  int     dataID   = -1;
  int    *vertexIDs;
  int N = 3;
  int writeDataID = -1;
  int readDataID = -1;
  

  if (argc != 4) {
    printf("Usage: ./solverdummy configFile solverName meshName\n\n");
    printf("Parameter description\n");
    printf("  configurationFile: Path and filename of preCICE configuration\n");
    printf("  solverName:        SolverDummy participant name in preCICE configuration\n");
    printf("  meshName:          Mesh in preCICE configuration that carries read and write data\n");
    return 1;
  }

  const char *configFileName  = argv[1];
  const char *participantName = argv[2];
  const char *meshName        = argv[3];
  
  printf("DUMMY: Running solver dummy with preCICE config file \"%s\", participant name \"%s\", and mesh name \"%s\".\n",
         configFileName, participantName, meshName);

  const char *writeItCheckp = precicec_actionWriteIterationCheckpoint();
  const char *readItCheckp  = precicec_actionReadIterationCheckpoint();

  precicec_createSolverInterface(participantName, configFileName, solverProcessIndex, solverProcessSize);

  meshID = precicec_getMeshID(meshName);

  if (strcmp(participantName,"SolverOne") == 0){
    writeDataID = precicec_getDataID("dataOne",meshID);
    readDataID = precicec_getDataID("dataTwo",meshID);
  }
  if (strcmp(participantName,"SolverTwo") == 0){
    writeDataID = precicec_getDataID("dataTwo",meshID);
    readDataID = precicec_getDataID("dataOne",meshID);
  }

  dimensions = precicec_getDimensions();
  vertex     = malloc(N * dimensions * sizeof(double));
  readData   = malloc(N * dimensions * sizeof(double));
  writeData   = malloc(N * dimensions * sizeof(double));
  vertexIDs = malloc(N * sizeof(int));

  for (int i = 0; i < N; i++){
    for (int j = 0; j < dimensions; j++) {
      vertex[j + N*i] = i;
      readData[j + N*i] = i;
      writeData[j + N*i] = i;
    }
    vertexIDs[i] = i;
  }

  precicec_setMeshVertices(meshID, N, vertex, vertexIDs);

  free(vertex);

  dt = precicec_initialize();

  while (precicec_isCouplingOngoing()) {

    precicec_readBlockVectorData(readDataID, N, vertexIDs, readData);
    printf("DUMMY: Reading iteration checkpoint \n");
    precicec_markActionFulfilled(readItCheckp);

    for (int i = 0; i < N*dimensions; i++){
      writeData[i] = readData[i] + 1;
    }

    precicec_writeBlockVectorData(writeDataID, N, vertexIDs, writeData);
    printf("DUMMY: Writing iteration checkpoint \n");
    precicec_markActionFulfilled(writeItCheckp);

    printf("DUMMY: Advancing in time \n");
    dt = precicec_advance(dt);
  
  }

  precicec_finalize();
  free(writeData);
  free(readData);
  free(vertexIDs);
  printf("DUMMY: Closing C solver dummy... \n");

  return 0;
}

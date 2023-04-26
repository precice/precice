#include "precice/SolverInterfaceC.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
  double      dt                 = 0.0;
  int         solverProcessIndex = 0;
  int         solverProcessSize  = 1;
  int         dimensions         = -1;
  double *    vertices;
  double *    readData;
  double *    writeData;
  int *       vertexIDs;
  int         numberOfVertices = 3;
  int         writeDataID      = -1;
  int         readDataID       = -1;
  const char *meshName;
  const char *writeDataName;
  const char *readDataName;

  if (argc != 3) {
    printf("The solverdummy was called with an incorrect number of arguments. Usage: ./solverdummy configFile solverName\n\n");
    printf("Parameter description\n");
    printf("  configurationFile: Path and filename of preCICE configuration\n");
    printf("  solverName:        SolverDummy participant name in preCICE configuration\n");
    return 1;
  }

  const char *configFileName  = argv[1];
  const char *participantName = argv[2];

  printf("DUMMY: Running solver dummy with preCICE config file \"%s\" and participant name \"%s\".\n",
         configFileName, participantName);

  precicec_createSolverInterface(participantName, configFileName, solverProcessIndex, solverProcessSize);

  if (strcmp(participantName, "SolverOne") == 0) {
    writeDataName = "Data-One";
    readDataName  = "Data-Two";
    meshName      = "SolverOne-Mesh";
  }
  if (strcmp(participantName, "SolverTwo") == 0) {
    writeDataName = "Data-Two";
    readDataName  = "Data-One";
    meshName      = "SolverTwo-Mesh";
  }

  dimensions = precicec_getDimensions();
  vertices   = malloc(numberOfVertices * dimensions * sizeof(double));
  readData   = malloc(numberOfVertices * dimensions * sizeof(double));
  writeData  = malloc(numberOfVertices * dimensions * sizeof(double));
  vertexIDs  = malloc(numberOfVertices * sizeof(int));

  for (int i = 0; i < numberOfVertices; i++) {
    for (int j = 0; j < dimensions; j++) {
      vertices[j + dimensions * i]  = i;
      readData[j + dimensions * i]  = i;
      writeData[j + dimensions * i] = i;
    }
  }

  precicec_setMeshVertices(meshName, numberOfVertices, vertices, vertexIDs);

  free(vertices);

  if (precicec_requiresInitialData()) {
    printf("DUMMY: Writing initial data\n");
  }

  precicec_initialize();

  while (precicec_isCouplingOngoing()) {

    if (precicec_requiresWritingCheckpoint()) {
      printf("DUMMY: Writing iteration checkpoint \n");
    }

    dt = precicec_getMaxTimeStepSize();
    precicec_readBlockVectorData(meshName, readDataName, numberOfVertices, vertexIDs, dt, readData);

    for (int i = 0; i < numberOfVertices * dimensions; i++) {
      writeData[i] = readData[i] + 1;
    }

    precicec_writeBlockVectorData(meshName, writeDataName, numberOfVertices, vertexIDs, writeData);

    precicec_advance(dt);

    if (precicec_requiresReadingCheckpoint()) {
      printf("DUMMY: Reading iteration checkpoint \n");
    } else {
      printf("DUMMY: Advancing in time \n");
    }
  }

  precicec_finalize();
  free(writeData);
  free(readData);
  free(vertexIDs);
  printf("DUMMY: Closing C solver dummy... \n");

  return 0;
}

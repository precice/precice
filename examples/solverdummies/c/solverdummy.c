#include "precice/SolverInterfaceC.h"

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
  double  dt                 = 0.0;
  int     solverProcessIndex = 0;
  int     solverProcessSize  = 1;
  int     dimensions         = -1;
  double *vertex;
  int     meshID   = -1;
  int     dataID   = -1;
  int     vertexID = -1;

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

  dimensions = precicec_getDimensions();
  vertex     = malloc(dimensions * sizeof(double));

  for (int i = 0; i < dimensions; i++) {
    vertex[i] = 0.0;
  }

  vertexID = precicec_setMeshVertex(meshID, vertex);
  free(vertex);

  dt = precicec_initialize();

  while (precicec_isCouplingOngoing()) {

    if (precicec_isActionRequired(writeItCheckp)) {
      printf("DUMMY: Writing iteration checkpoint");
      precicec_markActionFulfilled(writeItCheckp);
    }

    dt = precicec_advance(dt);

    if (precicec_isActionRequired(readItCheckp)) {
      precicec_markActionFulfilled(readItCheckp);
      printf("DUMMY: Reading iteration checkpoint");
    } else {
      printf("DUMMY: Advancing in time");
    }
  }

  precicec_finalize();
  printf("DUMMY: Closing C solver dummy...");

  return 0;
}

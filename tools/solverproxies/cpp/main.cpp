#include <iostream>
#include <sstream>
#include <cstdlib>
//#include <ctime>
#include "../../../src/precice/SolverInterface.hpp"
#include "mpi.h"

/**
 * @brief For printing to the command line.
 *
 * @param message  An input stream such as: "Blabla = " << variable << ", ..."
 */
#define PRINT(message) \
  { \
    std::ostringstream conv; \
    conv << "(" << comm_rank << "/" << comm_size << ") "; \
    conv << message; \
    std::cout << conv.str() << std::endl; \
  }

void printData (const std::vector<double>& data)
{
  std::cout << "Received data = " << data[0];
  for (size_t i=1; i < data.size(); i++){
    std::cout << ", " << data[i];
  }
  std::cout << std::endl;
}


int main (int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  int comm_rank = -1;
  int comm_size = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  PRINT("Starting solver dummy...");
  using namespace precice;
  using namespace precice::constants;

  if (argc == 1 || (argc != 3 && argc != 7)){
    PRINT("Usage: ./solverproxy configFile proxyName [meshName readDataName writeDataName computationTimeInSeconds]");
    PRINT("");
    PRINT("Parameters in [] are optional, but are needed together.");
    PRINT("If optional parameters are given, the proxy writes and reads data.");
    PRINT("");
    PRINT("Parameter description");
    PRINT("  configurationFile: Path and filename of preCICE configuration");
    PRINT("  proxyName:         Proxy participant name in preCICE configuration");
    PRINT("  meshName:          Mesh in preCICE configuration that carries read and write data");
    PRINT("  readDataName:      Data in preCICE config. that is read by this proxy");
    PRINT("  writeDataName:     Data in preCICE config. that is written by this proxy");
    PRINT("  computationTimeInSeconds: Time waited by proxy in every computation cycle");
    return 1;
  }
  std::string configFileName(argv[1]);
  std::string proxyName(argv[2]);
  bool readWriteData = false;
  std::string meshName;
  std::string readDataName;
  std::string writeDataName;
  double computationTime = 0.0;
  if (argc > 3){
    readWriteData = true;
    meshName = argv[3];
    readDataName = argv[4];
    writeDataName = argv[5];
    computationTime = atof(argv[6]);
  }

  SolverInterface interface(proxyName, comm_rank, comm_size);
  interface.configure(configFileName);

  double computedTime = 0.0;
  int computedTimeSteps = 0;

  int meshID = -1;
  int readDataID = -1;
  int writeDataID = -1;
  int dimensions = -1;
  if (readWriteData){
    meshID = interface.getMeshID(meshName);
    readDataID = interface.getDataID(readDataName);
    writeDataID = interface.getDataID(writeDataName);
    dimensions = interface.getDimensions();
  }

  if (interface.isActionRequired(actionReadSimulationCheckpoint())){
    interface.fulfilledAction(actionReadSimulationCheckpoint());
  }

  double dt = interface.initialize();

  int dataSize = -1;
  double* data = NULL;
  int* dataIndices = NULL;
  if (readWriteData){
    dataSize = interface.getMeshVertexSize(meshID);
    data = new double[dataSize*dimensions];
    dataIndices = new int[dataSize];
    for (int i=0; i < dataSize; i++){
      dataIndices[i] = i;
    }
  }

  double mpi_start_time = MPI_Wtime();
  double mpi_compute_time = 0.0;

  if (readWriteData && interface.isReadDataAvailable()){
    interface.readBlockVectorData(readDataID, dataSize, dataIndices, data);
  }

  while (interface.isCouplingOngoing()){
    // When an implicit coupling scheme is used, checkpointing is required
    if (interface.isActionRequired(actionWriteIterationCheckpoint())){
      PRINT(actionWriteIterationCheckpoint());
      interface.fulfilledAction(actionWriteIterationCheckpoint());
    }

    // Wait for computedTime milliseconds
    PRINT("Computing for " << computationTime << " seconds ...");
    double mpi_compute_start = MPI_Wtime();
    while (MPI_Wtime() - mpi_compute_start < computationTime) {}
    double mpi_compute_end = MPI_Wtime();
    mpi_compute_time += mpi_compute_end - mpi_compute_start;
    PRINT("...done");

    computedTime += dt;
    computedTimeSteps++;

    if (readWriteData){
      interface.writeBlockVectorData(writeDataID, dataSize, dataIndices, data);
    }

    dt = interface.advance(dt);

    if (readWriteData && interface.isReadDataAvailable()){
      interface.readBlockVectorData(readDataID, dataSize, dataIndices, data);
    }

    if (interface.isActionRequired(actionWriteSimulationCheckpoint())){
      interface.fulfilledAction(actionWriteSimulationCheckpoint());
    }

    if (interface.isActionRequired(actionReadIterationCheckpoint())){
      PRINT("Loading checkpoint");
      interface.fulfilledAction(actionReadIterationCheckpoint());
    }
    else {
      PRINT("Advancing in time");
    }

    PRINT("Computed time = " << computedTime
          << ", computed timesteps = " << computedTimeSteps);
  }

  double mpi_end_time = MPI_Wtime();
  double mpi_overall_time = mpi_end_time - mpi_start_time;

  PRINT("Time spent computing: " << mpi_compute_time);
  PRINT("Overall time in main computation loop: " << mpi_overall_time);
  PRINT("Ratio Computing/Overall time: " << mpi_compute_time / mpi_overall_time);

  interface.finalize();
  PRINT("Exiting SolverDummy");

  // MPI is finalized in preCICE already

  return 0;
}

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <stdio.h>
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
    conv << "(" << commRank << "/" << commSize << ") "; \
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
  std::cout << "Starting solver dummy..." << std::endl;
  MPI_Init(&argc, &argv);
  int commRank = -1;
  int commSize = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
  MPI_Comm_size(MPI_COMM_WORLD, &commSize);
  PRINT("Rank = " << commRank << ", size = " << commSize);

  using namespace precice;
  using namespace precice::constants;

  if (argc == 1 || (argc != 3 && argc != 8)){
    PRINT("Usage: ./solverproxy configFile proxyName [meshName readDataName writeDataName computationTimeInSeconds N]");
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
    PRINT("  N:                 Number of vertices");
    return 1;
  }
  std::string configFileName(argv[1]);
  std::string proxyName(argv[2]);
  bool readWriteData = false;
  std::string meshName;
  std::string readDataName;
  std::string writeDataName;
  double computationTime = 0.0;
  int N = -1;
  if (argc > 3){
    readWriteData = true;
    PRINT("Reading and writing data");
    meshName = argv[3];
    readDataName = argv[4];
    writeDataName = argv[5];
    computationTime = atof(argv[6]);
    N = atoi(argv[7]);
  }

  SolverInterface interface(proxyName, commRank, commSize);
  interface.configure(configFileName);

  double computedTime = 0.0;
  int computedTimeSteps = 0;

  int meshID = -1;
  int readDataID = -1;
  int writeDataID = -1;
  int dimensions = -1;
  if (readWriteData){
    meshID = interface.getMeshID(meshName);
    readDataID = interface.getDataID(readDataName, meshID);
    writeDataID = interface.getDataID(writeDataName, meshID);
    dimensions = interface.getDimensions();
  }

  if (interface.isActionRequired(actionReadSimulationCheckpoint())){
    interface.fulfilledAction(actionReadSimulationCheckpoint());
  }

  int dataSize = -1;
  double* data = NULL;
  int* dataIndices = NULL;
  if (readWriteData){
    int parallelChunk = N / commSize;
    PRINT("parallelChunk = " << parallelChunk);
    int omittedPart = N - parallelChunk*commSize;
    PRINT("omittedPart = " << omittedPart);

    int addon = 0;
    if (commRank < omittedPart) addon = 1;

    dataSize = parallelChunk + addon;
    data = new double[dataSize*dimensions];
    dataIndices = new int[dataSize];

    int startIndex = 0;
    for (int i=0; i < commRank; i++){
      if (i < omittedPart){
        startIndex += parallelChunk + 1;
      }
      else {
        startIndex += parallelChunk;
      }
    }

    for ( int i=0; i < dataSize; i++){
      double vertex[3];
      vertex[0] = (startIndex + i) * 1.0;
      vertex[1] = 0.0;
      vertex[2] = 0.0;
      dataIndices[i] = interface.setMeshVertex(meshID, vertex);
    }
  }

  double dt = interface.initialize();

  double mpi_start_time = MPI_Wtime();
  double mpi_compute_time = 0.0;

  double mpi_advance_time = 0.0;
  double mpi_read_time = 0.0;
  double mpi_write_time = 0.0;

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

    double mpi_write_time_start = MPI_Wtime();
    if (readWriteData){
      data[dataSize*dimensions-1] = computedTimeSteps;
      interface.writeBlockVectorData(writeDataID, dataSize, dataIndices, data);
    }
    double mpi_write_time_end = MPI_Wtime();
    mpi_write_time += mpi_write_time_end - mpi_write_time_start;

    double mpi_advance_time_start = MPI_Wtime();
    dt = interface.advance(dt);
    double mpi_advance_time_end = MPI_Wtime();
    mpi_advance_time += mpi_advance_time_end - mpi_advance_time_start;

    double mpi_read_time_start = MPI_Wtime();
    if (readWriteData && interface.isReadDataAvailable()){
      interface.readBlockVectorData(readDataID, dataSize, dataIndices, data);
      PRINT("data = " << data[dataSize*dimensions-1]);
    }
    double mpi_read_time_end = MPI_Wtime();
    mpi_read_time += mpi_read_time_end - mpi_read_time_start;

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
  PRINT("Time spent reading: " << mpi_read_time);
  PRINT("Time spent writing: " << mpi_write_time);
  PRINT("Time spent advancing: " << mpi_advance_time);
  PRINT("Sum computing + reading + writing + advancing = "
        << mpi_compute_time + mpi_read_time + mpi_write_time + mpi_advance_time);
  PRINT("MPI Overall time in main computation loop: " << mpi_overall_time);
  PRINT("Ratio Computing/Overall time: " << mpi_compute_time / mpi_overall_time);

  if (commSize > 1){
    double averageComputeTimes = 0.0;
    double maxTime = 0.0;
    MPI_Reduce(&mpi_overall_time, &averageComputeTimes, 1, MPI_DOUBLE,
               MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&mpi_overall_time, &maxTime, 1, MPI_DOUBLE,
               MPI_MAX, 0, MPI_COMM_WORLD);
    averageComputeTimes /= (double)commSize;

    if (commRank == 0){
      PRINT("Average overall time spent = " << averageComputeTimes);
      PRINT("Maximum overall time spent = " << maxTime);
    }
  }

  interface.finalize();
  PRINT("Exiting SolverDummy");

  // MPI is finalized in preCICE already

  return 0;
}

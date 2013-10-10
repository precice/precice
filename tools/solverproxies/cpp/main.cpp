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

  if (argc < 3 || argc > 4){
    PRINT("Usage: ./solverproxy configurationFileName proxyName [computation time in seconds]");
    return 1;
  }
  std::string configFileName(argv[1]);
  std::string dummyName(argv[2]);
  double computationTime = 0.0;
  if (argc == 4){
    computationTime = atof(argv[3]);
  }

  SolverInterface interface(dummyName, comm_rank, comm_size);
  interface.configure(configFileName);

  double computedTime = 0.0;
  int computedTimeSteps = 0;

  if (interface.isActionRequired(actionReadSimulationCheckpoint())){
    interface.fulfilledAction(actionReadSimulationCheckpoint());
  }

  double dt = interface.initialize();

  double mpi_start_time = MPI_Wtime();
  double mpi_compute_time = 0.0;

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

    dt = interface.advance(dt);

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

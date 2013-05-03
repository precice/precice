#include <iostream>
#include <sstream>
#include "../../src/precice/SolverInterface.hpp"
//#include "mpi.h"

/**
 * @brief For printing to the command line.
 *
 * @param message  An input stream such as: "Blabla = " << variable << ", ..."
 */
#define PRINT(message) \
  { \
    std::ostringstream conv; \
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
  PRINT("Starting solver dummy...");
  using namespace precice;
  using namespace precice::constants;

  if (argc != 3){
    PRINT("Usage: ./solverdummy configurationFileName dummyName");
    return 1;
  }
  std::string configFileName(argv[1]);
  std::string dummyName(argv[2]);

  SolverInterface interface(dummyName , 0 ,1);
  interface.configure(configFileName);

  double computedTime = 0.0;
  int computedTimeSteps = 0;

  if (interface.isActionRequired(actionReadSimulationCheckpoint())){
    interface.fulfilledAction(actionReadSimulationCheckpoint());
  }

  double dt = interface.initialize();

  while (interface.isCouplingOngoing()){
    // When an implicit coupling scheme is used, checkpointing is required
    if (interface.isActionRequired(actionWriteIterationCheckpoint())){
      PRINT(actionWriteIterationCheckpoint());
      interface.fulfilledAction(actionWriteIterationCheckpoint());
    }

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

  interface.finalize();
  PRINT("Exiting SolverDummy");

  return 0;
}

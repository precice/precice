#include "utils/Parallel.hpp"
#include "utils/Petsc.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"
#include "precice/config/Configuration.hpp"
#include <iostream>
#include "logging/Logger.hpp"

void printUsage()
{
  std::cout << "Usage:" << std::endl << std::endl;
  std::cout << "Run server (deprecated)  :  binprecice server ParticipantName ConfigurationName [LogConfFile]" << std::endl;
  std::cout << "Print XML reference      :  binprecice xml" << std::endl;
  std::cout << "Print DTD for XML config :  binprecice dtd" << std::endl;
}

int main ( int argc, char** argv )
{
  bool runServer = false;
  bool runHelp = false;
  bool runDtd = false;
  bool hasLogConfFile = false;

  bool wrongParameters = true;

  if (argc >= 2) {
    std::string action(argv[1]);
    if ( action == "dtd" ) {
      wrongParameters = false;
      runDtd = true;
    }
    if ( action == "xml" ) {
      wrongParameters = false;
      runHelp = true;
    }
    if ( action == "server" and argc >= 4 ) {
      wrongParameters = false;
      runServer = true;
      if (argc >= 5){
        hasLogConfFile = true;
      }
    }
  }

  if (wrongParameters) {
    printUsage();
    return 1;
  }

  if (hasLogConfFile){
    precice::logging::setupLogging(argv[3]);
  } else {
    precice::logging::setupLogging();
  }

  precice::utils::Parallel::initializeMPI(&argc, &argv);
  precice::logging::setMPIRank(precice::utils::Parallel::getProcessRank());

  precice::utils::Petsc::initialize(&argc, &argv);

  if ( runServer ){
    assertion(not runHelp);
    std::cout << "PreCICE running server..." << std::endl;
    std::string participantName ( argv[2] );
    std::string configFile ( argv[3] );
    std::cout << "  Participant = " << participantName << std::endl;
    std::cout << "  Configuration = " << configFile << std::endl;
    int size = precice::utils::Parallel::getCommunicatorSize();
    if ( size != 1 ){
      std::cerr << "Server can be run with only one process!" << std::endl;
    }
    precice::impl::SolverInterfaceImpl server ( participantName, 0, 1, true );
    server.configure(configFile);
    server.runServer();
    std::cout << std::endl << std::endl << "...finished running server" << std::endl;
  }
  else if (runHelp){
    assertion(not runServer);
    precice::config::Configuration config;
    std::cout << config.getXMLTag().printDocumentation(0) << std::endl << std::endl;
  }
  else if (runDtd) {
	assertion(not runServer);
    precice::config::Configuration config;
    std::cout << config.getXMLTag().printDTD(true) << std::endl << std::endl;
  }
  else {
    assertion ( false );
  }
  precice::utils::Petsc::finalize();
  //precice::utils::Parallel::synchronizeProcesses();
  //std::cout << "close: " << precice::utils::Parallel::getProcessRank() << std::endl;
  precice::utils::Parallel::finalizeMPI();
  //std::cout << "done" << std::endl;
  return 0;
}

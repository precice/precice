#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
/*
 * ArgumentSetFabric.cpp
 *
 *  Created on: 21.07.2011
 *      Author: michael
 */

#include "tarch/argument/ArgumentSetFabric.h"
#include "tarch/Assertions.h"
#include <sstream>

using namespace tarch::argument;

tarch::logging::Log ArgumentSetFabric::_log("ArgumentSetFabric()");
std::vector<ArgumentSet> ArgumentSetFabric::_argumentSets;
bool ArgumentSetFabric::_isInitialized(false);

ArgumentSetFabric::ArgumentSetFabric() {
}

ArgumentSetFabric::~ArgumentSetFabric() {
  // TODO Auto-generated destructor stub
}

void ArgumentSetFabric::print() {
  std::stringstream ss;
  print(ss);
  _log.info("print",ss.str());
}

void ArgumentSetFabric::print(std::stringstream& ss) {
  for (std::vector<ArgumentSet>::const_iterator it = _argumentSets.begin(); it!=_argumentSets.end(); ++it) {
    (*it).print(ss);
    ss <<std::endl;
  }
}

void ArgumentSetFabric::printDefaultArguments(std::stringstream& ss) {
  for (std::vector<ArgumentSet>::const_iterator it = _argumentSets.begin(); it!=_argumentSets.end(); ++it) {
    (*it).printDefaultArguments(ss);
    ss <<std::endl;
  }
}

void ArgumentSetFabric::printDefaultArguments() {
  std::stringstream ss;
  printDefaultArguments(ss);
  _log.info("print",ss.str());
}

void ArgumentSetFabric::printAllPossibleConfigurationsAndExit() {
  std::stringstream ss;
  ss << "ERROR! ArgumentSet not found!! Available configurations:\n";
  print(ss);
  assertionFail(ss.str());
  _log.error("getArgumentSet",ss.str());
  exit(ASSERTION_EXIT_CODE);
}

ArgumentSet ArgumentSetFabric::getSpecificSet(unsigned int argc, char* argv[]) {
  if(argc<2) {
    printAllPossibleConfigurationsAndExit();
  }

  std::stringstream ss;
  for (std::vector<ArgumentSet>::iterator it = _argumentSets.begin(); it!=_argumentSets.end(); ++it) {
    if((*it).isArgumentSetName(argv[1])) {
#ifdef Debug
      ss << "\nArgumentSet ";
      (*it).printDefaultArguments(ss);
      ss <<"was selected.";
      _log.debug("getArgumentSet", ss.str());
#endif
      (*it).initialize(argc, argv);
      return *it;
    }
  }
  printAllPossibleConfigurationsAndExit();
  return _argumentSets[0];
}

//ArgumentSetFabric::ArgumentSetType ArgumentSetFabric::getCurrentSetting(unsigned int argc, char* argv[]) {
//  if(strcmp(argv[1], "dp")==0&&argc==5){
//      return ArgumentSetFabric::DP;
//  }
//  else if(strcmp(argv[1],"ls-regular")==0&&argc==3) {
//      return ArgumentSetFabric::LS_REGULAR;
//  }
//  else if(strcmp(argv[1],"ls-variant")==0&&argc==4) {
//      return ArgumentSetFabric::LS_VARIANT;
//  }
//  else if(strcmp(argv[1],"run-tests")==0&&argc==2) {
//      return ArgumentSetFabric::RUN_TESTS;
//  }
//  else{
//      return ArgumentSetFabric::UNDEFINED;
//  }
//}

ArgumentSet ArgumentSetFabric::getArgumentSet(unsigned int argc, char* argv[]) {
  if(!_isInitialized) {
    _log.error("getArgumentSet","ArgumentSetFabric is not initialized. You have to add Argumentsets to initialize it before "
        "trying to receiving a configuration.");
  }
#ifdef Debug
  std::stringstream ss;
  ss << "The configured argument sets are: " <<std::endl;
  printDefaultArguments(ss);
  _log.debug("getArgumentSet", ss.str());
#endif

  /*
   * Find user specified argument set
   */
  ArgumentSet argumentSet = getSpecificSet(argc, argv);
  _log.debug("getArgumentSet()","Argument set successfully returned.");
  return argumentSet;
}

void ArgumentSetFabric::addArgumentSet(ArgumentSet argumentSet){
  _isInitialized = true;
  _argumentSets.push_back(argumentSet);
}

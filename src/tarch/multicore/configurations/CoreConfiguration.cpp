#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/multicore/configurations/CoreConfiguration.h"
#include "utils/assertion.hpp"


precice::logging::Logger tarch::multicore::configurations::CoreConfiguration::_log( "tarch::multicore::configurations::CoreConfiguration" );


tarch::multicore::configurations::CoreConfiguration::CoreConfiguration():
  _numberOfThreads(-1),
  _hasParsed(false) {
}


tarch::multicore::configurations::CoreConfiguration::~CoreConfiguration() {
}


std::string tarch::multicore::configurations::CoreConfiguration::getTag() const {
  return "multicore";
}


void tarch::multicore::configurations::CoreConfiguration::parseSubtag(precice::xml::Parser::CTag *pTag) {
  _hasParsed = true;
  
  std::cout << "CoreConfiguration" << std::endl;

  if ( pTag->m_aAttributes.find("threads") == pTag->m_aAttributes.end() ) {
    WARN("attribute \"threads\" missing within tag <" + getTag() + ">");
    _numberOfThreads = -1;
  }
  else if ( std::string("*") == pTag->m_aAttributes["threads"]) ) {
    #ifdef SharedCobra
    WARN("Cobra does not support a default number of threads. Please specify valid number of threads");
    _numberOfThreads = -1;
    #else
    DEBUG("parseSubtag(...)", "use default number of threads");
    _numberOfThreads = 0;
    #endif
  }
  else {
    _numberOfThreads = std::stoi(pTag->m_aAttributes["threads"]);
  }
}

bool tarch::multicore::configurations::CoreConfiguration::isValid() const {
  if (!_hasParsed) {
    assertion( _numberOfThreads<0 );
    WARN("tag <" + getTag() + "> is missing" );
  }
  std::cout << "CoreConfiguration" << std::endl;
  return _numberOfThreads >=0;
}


void tarch::multicore::configurations::CoreConfiguration::toXML(std::ostream& out) const {
  out << "<!--" << std::endl
      << "  Configures the multithreading. This tag has only one attribute and" << std::endl
      << "  the attribute's value is either a natural number of a start to denote" << std::endl
      << "  that the system should use the default settings." << std::endl
      << "  -->" << std::endl;
  if (_numberOfThreads>0) {
    out << "<" << getTag() << " threads=\"" << _numberOfThreads << "\" />";
  }
  else {
    out << "<" << getTag() << " threads=\"*\" />";
  }
}


int tarch::multicore::configurations::CoreConfiguration::getNumberOfThreads() const {
  assertion( isValid() );
  return _numberOfThreads;
}

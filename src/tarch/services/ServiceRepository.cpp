#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/services/ServiceRepository.h"
#include "tarch/Assertions.h"

#include <sstream>


tarch::services::ServiceRepository::ServiceRepository():
  _services(),
  _serviceNames() {
}


tarch::services::ServiceRepository::~ServiceRepository() {
  _services.clear();
}


tarch::services::ServiceRepository& tarch::services::ServiceRepository::getInstance() {
  static tarch::services::ServiceRepository singleton;
  return singleton;
}


void tarch::services::ServiceRepository::addService( Service* const service, const std::string& name ) {
  assertion( service!=nullptr );
  _services.push_back( service );
  _serviceNames.push_back( name );
}


void tarch::services::ServiceRepository::receiveDanglingMessages() {
  for (auto & elem : _services) {
    (elem)->receiveDanglingMessages();
  }
}


std::string tarch::services::ServiceRepository::getListOfRegisteredServices() const {
  std::ostringstream result;
  for (const auto & elem : _serviceNames) {
    result << " " << elem;
  }
  return result.str();
}


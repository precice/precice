// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_SERVICE_SERVICE_FACTORY_H_
#define _TARCH_SERVICE_SERVICE_FACTORY_H_

#ifdef Parallel
#include <mpi.h>
#endif
#include <string>
#include "tarch/services/Service.h"


namespace tarch {
  namespace services {
    template <class ServiceName>
    class ServiceFactory;
  }
}


template <class ServiceName>
class tarch::services::ServiceFactory {
  public:
    ServiceFactory(const std::string& serviceName);
    ~ServiceFactory();
};


#include "tarch/services/ServiceFactory.cpph"

#endif

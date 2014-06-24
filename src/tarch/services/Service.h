// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_SERVICES_SERVICE_H_
#define _TARCH_SERVICES_SERVICE_H_

#define registerService(name) \
  static tarch::services::ServiceFactory<name> thisServiceFactoryInstance( #name );



namespace tarch {
  namespace services {
    class Service;
  }
}


/**
 * Service Interface
 *
 * A service is a singleton running the background which is available always.
 * If it is shut down, all calls to the service becomes nop, i.e. it is a
 * robust thing. On a parallel machine, services are also do exist only once,
 * i.e. they are global singletons. Therefore, you shall never access a service
 * directly, but always via a proxy/accessor. This proxy either invokes the
 * service directly (if the services are on the local node) or it creates an
 * MPI call to the real service. Consequently, services have to poll the MPI
 * queue regulary. To enable them to do so, all services register at the
 * service repository, and this service repository has a poll operation.
 * Services also have to be protected by semaphores in a multicore environment.
 *
 * If you wanna write a service yourself, you have to do two things:
 * - Implement this interface with a class that is a singleton, i.e. it only
 *   has private constructors and a static getInstance() operation.
 * - Add a statement \code
#include "tarch/services/ServiceFactory.h"
registerService(full-qualified-class-name-of-service)
\endcode
 *   to your implementation file.
 *
 * The latter creates the singleton instance and registers the service at the
 * service repository.
 *
 * @author Tobias Weinzierl
 */
class tarch::services::Service {
  public:
    virtual ~Service() {};

    virtual void receiveDanglingMessages() = 0;
};

#endif

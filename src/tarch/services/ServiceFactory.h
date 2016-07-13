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

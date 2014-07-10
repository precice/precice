// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_CONFIGURATION_TOP_LEVEL_CONFIGURATION_FACTORY_H_
#define _TARCH_CONFIGURATION_TOP_LEVEL_CONFIGURATION_FACTORY_H_

#ifdef Parallel
#include <mpi.h>
#endif

#include "tarch/configuration/TopLevelConfiguration.h"
#include "tarch/configuration/BaseTopLevelConfigurationFactory.h"


namespace tarch {
  namespace configuration {
    template <class TopLevelConfigurationType>
    class TopLevelConfigurationFactory;
  }
}

/**
 * Initializes a TopLevelConfiguration class. This is done by
 * the init() method.
 *
 * Note:former versions used to initialize it by the constructor,
 * but this was changed to avoid segmentation faults because of
 * non-initialized static variables.
 */
template <class TopLevelConfigurationType>
class tarch::configuration::TopLevelConfigurationFactory
: public tarch::configuration::BaseTopLevelConfigurationFactory
{
  private:
    TopLevelConfiguration* _configuration;
  public:
    /**
     * Constructor for TopLevelConfigurationFactory.
     *
     * When the constructor is called, it will register this class
     * in the configuration registry
     */
    TopLevelConfigurationFactory();

    virtual ~TopLevelConfigurationFactory();

    /**
     * Creates an instance of the template argument, and adds it
     * to the ConfigurationRegistry.
     */
    virtual void init();
};


#include "tarch/configuration/TopLevelConfigurationFactory.cpph"

#endif


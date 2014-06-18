// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_MULTICORE_TBB_CONFIGURATION_CORE_CONFIGURATION_H_
#define _TARCH_MULTICORE_TBB_CONFIGURATION_CORE_CONFIGURATION_H_

#include "tarch/configuration/Configuration.h"
#include "tarch/logging/Log.h"
#include <string>


namespace tarch {
  namespace multicore {
    namespace configurations {
      class CoreConfiguration;
    }
  }
}


/**
 * Core Configuration
 *
 * @author Tobias Weinzierl
 */
class tarch::multicore::configurations::CoreConfiguration: public tarch::configuration::Configuration {
  private:
    /**
     * Log device.
     */
    static tarch::logging::Log _log;

    int _numberOfThreads;

    bool _hasParsed;

  public:
    CoreConfiguration();
    virtual ~CoreConfiguration();

    virtual std::string getTag() const;

    virtual void parseSubtag( tarch::irr::io::IrrXMLReader* xmlReader );

    virtual bool isValid() const;

    virtual void toXML(std::ostream& out) const;

    /**
     * @return Number of threads or 0 if you should use the standard number.
     */
    int getNumberOfThreads() const;
};


#endif

// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _PLOTTER_CONFIGURATION_PLOTTERCONFIGURATION_H_
#define _PLOTTER_CONFIGURATION_PLOTTERCONFIGURATION_H_


#include "tarch/configuration/Configuration.h"
#include "tarch/logging/Log.h"



#include <string>

namespace tarch {
  namespace plotter {
    namespace griddata {
      namespace unstructured {
        namespace configuration {
          class UnstructuredGridPlotterConfiguration;
        }
      }
    }
  }
}



/**
 * Plotter Configuration
 *
 * Represents the configuration for one plotter for unstructered grids.
 *
 * @author Tobias Weinzierl
 */
class tarch::plotter::griddata::unstructured::configuration::UnstructuredGridPlotterConfiguration:
  public tarch::configuration::Configuration {
  private:
    /**
     * Logging device.
     */
    static tarch::logging::Log _log;

    const static std::string PLOTTER_IDENTIFIER_VTK_TEXTFILE;
    const static std::string ATTRIBUTE_PATH;
    const static std::string ATTRIBUTE_FILENAME;

  public:
    UnstructuredGridPlotterConfiguration();
    virtual ~UnstructuredGridPlotterConfiguration();

    virtual void parseSubtag( tarch::irr::io::IrrXMLReader* xmlReader );
    virtual std::string getTag() const;
    virtual bool isValid() const;
    virtual std::string toXML() const;

    /**
     * Factory mechanism. You are responsible to delete the instance afterwards.
     */
    tarch::plotter::griddata::unstructured::UnstructuredGridWriter* createPlotter() const;
};

#endif

// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_CONFIGURATION_CONFIGURATION_H_
#define _TARCH_CONFIGURATION_CONFIGURATION_H_

#ifdef Parallel
#include <mpi.h>
#endif

#include "tarch/irr/XML.h"
#include <sstream>


namespace tarch {
  namespace configuration {
    class Configuration;
  }
}

/**
 * Abstract supertype of all interfaces
 *
 * @author Tobias Weinzierl
 * @version $Revision: 1.21 $
 */
class tarch::configuration::Configuration {
  public:
    /**
     * Destructor
     */
    virtual ~Configuration() {}

    /**
     * Parse a (sub)tag.
     *
     * Parse a tag corresponding to this configuration within the xml file.
     * The argument passed (the reader) may not be 0 and the reader's state has
     * to be 'I've just read a tag whose name equals the TAG of this class'. The
     * operation terminates as soon as the reader passed has read the closing
     * TAG-Tag. On any syntactic error the operation terminates immediately, and
     * any successing isValid() call fails. If something is wrong within the
     * configuration, parseSubtag() should write a (detailed) error description.
     *
     * @param xmlReader Reader to be used.
     */
    virtual void parseSubtag( tarch::irr::io::IrrXMLReader* xmlReader ) = 0;

    /**
     * Return name of xml tag that is associated to the configuration.
     */
    virtual std::string getTag() const = 0;

    /**
     * Is config valid?
     *
     * This operation usually fails, if
     *
     * - parseSubtag() hasn't been called, i.e. configuration has not been
     *   used, or
     * - parseSubtag() failed due to a wrong file.
     *
     * If a tag ain't optional and parseSubtag() was not called (first case)
     */
    virtual bool isValid() const = 0;

    /**
     * Create xml representation of configuration. This description also should
     * comprise some documentation on the required fields, possible attribute
     * values, and subtags. So, if you invoke toXML() on an invalid
     * (empty) configuration, it gives you a description of a dummy xml file.
     *
     * @param out Output stream to write the xml description to.
     */
    virtual void toXML(std::ostream& out) const = 0;
};


#endif


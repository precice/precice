#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/logging/configurations/LogFilterConfiguration.h"
#include "utils/assertion.hpp"

#include <cstring>
#include <sstream>


tarch::logging::Log tarch::logging::configurations::LogFilterConfiguration::_log("tarch::logging::configurations::LogFilterConfiguration");


tarch::logging::configurations::LogFilterConfiguration::LogFilterConfiguration():
  _isValid(true) {
}


tarch::logging::configurations::LogFilterConfiguration::~LogFilterConfiguration() {
}


std::string tarch::logging::configurations::LogFilterConfiguration::getTag() const {
  return "log-filter";
}


void tarch::logging::configurations::LogFilterConfiguration::parseSubtag( tarch::irr::io::IrrXMLReader* xmlReader ) {
  assertion( xmlReader != nullptr );

  if (xmlReader->getAttributeValue("target")==nullptr) {
    _log.error("parseSubtag(...)", "attribute \"target\" missing within tag <" + getTag() + ">");
    _isValid = false;
    return;
  }

  if (xmlReader->getAttributeValue("switch")==nullptr) {
    _log.error("parseSubtag(...)", "attribute \"switch\" missing within tag <" + getTag() + ">");
    _isValid = false;
    return;
  }

  if (xmlReader->getAttributeValue("component")==nullptr) {
    _log.error("parseSubtag(...)", "attribute \"component\" missing within tag <" + getTag() + ">");
    _isValid = false;
    return;
  }

  if (strcmp("debug", xmlReader->getAttributeValue("target")) &&
      strcmp("info", xmlReader->getAttributeValue("target")) &&
      strcmp("", xmlReader->getAttributeValue("target")))
  {
  	_log.error("parseSubtag(...)", "only value \"debug\", \"info\", or \"\" allowed for \"target\"-attribute within tag <" + getTag() + ">");
  	_isValid = false;
  	return;
  }

  tarch::logging::CommandLineLogger::FilterListEntry newEntry;
  newEntry._targetName = xmlReader->getAttributeValue("target");

  if (!strcmp("on", xmlReader->getAttributeValue("switch"))) {
    newEntry._isBlackEntry = false;
  }
  else if (!strcmp("off", xmlReader->getAttributeValue("switch"))) {
    newEntry._isBlackEntry = true;
  }
  else {
    _log.error("parseSubtag(...)", "only value \"on\" or \"off\" allowed for \"switch\"-attribute within tag <" + getTag() + ">");
    _isValid = false;
    return;
  }

  if ( (xmlReader->getAttributeValue("rank")!=nullptr) && strcmp("*",xmlReader->getAttributeValue("rank")) ) {
  	newEntry._rank = xmlReader->getAttributeValueAsInt("rank");
  }
  else {
    newEntry._rank=-1;
  }

  newEntry._namespaceName = xmlReader->getAttributeValue("component");

  if ( _filterList.count(newEntry)!=0 ) {
    logError( "parseSubtag(...)", "tried to insert " << newEntry.toString() << " multiple times");
    _isValid = false;
  }
  else {
    _filterList.insert( newEntry );
  }
}


bool tarch::logging::configurations::LogFilterConfiguration::isValid() const {
  return _isValid;
}


tarch::logging::CommandLineLogger::FilterList tarch::logging::configurations::LogFilterConfiguration::getFilterList() const {
  return _filterList;
}


void tarch::logging::configurations::LogFilterConfiguration::toXML(std::ostream& out) const {
  out << "<!--" << std::endl
      << "  This is the configuration tag corresponding to tarch::logging::configurations::LogFilterConfiguration. " << std::endl
      << "  A log filter tag has to contain the following attributes " << std::endl
      << std::endl
      << "    | attribute name | semantics | allowed values " << std::endl
      << "    | target         | Log level to write. | At the time, only debug is supported." << std::endl
      << "    | component      | namespace/class whose outputs are to be written. | Any namespace/class. " << std::endl
      << "    | rank           | Rank of the node from where the message has to be written. | Any positive number of * for all nodes. " << std::endl
      << "    | switch         | Block or let pass. | on or off " << std::endl
      << std::endl
      << "  and there always might be several of them. They define the filter list " << std::endl
      << "  entries of the logger. The tag may not contain any subtags. " << std::endl
      << "  -->" << std::endl;
  if ( isValid() && !_filterList.empty() ) {
    for (const auto & elem : _filterList) {
      out << "<" + getTag();
      out << " target=\"debug\"";
      out << " component=\"" << elem._namespaceName << "\"";
      out << " rank=\""      << elem._rank          << "\"";
      out << " switch=\"off\"";
      out << "/>\n";
    }
  }
  else {
    out << "<" + getTag();
    out << " target=\"debug\"";
    out << " component=\"a class/namespace name\"";
    out << " rank=\"14\"";
    out << " switch=\"off\"";
    out << "/>\n";
  }
}

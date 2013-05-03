// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "LogFilterConfiguration.hpp"
#include "utils/xml/ValidatorEquals.hpp"
#include "utils/xml/ValidatorOr.hpp"

namespace precice {
namespace config {

tarch::logging::Log LogFilterConfiguration::
  _log("precice::config::LogFilterConfiguration");

LogFilterConfiguration:: LogFilterConfiguration
(
  utils::XMLTag& parent )
:
  TAG("log-filter"),
  ATTR_TARGET("target"),
  ATTR_COMPONENT("component"),
  ATTR_SWITCH("switch"),
  //_isValid(false),
  _filterList()
{
  using namespace utils;
  std::string doc;
  XMLTag tag(*this, TAG, XMLTag::OCCUR_ARBITRARY);
  doc = "Filters debug and/or info messages written by preCICE. ";
  doc += "Debug messages are only written when preCICE is compiled in debug mode.";
  tag.setDocumentation(doc);

  XMLAttribute<std::string> attrTarget(ATTR_TARGET);
  doc = "Specifies the message type/s to be filtered. ";
  doc += "An empty string means debug and info.";
  attrTarget.setDocumentation(doc);
  ValidatorEquals<std::string> validDebug("debug");
  ValidatorEquals<std::string> validInfo("info");
  ValidatorEquals<std::string> validAll("");
  attrTarget.setValidator(validDebug || validInfo || validAll);
  tag.addAttribute(attrTarget);

  XMLAttribute<std::string> attrComponent(ATTR_COMPONENT);
  doc = "Specifies the component to be filtered, which is visible when log-trace";
  doc += " in the <log-output> tag is set to true. It is possible to specify ";
  doc += "only a part of the full component name (e.g. just the namespace) to ";
  doc += "address all components starting with that part. The empty string means";
  doc += "all components.";
  attrComponent.setDocumentation(doc);
  tag.addAttribute(attrComponent);

  XMLAttribute<bool> attrSwitch(ATTR_SWITCH);
  doc = "If switch is set to \"on\", the messages are enabled, if set to \"off\"";
  doc += " messages are disabled. By default, all messages are disabled.";
  attrSwitch.setDocumentation(doc);
  tag.addAttribute(attrSwitch);

  parent.addSubtag(tag);
}

//void tarch::logging::configurations::LogFilterConfiguration::parseSubtag( tarch::irr::io::IrrXMLReader* xmlReader ) {
//  assertion( xmlReader != 0 );
//
//  if (xmlReader->getAttributeValue("target")==0) {
//    _log.error("parseSubtag(...)", "attribute \"target\" missing within tag <" + getTag() + ">");
//    _isValid = false;
//    return;
//  }
//
//  if (xmlReader->getAttributeValue("switch")==0) {
//    _log.error("parseSubtag(...)", "attribute \"switch\" missing within tag <" + getTag() + ">");
//    _isValid = false;
//    return;
//  }
//
//  if (xmlReader->getAttributeValue("component")==0) {
//    _log.error("parseSubtag(...)", "attribute \"component\" missing within tag <" + getTag() + ">");
//    _isValid = false;
//    return;
//  }
//
//  if (strcmp("debug", xmlReader->getAttributeValue("target")) &&
//      strcmp("info", xmlReader->getAttributeValue("target")) &&
//      strcmp("", xmlReader->getAttributeValue("target")))
//  {
//  	_log.error("parseSubtag(...)", "only value \"debug\", \"info\", or \"\" allowed for \"target\"-attribute within tag <" + getTag() + ">");
//  	_isValid = false;
//  	return;
//  }
//
//  tarch::logging::CommandLineLogger::FilterListEntry newEntry;
//  newEntry._targetName = xmlReader->getAttributeValue("target");
//
//  if (!strcmp("on", xmlReader->getAttributeValue("switch"))) {
//    newEntry._isBlackEntry = false;
//  }
//  else if (!strcmp("off", xmlReader->getAttributeValue("switch"))) {
//    newEntry._isBlackEntry = true;
//  }
//  else {
//    _log.error("parseSubtag(...)", "only value \"on\" or \"off\" allowed for \"switch\"-attribute within tag <" + getTag() + ">");
//    _isValid = false;
//    return;
//  }
//
//  if ( (xmlReader->getAttributeValue("rank")!=0) && strcmp("*",xmlReader->getAttributeValue("rank")) ) {
//  	newEntry._rank = xmlReader->getAttributeValueAsInt("rank");
//  }
//  else {
//    newEntry._rank=-1;
//  }
//
//  newEntry._namespaceName = xmlReader->getAttributeValue("component");
//
//  if ( _filterList.count(newEntry)!=0 ) {
//    logError( "parseSubtag(...)", "tried to insert " << newEntry.toString() << " multiple times");
//    _isValid = false;
//  }
//  else {
//    _filterList.insert( newEntry );
//  }
//}


//bool LogFilterConfiguration:: isValid() const
//{
//  return _isValid;
//}

tarch::logging::CommandLineLogger::FilterList LogFilterConfiguration:: getFilterList() const
{
  return _filterList;
}

void LogFilterConfiguration:: xmlTagCallback
(
  utils::XMLTag& tag )
{
  preciceTrace1("xmlTagCallback()", tag.getFullName());
  if (tag.getName() == TAG){
    tarch::logging::CommandLineLogger::FilterListEntry newEntry;
    newEntry._targetName = tag.getStringAttributeValue(ATTR_TARGET);
    newEntry._isBlackEntry = not tag.getBooleanAttributeValue(ATTR_SWITCH);
    newEntry._rank = -1; // no ranks configured
    newEntry._namespaceName = tag.getStringAttributeValue(ATTR_COMPONENT);
    _filterList.insert(newEntry);
  }
}

void LogFilterConfiguration:: xmlEndTagCallback
(
  utils::XMLTag& tag )
{
  preciceTrace1("xmlEndTagCallback()", tag.getFullName());
}

//void tarch::logging::configurations::LogFilterConfiguration::toXML(std::ostream& out) const {
//  out << "<!--" << std::endl
//      << "  This is the configuration tag corresponding to tarch::logging::configurations::LogFilterConfiguration. " << std::endl
//      << "  A log filter tag has to contain the following attributes " << std::endl
//      << std::endl
//      << "    | attribute name | semantics | allowed values " << std::endl
//      << "    | target         | Log level to write. | At the time, only debug is supported." << std::endl
//      << "    | component      | namespace/class whose outputs are to be written. | Any namespace/class. " << std::endl
//      << "    | rank           | Rank of the node from where the message has to be written. | Any positive number of * for all nodes. " << std::endl
//      << "    | switch         | Block or let pass. | on or off " << std::endl
//      << std::endl
//      << "  and there always might be several of them. They define the filter list " << std::endl
//      << "  entries of the logger. The tag may not contain any subtags. " << std::endl
//      << "  -->" << std::endl;
//  if ( isValid() && !_filterList.empty() ) {
//    for (tarch::logging::CommandLineLogger::FilterList::const_iterator p = _filterList.begin(); p!=_filterList.end(); p++ ) {
//      out << "<" + getTag();
//      out << " target=\"debug\"";
//      out << " component=\"" << p->_namespaceName << "\"";
//      out << " rank=\""      << p->_rank          << "\"";
//      out << " switch=\"off\"";
//      out << "/>\n";
//    }
//  }
//  else {
//    out << "<" + getTag();
//    out << " target=\"debug\"";
//    out << " component=\"a class/namespace name\"";
//    out << " rank=\"14\"";
//    out << " switch=\"off\"";
//    out << "/>\n";
//  }
//}

}} // namespace precice, config

// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "LogOutputFormatConfiguration.hpp"

namespace precice {
namespace config {

tarch::logging::Log LogOutputFormatConfiguration:: _log("precice::config::LogOutputFormatConfiguration");

LogOutputFormatConfiguration:: LogOutputFormatConfiguration
(
  utils::XMLTag& parent )
:
  TAG("log-output"),
  ATTR_COLUMN_SEPARATOR("column-separator"),
  ATTR_LOG_TIMESTAMP("log-time-stamp"),
  ATTR_LOG_TIMESTAMP_HUMAN_READABLE("log-time-stamp-human-readable"),
  ATTR_LOG_MACHINE_NAME("log-machine-name"),
  ATTR_LOG_MESSAGE_TYPE("log-message-type"),
  ATTR_LOG_TRACE("log-trace"),
  _logColumnSeparator(""),
  _logTimeStamp(false),
  _logTimeStampHumanReadable(false),
  _logMachineName(false),
  _logMessageType(false),
  _logTrace(false)
  //_isValid(false)
{
  using namespace utils;
  std::string doc;
  XMLTag tag(*this, TAG, XMLTag::OCCUR_NOT_OR_ONCE);
  doc = "Specifies the format of messages output by preCICE.";
  tag.setDocumentation(doc);

  XMLAttribute<std::string> attrColumnSeparator(ATTR_COLUMN_SEPARATOR);
  doc = "Defines the speparating characters for the content types of a message.";
  attrColumnSeparator.setDocumentation(doc);
  tag.addAttribute(attrColumnSeparator);

  XMLAttribute<bool> attrLogTimestamp(ATTR_LOG_TIMESTAMP);
  doc = "If set to \"on\", the current time is added to a message.";
  attrLogTimestamp.setDocumentation(doc);
  tag.addAttribute(attrLogTimestamp);

  XMLAttribute<bool> attrLogTimestampHR(ATTR_LOG_TIMESTAMP_HUMAN_READABLE);
  doc = "If set to \"on\", the time is printed in a human-readable format";
  attrLogTimestampHR.setDocumentation(doc);
  tag.addAttribute(attrLogTimestampHR);

  XMLAttribute<bool> attrLogMachineName(ATTR_LOG_MACHINE_NAME);
  doc = "If set to on, the machine name is added to a message.";
  attrLogMachineName.setDocumentation(doc);
  tag.addAttribute(attrLogMachineName);

  XMLAttribute<bool> attrLogMessageType(ATTR_LOG_MESSAGE_TYPE);
  doc = "If set to on, the type of a message (debug, info, ...) is added to it.";
  attrLogMessageType.setDocumentation(doc);
  tag.addAttribute(attrLogMessageType);

  XMLAttribute<bool> attrLogTrace(ATTR_LOG_TRACE);
  doc = "If set to on, the component of a message, i.e. the full name of the place";
  doc += " the message was written from is added to it.";
  attrLogTrace.setDocumentation(doc);
  tag.addAttribute(attrLogTrace);

  parent.addSubtag(tag);
}

//std::string tarch::logging::configurations::LogOutputFormatConfiguration::getTag() const
//{
//  return "log-output";
//}


//void tarch::logging::configurations::LogOutputFormatConfiguration::parseSubtag( tarch::irr::io::IrrXMLReader* _xmlReader )
//{
//  assertion( _xmlReader != 0 );
//
//  _isValid   = true;
//  _hasParsed = true;
//
//  if (_xmlReader->getAttributeValue("log-time-stamp") != 0) {
//    _logTimeStamp = _xmlReader->getAttributeValueAsBool("log-time-stamp");
//  }
//  else {
//    _log.error( "parseSubtag(...)", "attribute \"log-time-stamp\" missing within tag " + getTag() );
//    _isValid = false;
//  }
//
//  if (_xmlReader->getAttributeValue("log-time-stamp-human-readable") != 0) {
//    _logTimeStampHumanReadable = _xmlReader->getAttributeValueAsBool("log-time-stamp-human-readable");
//  }
//  else {
//    _log.error( "parseSubtag(...)", "attribute \"log-time-stamp-human-readable\" missing within tag " + getTag() );
//    _isValid = false;
//  }
//
//  if (_xmlReader->getAttributeValue("log-machine-name") != 0) {
//    _logMachineName = _xmlReader->getAttributeValueAsBool("log-machine-name");
//  }
//  else {
//    _log.error( "parseSubtag(...)", "attribute \"log-machine-name\" missing within tag " + getTag() );
//    _isValid = false;
//  }
//
//  if (_xmlReader->getAttributeValue("log-message-type") != 0) {
//    _logMessageType = _xmlReader->getAttributeValueAsBool("log-message-type");
//  }
//  else {
//    _log.error( "parseSubtag(...)", "attribute \"log-message-type\" missing within tag " + getTag() );
//    _isValid = false;
//  }
//
//  if (_xmlReader->getAttributeValue("log-trace") != 0) {
//    _logTrace = _xmlReader->getAttributeValueAsBool("log-trace");
//  }
//  else {
//    _log.error( "parseSubtag(...)", "attribute \"log-trace\" missing within tag " + getTag() );
//    _isValid = false;
//  }
//}


//void tarch::logging::configurations::LogOutputFormatConfiguration::toXML(std::ostream& out) const {
//  out << "<!--" << std::endl
//      << "  This is the configuration tag corresponding to tarch::logging::configurations::LogOutputFormatConfiguration. " << std::endl
//      << "  The xml tag has to contain the following attributes: " << std::endl
//      << std::endl
//      << "    | attribute name | semantics | allowed values " << std::endl
//      << "    | column-separator | What string to use to separate different columns. | I typically use a blank. " << std::endl
//      << "    | log-time-stamp   | Shall log add a time stamp as machine counter. | yes or no" << std::endl
//      << "    | log-time-stamp-human-readable | Shall log add a time stamp that is human readable (hh::mm::ss). | yes or no" << std::endl
//      << "    | log-machine-name | Shall log add the node's name. | yes or no " << std::endl
//      << "    | log-message-type | Shall log add the message type (debug, info, warning, and so forth). | yes or no " << std::endl
//      << "    | log-trace        | Shall log add a trace, i.e. who printed the message. | yes or no " << std::endl
//      << std::endl
//      << "  The xml tag my not contain any subtags. " << std::endl
//      << "  -->" << std::endl;
//  out << "<" + getTag() << " ";
//  out << "column-separator=\"" << _logColumnSeparator << "\" ";
//  out << "log-time-stamp=\"" << _logTimeStamp << "\" ";
//  out << "log-time-stamp-human-readable=\"" << _logTimeStampHumanReadable << "\" ";
//  out << "log-machine-name=\"" << _logMachineName << "\" ";
//  out << "log-message-type=\"" << _logMessageType << "\" ";
//  out << "log-trace=\"" << _logTrace << "\" ";
//  out << "/>" << std::endl;
//}

//bool tarch::logging::configurations::LogOutputFormatConfiguration::hasParsed() const {
//  return _hasParsed;
//}


//bool LogOutputFormatConfiguration:: isValid() const
//{
//  return _isValid;
//}


const std::string& LogOutputFormatConfiguration:: getLogColumnSeparator() const
{
  //assertion( isValid() );
  return _logColumnSeparator;
}


bool LogOutputFormatConfiguration:: getLogTimeStamp() const
{
  //assertion( isValid() );
  return _logTimeStamp;
}


bool LogOutputFormatConfiguration:: getLogTimeStampHumanReadable() const
{
  //assertion( isValid() );
  return _logTimeStampHumanReadable;
}


bool LogOutputFormatConfiguration:: getLogMachineName() const
{
  //assertion( isValid() );
  return _logMachineName;
}


bool LogOutputFormatConfiguration:: getLogMessageType() const
{
  //assertion( isValid() );
  return _logMessageType;
}


bool LogOutputFormatConfiguration:: getLogTrace() const
{
  //assertion( isValid() );
  return _logTrace;
}

void LogOutputFormatConfiguration:: xmlTagCallback
(
  utils::XMLTag& tag )
{
  if (tag.getName() == TAG){
    _logColumnSeparator = tag.getStringAttributeValue(ATTR_COLUMN_SEPARATOR);
    _logTimeStamp = tag.getBooleanAttributeValue(ATTR_LOG_TIMESTAMP);
    _logTimeStampHumanReadable = tag.getBooleanAttributeValue(ATTR_LOG_TIMESTAMP_HUMAN_READABLE);
    _logMachineName = tag.getBooleanAttributeValue(ATTR_LOG_MACHINE_NAME);
    _logMessageType = tag.getBooleanAttributeValue(ATTR_LOG_MESSAGE_TYPE);
    _logTrace = tag.getBooleanAttributeValue(ATTR_LOG_TRACE);
    //_isValid = true;
  }
}

void LogOutputFormatConfiguration:: xmlEndTagCallback
(
  utils::XMLTag& callingTag )
{
}

}} // precice, config

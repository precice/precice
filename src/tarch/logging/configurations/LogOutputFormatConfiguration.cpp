#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/logging/configurations/LogOutputFormatConfiguration.h"
#include "tarch/Assertions.h"


tarch::logging::Log tarch::logging::configurations::LogOutputFormatConfiguration::_log("tarch::logging::configurations::LogOutputFormatConfiguration");


tarch::logging::configurations::LogOutputFormatConfiguration::LogOutputFormatConfiguration():
  _isValid(false),
  _hasParsed(false) {
}


tarch::logging::configurations::LogOutputFormatConfiguration::~LogOutputFormatConfiguration() {
}


std::string tarch::logging::configurations::LogOutputFormatConfiguration::getTag() const {
  return "log-output";
}


void tarch::logging::configurations::LogOutputFormatConfiguration::parseSubtag( tarch::irr::io::IrrXMLReader* _xmlReader ) {
  assertion( _xmlReader != 0 );

  _isValid   = true;
  _hasParsed = true;

  if (_xmlReader->getAttributeValue("column-separator") != 0) {
    _logColumnSeparator = _xmlReader->getAttributeValue("column-separator");
  }
  else {
    _log.error( "parseSubtag(...)", "attribute \"column-separator\" missing within tag " + getTag() );
    _isValid = false;
  }

  if (_xmlReader->getAttributeValue("log-time-stamp") != 0) {
    _logTimeStamp = _xmlReader->getAttributeValueAsBool("log-time-stamp");
  }
  else {
    _log.error( "parseSubtag(...)", "attribute \"log-time-stamp\" missing within tag " + getTag() );
    _isValid = false;
  }

  if (_xmlReader->getAttributeValue("log-output-file") != 0) {
    _logOutputFile = _xmlReader->getAttributeValue("log-output-file");
  }
  else {
    _log.error( "parseSubtag(...)", "attribute \"log-output-file\" missing within tag " + getTag() + ". Can be empty string, but has to exist" );
    _isValid = false;
  }

  if (_xmlReader->getAttributeValue("log-time-stamp-human-readable") != 0) {
    _logTimeStampHumanReadable = _xmlReader->getAttributeValueAsBool("log-time-stamp-human-readable");
  }
  else {
    _log.error( "parseSubtag(...)", "attribute \"log-time-stamp-human-readable\" missing within tag " + getTag() );
    _isValid = false;
  }

  if (_xmlReader->getAttributeValue("log-machine-name") != 0) {
    _logMachineName = _xmlReader->getAttributeValueAsBool("log-machine-name");
  }
  else {
    _log.error( "parseSubtag(...)", "attribute \"log-machine-name\" missing within tag " + getTag() );
    _isValid = false;
  }

  if (_xmlReader->getAttributeValue("log-message-type") != 0) {
    _logMessageType = _xmlReader->getAttributeValueAsBool("log-message-type");
  }
  else {
    _log.error( "parseSubtag(...)", "attribute \"log-message-type\" missing within tag " + getTag() );
    _isValid = false;
  }

  if (_xmlReader->getAttributeValue("log-trace") != 0) {
    _logTrace = _xmlReader->getAttributeValueAsBool("log-trace");
  }
  else {
    _log.error( "parseSubtag(...)", "attribute \"log-trace\" missing within tag " + getTag() );
    _isValid = false;
  }
}


void tarch::logging::configurations::LogOutputFormatConfiguration::toXML(std::ostream& out) const {
  out << "<!--" << std::endl
      << "  This is the configuration tag corresponding to tarch::logging::configurations::LogOutputFormatConfiguration. " << std::endl
      << "  The xml tag has to contain the following attributes: " << std::endl
      << std::endl
      << "    | attribute name | semantics | allowed values " << std::endl
      << "    | column-separator | What string to use to separate different columns. | I typically use a blank. " << std::endl
      << "    | log-time-stamp   | Shall log add a time stamp as machine counter. | yes or no" << std::endl
      << "    | log-time-stamp-human-readable | Shall log add a time stamp that is human readable (hh::mm::ss). | yes or no" << std::endl
      << "    | log-machine-name | Shall log add the node's name. | yes or no " << std::endl
      << "    | log-message-type | Shall log add the message type (debug, info, warning, and so forth). | yes or no " << std::endl
      << "    | log-trace        | Shall log add a trace, i.e. who printed the message. | yes or no " << std::endl
      << std::endl
      << "  The xml tag my not contain any subtags. " << std::endl
      << "  -->" << std::endl;
  out << "<" + getTag() << " ";
  out << "column-separator=\"" << _logColumnSeparator << "\" ";
  out << "log-time-stamp=\"" << _logTimeStamp << "\" ";
  out << "log-time-stamp-human-readable=\"" << _logTimeStampHumanReadable << "\" ";
  out << "log-machine-name=\"" << _logMachineName << "\" ";
  out << "log-message-type=\"" << _logMessageType << "\" ";
  out << "log-trace=\"" << _logTrace << "\" ";
  out << "/>" << std::endl;
}

bool tarch::logging::configurations::LogOutputFormatConfiguration::hasParsed() const {
  return _hasParsed;
}


bool tarch::logging::configurations::LogOutputFormatConfiguration::isValid() const {
  if (!_hasParsed) {
    assertion( !_isValid );
    logError( "isValid()", "tag <" + getTag() + "> is missing" );
  }
  return _isValid;
}


std::string tarch::logging::configurations::LogOutputFormatConfiguration::getLogColumnSeparator() const {
  assertion( isValid() );
  return _logColumnSeparator;
}


bool tarch::logging::configurations::LogOutputFormatConfiguration::getLogTimeStamp() const {
  assertion( isValid() );
  return _logTimeStamp;
}


bool tarch::logging::configurations::LogOutputFormatConfiguration::getLogTimeStampHumanReadable() const {
  assertion( isValid() );
  return _logTimeStampHumanReadable;
}


bool tarch::logging::configurations::LogOutputFormatConfiguration::getLogMachineName() const {
  assertion( isValid() );
  return _logMachineName;
}


bool tarch::logging::configurations::LogOutputFormatConfiguration::getLogMessageType() const  {
  assertion( isValid() );
  return _logMessageType;
}


std::string tarch::logging::configurations::LogOutputFormatConfiguration::getLogFileName() const {
  assertion( isValid() );

  if (_logOutputFile.empty()) {
    return "";
  }
  else {
    std::ostringstream outputFileName;
    #ifdef Parallel
    outputFileName << "rank-"
                   << tarch::parallel::Node::getInstance().getRank()
                   << "-";
    #endif
    outputFileName << _logOutputFile;

    return outputFileName.str();
  }
}


bool tarch::logging::configurations::LogOutputFormatConfiguration::getLogTrace() const {
  assertion( isValid() );
  return _logTrace;
}

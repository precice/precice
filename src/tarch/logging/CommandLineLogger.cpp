#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/logging/CommandLineLogger.h"

#ifdef Parallel
#include "tarch/parallel/Node.h"
#endif

#include "tarch/Assertions.h"

#include "tarch/multicore/Lock.h"

#include <sstream>
#include <fstream>


#include "tarch/services/ServiceFactory.h"
registerService(tarch::logging::CommandLineLogger)



tarch::logging::Log tarch::logging::CommandLineLogger::_log( "tarch::logging::CommandLineLogger" );


std::string::size_type tarch::logging::CommandLineLogger::_indent = 0;


const std::string::size_type tarch::logging::CommandLineLogger::NumberOfIndentSpaces         = 2;
const std::string::size_type tarch::logging::CommandLineLogger::NumberOfTraceColumnSpaces    = 55;
const std::string::size_type tarch::logging::CommandLineLogger::NumberOfStandardColumnSpaces = 12;


tarch::logging::CommandLineLogger::FilterListEntry::FilterListEntry( const std::string& targetName, bool isBlackListEntry ):
  _targetName(targetName),
  _rank(-1),
  _namespaceName(""),
  _isBlackEntry(isBlackListEntry) {
  assertion1( targetName==std::string("info") || targetName==std::string("debug") || targetName==std::string(""), targetName );
}


tarch::logging::CommandLineLogger::FilterListEntry::FilterListEntry( const std::string& targetName, int rank, const std::string& className, bool isBlackListEntry ):
  _targetName(targetName),
  _rank(rank),
  _namespaceName(className),
  _isBlackEntry(isBlackListEntry) {
  assertion1( targetName==std::string("info") || targetName==std::string("debug") || targetName==std::string(""), targetName );
}


std::string tarch::logging::CommandLineLogger::FilterListEntry::toString() const {
  std::ostringstream msg;
  msg << "(";
  if (_targetName=="") {
    msg << "target:*,";
  }
  else {
    msg << _targetName << ",";
  }
  if (_namespaceName=="") {
    msg << "namespace:*,";
  }
  else {
    msg << _namespaceName << ",";
  }
  #ifdef Parallel
  msg << "rank:";
  if (_rank!=-1) {
    msg << _rank;
  }
  else {
    msg << "*";
  }
  msg << ",";
  #endif
  if (_isBlackEntry) {
    msg << "blacklist-entry";
  }
  else {
    msg << "whitelist-entry";
  }
  msg << ")";
  return msg.str();
}


bool tarch::logging::CommandLineLogger::FilterListEntry::operator<(const FilterListEntry& b) const {
  return _rank < b._rank || _namespaceName < b._namespaceName || _targetName < b._targetName;
}


bool tarch::logging::CommandLineLogger::FilterListEntry::operator==(const FilterListEntry& b) const {
  return _rank == b._rank || _namespaceName == b._namespaceName || _targetName == b._targetName;
}


tarch::logging::CommandLineLogger::CommandLineLogger():
  _outputStream(nullptr) {
  configureOutputStreams();
  setLogColumnSeparator();
  setLogTimeStamp();
  setLogTimeStampHumanReadable();
  setLogMachineName();
  setLogMessageType();
  setLogTrace();

  #ifdef Debug
  addFilterListEntry( FilterListEntry( "debug", false )  );
  #endif
  addFilterListEntry( FilterListEntry( "info", false )  );
}


tarch::logging::CommandLineLogger& tarch::logging::CommandLineLogger::getInstance() {
  static CommandLineLogger singleton;
  return singleton;
}


std::ostream& tarch::logging::CommandLineLogger::out() {
  if (_outputStream==nullptr) {
    return std::cout;
  }
  else {
    return *_outputStream;
  }
}


void tarch::logging::CommandLineLogger::configureOutputStreams() {
  out().setf( std::ios_base::scientific, std::ios_base::floatfield );
  out().precision(20);
  std::cerr.setf( std::ios_base::scientific, std::ios_base::floatfield );
  std::cerr.precision(20);
}


tarch::logging::CommandLineLogger& tarch::logging::CommandLineLogger::operator=(const CommandLineLogger& rhs) {
  return *this;
}


tarch::logging::CommandLineLogger::CommandLineLogger(const CommandLineLogger& param) {
}


tarch::logging::CommandLineLogger::~CommandLineLogger() {
  if (_outputStream!=nullptr) {
    delete _outputStream;
    _outputStream = nullptr;
  }
}


std::string tarch::logging::CommandLineLogger::addSeparators(std::string::size_type spaces, std::string message) const {
  if ( message.size() > 0 ) {
    while (message.size() < spaces) {
      message += " ";
    }
    message = getLogColumnSeparator() + message;
  }

  return message;
}


std::string tarch::logging::CommandLineLogger::constructMessageString(
  std::string          messageType,
  const long int&      timestampMS,
  std::string          timestampHumanReadable,
  std::string          machineName,
  std::string          trace,
  const std::string&   message
) {
  std::string prefix = "";
  for (unsigned int i=0; i<_indent; i++ ) prefix += " ";

  std::string result;

  if ( getLogTimeStamp() && timestampMS!=0) {
    std::ostringstream timeStampString;
    timeStampString << timestampMS;
    result += addSeparators(NumberOfStandardColumnSpaces,timeStampString.str() );
  }

  if ( getLogTimeStampHumanReadable() ) {
    result += addSeparators(NumberOfStandardColumnSpaces,timestampHumanReadable);
  }

  if ( getLogMachineName() ) {
    result += addSeparators(NumberOfStandardColumnSpaces,machineName);
  }

  if ( getLogMessageType() ) {
    result += addSeparators(NumberOfStandardColumnSpaces,messageType);
  }

  if ( getLogTrace() ) {
    result += addSeparators(NumberOfTraceColumnSpaces,trace);
  }

  result += addSeparators(NumberOfStandardColumnSpaces,prefix + message);

  result += "\n";

  return result;
}


void tarch::logging::CommandLineLogger::debug(const long int& timestampMS, const std::string& timestampHumanReadable, const std::string& machineName, const std::string& trace, const std::string& message) {
  if (writeDebug(trace)) {
    #if !defined(Debug)
    assertion(false);
    #endif

    std::string outputMessage = constructMessageString(
      "debug",
      timestampMS,
      timestampHumanReadable,
      machineName,
      trace,
      message
    );

    tarch::multicore::Lock lockCout( _semaphore );
    out() << outputMessage;
    out().flush();
  }
}


void tarch::logging::CommandLineLogger::info(const long int& timestampMS, const std::string& timestampHumanReadable, const std::string& machineName, const std::string& trace, const std::string& message) {
  if (writeInfo(trace)) {
    std::string outputMessage = constructMessageString(
      "info",
      timestampMS,
      timestampHumanReadable,
      machineName,
      trace,
      message
    );

    tarch::multicore::Lock lockCout( _semaphore );
    out() << outputMessage;
    if (&out()!=&std::cout) {
      std::cout << outputMessage;
    }
  }
}


void tarch::logging::CommandLineLogger::warning(const long int& timestampMS, const std::string& timestampHumanReadable, const std::string& machineName, const std::string& trace, const std::string& message) {
  
  if (writeWarning(trace)) {
    std::string outputMessage = constructMessageString(
      "warning",
      timestampMS,
      timestampHumanReadable,
      machineName,
      trace,
      message
    );

    tarch::multicore::Lock lockCout( _semaphore );
    out().flush();
    #ifdef CompilerCLX
    if(&out()!=&std::cout) {
      std::cout << outputMessage;
      std::cout.flush();
    }
    #else
    std::cerr << outputMessage;
    std::cerr.flush();
    #endif
  }
}


void tarch::logging::CommandLineLogger::error(const long int& timestampMS, const std::string& timestampHumanReadable, const std::string& machineName, const std::string& trace, const std::string& message) {
  if ( writeError(trace) ) {
    std::string outputMessage = constructMessageString(
      "error",
      timestampMS,
      timestampHumanReadable,
      machineName,
      trace,
      message
    );

    tarch::multicore::Lock lockCout( _semaphore );
    out().flush();
    #ifdef CompilerCLX
    if(&out()!=&std::cout) {
      std::cout << outputMessage;
      std::cout.flush();
    }
    #else
    std::cerr << outputMessage;
    std::cerr.flush();
    #endif
  }
}




void tarch::logging::CommandLineLogger::indent( bool indent, const std::string& trace, const std::string& message ) {
  #ifdef Debug

  tarch::multicore::Lock lockCout( _semaphore );
  if (indent) {
    _indent+=NumberOfIndentSpaces;
     #if !defined(SharedMemoryParallelisation)
    _indentTraces.push(trace);
     #endif
  }
  else {
    #if !defined(SharedMemoryParallelisation)
    assertionEquals2(
      _indentTraces.top(),
      trace,
      message,
      indent
    );
    _indentTraces.pop();
    #endif
    assertion5(
      _indent >= NumberOfIndentSpaces,
      _indent, NumberOfIndentSpaces,
      "more logTraceOut calls than logTraceIn calls invoked before",
      trace, message
    );
    _indent-=NumberOfIndentSpaces;
  }
  #endif
}


void tarch::logging::CommandLineLogger::setLogFormat(const tarch::logging::configurations::LogOutputFormatConfiguration& configuration) {
  setLogFormat(
    configuration.getLogColumnSeparator(),
    configuration.getLogTimeStamp(),
    configuration.getLogTimeStampHumanReadable(),
    configuration.getLogMachineName(),
    configuration.getLogMessageType(),
    configuration.getLogTrace(),
    configuration.getLogFileName()
  );
}


void tarch::logging::CommandLineLogger::setLogFormat(
  const std::string& columnSeparator,
  bool logTimeStamp,
  bool logTimeStampHumanReadable,
  bool logMachineName,
  bool logMessageType,
  bool logTrace,
  const std::string&  outputLogFileName
) {
  _logColumnSeparator        = columnSeparator;
  _logTimeStamp              = logTimeStamp;
  _logTimeStampHumanReadable = logTimeStampHumanReadable;
  _logMachineName            = logMachineName;
  _logMessageType            = logMessageType;
  _logTrace                  = logTrace;

  if (!outputLogFileName.empty()) {
    _outputStream = new std::ofstream( outputLogFileName.c_str() );
  }
}


void tarch::logging::CommandLineLogger::addFilterListEntry( const FilterListEntry& entry) {
  if (_filterlist.count(entry)!=0) {
    logError( "addFilterListEntry(...)", "tried to insert " << entry.toString() << " multiple times");
  }
  else {
    _filterlist.insert( entry );
  }
}


void tarch::logging::CommandLineLogger::addFilterListEntries( const FilterList&    entries) {
  for (FilterList::const_iterator p = entries.begin(); p!=entries.end(); p++ ) {
    addFilterListEntry(*p);
  }
}


void tarch::logging::CommandLineLogger::clearFilterList() {
  _filterlist.clear();
}


bool tarch::logging::CommandLineLogger::filterOut(
  const std::string& targetName,
  const std::string& className,
  int                rank
) {
  bool result    = true;
  bool foundRule = false;
  int lengthActive = 0;
  for (FilterList::const_iterator p = _filterlist.begin(); p!=_filterlist.end(); p++ ) {
    int length = static_cast<int>(p->_namespaceName.size());
    if ( length >= lengthActive ) {
      #ifdef Parallel
      if (
        (targetName.find(p->_targetName, 0)==0) &&
        (className.find( p->_namespaceName, 0)==0) &&
        (p->_rank == -1 || p->_rank == tarch::parallel::Node::getInstance().getRank())
      ) {
        lengthActive = length;
        result       = p->_isBlackEntry;
        foundRule    = true;
      }
      #else
      if ( (targetName.find(p->_targetName, 0)==0) &&
           (className.find( p->_namespaceName, 0 ) == 0))
      {
        lengthActive = length;
        result       = p->_isBlackEntry;
        foundRule    = true;
      }
      #endif
    }
  }
  if (!foundRule) {
    logWarning( "filterOut(...)", "did not find filter rule for target \"" << targetName << "\" and class \"" << className << "\" on rank " << rank );
  }
  return result;
}




bool tarch::logging::CommandLineLogger::writeDebug(const std::string& className) {
  #ifdef Parallel
  int rank = -1;
  if (tarch::parallel::Node::getInstance().isInitialised()) {
    rank = tarch::parallel::Node::getInstance().getRank();
  }
  #else
  int rank = -1;
  #endif
  return !filterOut("debug",className,rank);
}


bool tarch::logging::CommandLineLogger::writeInfo(const std::string& className) {
  #ifdef Parallel
  int rank = -1;
  if (tarch::parallel::Node::getInstance().isInitialised()) {
    rank = tarch::parallel::Node::getInstance().getRank();
  }
  #else
  int rank = -1;
  #endif
  return !filterOut("info",className,rank);
}


bool tarch::logging::CommandLineLogger::writeWarning(const std::string& className) {
  return true;
}


bool tarch::logging::CommandLineLogger::writeError(const std::string& className) {
  return true;
}


std::string tarch::logging::CommandLineLogger::getLogColumnSeparator() const {
  return _logColumnSeparator;
}


bool tarch::logging::CommandLineLogger::getLogTimeStamp() const {
  return _logTimeStamp;
}


bool tarch::logging::CommandLineLogger::getLogTimeStampHumanReadable() const {
  return _logTimeStampHumanReadable;
}


bool tarch::logging::CommandLineLogger::getLogMachineName() const {
  return _logMachineName;
}


bool tarch::logging::CommandLineLogger::getLogMessageType() const {
  return _logMessageType;
}


bool tarch::logging::CommandLineLogger::getLogTrace() const {
  return _logTrace;
}


void tarch::logging::CommandLineLogger::setLogColumnSeparator( const std::string& separator ) {
  _logColumnSeparator = separator;
}


void tarch::logging::CommandLineLogger::setLogTimeStamp( bool value ) {
  _logTimeStamp = value;
}


void tarch::logging::CommandLineLogger::setLogTimeStampHumanReadable( bool value ) {
  _logTimeStampHumanReadable = value;
}


void tarch::logging::CommandLineLogger::setLogMachineName( bool value ) {
  _logMachineName = value;
}


void tarch::logging::CommandLineLogger::setLogMessageType( bool value ) {
  _logMessageType  = value;
}


void tarch::logging::CommandLineLogger::setLogTrace( bool value ) {
  _logTrace = value;
}


void tarch::logging::CommandLineLogger::receiveDanglingMessages() {
}


void tarch::logging::CommandLineLogger::printFilterListToWarningDevice() const {
  if (_filterlist.empty()) {
    logWarning( "printFilterListToWarningDevice()", "filter list is empty" );
  }
  else {
    logWarning( "printFilterListToWarningDevice()", "filter list is not empty and contains " << _filterlist.size() << " entries" );
    for (FilterList::const_iterator p = _filterlist.begin(); p!=_filterlist.end(); p++ ) {
      logWarning( "printFilterListToWarningDevice()", p->toString() );
    }
  }
}

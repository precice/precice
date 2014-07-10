// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LOGGING_COMMANDLINELOGGER_H_
#define _TARCH_LOGGING_COMMANDLINELOGGER_H_


#ifdef Parallel
#include <mpi.h>
#endif
#include <iostream>
#include <fstream>
#include <set>
#include <stack>

#include "tarch/logging/configurations/LogOutputFormatConfiguration.h"
#include "tarch/services/Service.h"

#include "tarch/multicore/BooleanSemaphore.h"
#include "tarch/multicore/MulticoreDefinitions.h"


namespace tarch {
  namespace logging {
    class CommandLineLogger;
  }
}

/**
 * Command Line Logger
 *
 * Standard log output device. Implements the LogHandle. Usually there's only
 * one instance of this class for only the Log type aggregates it in a static
 * way. For a detailed description of this compiler flags please browse through
 * the static consts of the class or the howto page.
 *
 * An interesting issue is the behaviour of the command line logger in the
 * parallel mode if debug is switched on. In this case the number of messages
 * soon becomes unreadable for any human. Thus, the command line logger
 * internally creates one output file rank.log for each rank and writes all the
 * information into this file. info, warnings and errors are written (in short form)
 * to std::cout / std::cerr, too.
 *
 * If the debug mode is not set, all information is written to std::cout.
 *
 * Note, that if multithreading is switched on (TBB-flag), the access to std::cout or the
 * logfile respectively is synchronized. The other methods are not thread-safe and should
 * be called only in safe context for defined behaviour.
 *
 * @version $Revision: 1.19 $
 * @author  Tobias Weinzierl, Wolfgang Eckhardt
 */
class tarch::logging::CommandLineLogger: public tarch::services::Service {
  public:
    /**
     * Represents one entry of the filter list. Syntax:
     *
     * <log-filter target="debug" component="component-name" rank="int|*" switch="on|off" />
     */
    struct FilterListEntry {
      /**
       * The message type target can be either "debug" or "info". Only messages
       * with corresponding target are filtered by the FilterListEntry.
       *
       * If left blank, all message types are targeted.
       */
      std::string _targetName;

      /**
       * Sometimes, one wants to block all log entries of one namespace of one
       * node (parallel case). In this case _rank holds the corresponding entry,
       * if all entries from one namespace have to be blocked, the value of
       * _rank equals -1.
       */
      int         _rank;

      /**
       * Name of the namespace that should not be logged.
       */
      std::string _namespaceName;

      /**
       * If true, filter list entry is a filter list entry, otherwise white list.
       */
      bool        _isBlackEntry;

      bool operator<(const FilterListEntry& b) const;
      bool operator==(const FilterListEntry& b) const;

      /**
       * Construct filter list entry for one target without any
       */
      FilterListEntry( const std::string& targetName="", bool isBlackListEntry=false );
      FilterListEntry( const std::string& targetName, int rank, const std::string& className, bool isBlackListEntry );

      std::string toString() const;
    };

    typedef std::set<FilterListEntry> FilterList;

    #if !defined(SharedMemoryParallelisation) && defined(Debug)
    std::stack<std::string>  _indentTraces;
    #endif
  private:
    static Log _log;

    FilterList  _filterlist;

    tarch::multicore::BooleanSemaphore _semaphore;

    /**
     * May not be const as it might write a warning itself
     */
    bool filterOut(
      const std::string& targetName,
      const std::string& className,
      int                rank
    );

    /**
     * Test for the column separator of a string output.
     */
    std::string    _logColumnSeparator;
    bool           _logTimeStamp;
    bool           _logTimeStampHumanReadable;
    bool           _logMachineName;
    bool           _logMessageType;
    bool           _logTrace;
    std::ostream*  _outputStream;


    /**
     * Declared private since assignment does not make sense for an output
     * class (output information mismatch).
     */
    CommandLineLogger& operator=(const CommandLineLogger& rhs);

    /**
     * Declared private since copying does not make sense for an output
     * class (output information mismatch).
     */
    CommandLineLogger(const CommandLineLogger& param);

    /**
     * Indent is supported only in debug mode.
     */
    static std::string::size_type _indent;

    /**
     * @todo
     */
    static const std::string::size_type NumberOfIndentSpaces;

    /**
     * @todo
     */
    static const std::string::size_type NumberOfStandardColumnSpaces;

    /**
     * The trace column is not formatted using only tabulators, but it uses
     * spaces to create a unique column size: First, TRACE_SPACE_NUMBER-n
     * spaces are added if n is the length of the original string, then a \\t
     * is appended. This message augmentation is done only, if the Debug-mode
     * is set.
     */
    static const std::string::size_type NumberOfTraceColumnSpaces;

    /**
     * Takes the message and adds spaces such that the entries are aligned like
     * in a table. Should either be passed NumberOfIndentSpaces or
     * NumberOfTraceSpaces.
     */
    std::string addSeparators(std::string::size_type spaces, std::string message) const;

    /**
     * Construct message string
     *
     * !!! Thread Safety
     *
     * The message string relies on the global field _indent. This one might
     * change throughout the execution of this method. However, I accept such a
     * behavior: Changing _indent throughout the message execution makes the
     * method add the wrong number of whitespaces in front of the message. That
     * is a 'bug' we can accept.
     */
    std::string constructMessageString(
      std::string          messageType,
      const long int&      timestampMS,
      std::string          timestampHumanReadable,
      std::string          machineName,
      std::string          trace,
      const std::string&   message
    );

    /**
     * Configures the output streams
     */
    void configureOutputStreams();

    bool writeDebug(const std::string& trace);
    bool writeInfo(const std::string& trace);
    bool writeWarning(const std::string& trace);
    bool writeError(const std::string& trace);

    /**
     * It's a singleton.
     */
    CommandLineLogger();

    std::string getLogColumnSeparator() const;
    bool        getLogTimeStamp() const;
    bool        getLogTimeStampHumanReadable() const;
    bool        getLogMachineName() const;
    bool        getLogMessageType() const;
    bool        getLogTrace() const;

    std::ostream& out();
  public:
    virtual ~CommandLineLogger();

    static CommandLineLogger& getInstance();

    /**
     * Add one filter list entry
     *
     * If you wanna switch on the logging globally, please add
     *
     * tarch::logging::CommandLineLogger::getInstance().addFilterListEntry(tarch::logging::CommandLineLogger::FilterListEntry());
     *
     * to your configuration.
     */
    void addFilterListEntry( const FilterListEntry& entry);
    void addFilterListEntries( const FilterList&    entries);
    void clearFilterList();

    void debug(      const long int& timestampMS, const std::string& timestampHumanReadable, const std::string& machineName, const std::string& trace, const std::string& message);
    void info(       const long int& timestampMS, const std::string& timestampHumanReadable, const std::string& machineName, const std::string& trace, const std::string& message);

    /**
     * Write Warning
     *
     * In the implementation, I call a flush on cout before I write to cerr.
     * Otherwise, the cerr messages might overtake cout. Before the operation
     * returns, it does a flush on cerr, too. Otherwise, the message might not
     * occur, i.e. the application might shut down before the message is flushed
     * to the terminal.
     */
    void warning(    const long int& timestampMS, const std::string& timestampHumanReadable, const std::string& machineName, const std::string& trace, const std::string& message);

    /**
     * Write Error
     *
     * In the implementation, I call a flush on cout before I write to cerr.
     * Otherwise, the cerr messages might overtake cout. Before the operation
     * returns, it does a flush on cerr, too. Otherwise, the message might not
     * occur, i.e. the application might shut down before the message is flushed
     * to the terminal.
     */
    void error(      const long int& timestampMS, const std::string& timestampHumanReadable, const std::string& machineName, const std::string& trace, const std::string& message);

    /**
     * Tells the logger to increment/decrement the indent.
     *
     * !!! Thread Safety
     *
     * _indent is a global static field shared by all threads. If we increment
     * or decrement it, this is first of all a read followed by a write.
     * Consequently data races could occur and the counter could become smaller
     * than zero. This ain't possible in the sequential code as each increment
     * is accompanied by a decrement. The following table illustrates the race:
     *
     * || value of _indent  || Thread 1 || Thread 2
     * |  2                 |  initial condition   |  initial condition
     * |                    |  enter indent(false) |  enter indent(true)
     * |                    |  fetch indent into register |  fetch indent into register
     * |                    |  register value -= 2 |  register value += 2
     * |  4                 |  is a little bit slower |  write back new value of indent
     * |  0                 |  write back new value of indent |
     *
     * To avoid this data race, I introduced a semaphore. This one could also
     * be implemented with TBB's atomic construct, e.g., but I prefer the
     * semaphor / critical section technique.
     *
     * @param trace    Needed in debug mode to be able to find out who called indent(false) without an indent(true)
     * @param message  Needed in debug mode to be able to find out who called indent(false) without an indent(true)
     */
    void indent( bool indent, const std::string& trace, const std::string& message );

    void setLogColumnSeparator( const std::string& separator = " ");
    void setLogTimeStamp( bool value = true );
    void setLogTimeStampHumanReadable( bool value = true );
    void setLogMachineName( bool value = true );
    void setLogMessageType( bool value = true );
    void setLogTrace( bool value = true );

    void setLogFormat(const tarch::logging::configurations::LogOutputFormatConfiguration& configuration);

    void setLogFormat(
      const std::string& columnSeparator,
      bool logTimeStamp,
      bool logTimeStampHumanReadable,
      bool logMachineName,
      bool logMessageType,
      bool logTrace,
      const std::string&  outputLogFileName
    );

    virtual void receiveDanglingMessages();

    void printFilterListToWarningDevice() const;
};

#endif

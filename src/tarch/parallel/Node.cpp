#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/parallel/Node.h"
#include "tarch/Assertions.h"
#include "tarch/services/ServiceRepository.h"

#include <sstream>
#include <cstdlib>

#include "tarch/compiler/CompilerSpecificSettings.h"

/**
 * For the machine name. If it doesn't work, switch it off in the file
 * CompilerSpecificSettings.h.
 */
#ifdef CompilerHasUTSName
#include <sys/utsname.h>
#endif


tarch::logging::Log tarch::parallel::Node::_log("tarch::parallel::Node");


bool tarch::parallel::Node::_initIsCalled = false;


int tarch::parallel::Node::reserveFreeTag(const std::string& fullQualifiedMessageName) {
  static int result = 0;
  result++;
  
  tarch::logging::Log _log("tarch::parallel::Node<static>");

  logDebug(
    "reserveFreeTag()",
    "assigned message " << fullQualifiedMessageName
     << " the free tag " << result
  );
  
  return result;
}


bool tarch::parallel::Node::isInitialised() const {
  return _initIsCalled;
}


void tarch::parallel::Node::ensureThatMessageQueuesAreEmpty( int fromRank, int tag ) {
  MPI_Status   status;
  int          flag;
  MPI_Iprobe(fromRank, tag, _communicator, &flag, &status);
  assertion3( flag==0, fromRank, tag, getRank() );
}


void tarch::parallel::Node::plotMessageQueues() {
  MPI_Status   status;
  int          flag;

  MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, _communicator, &flag, &status);
  if (flag==0) {
    _log.error("plotMessageQueues()", "there are no messages from any sender in MPI queue");
  }
  else {
    logError(
      "plotMessageQueues()",
      "there is still a message in queue "
      " from rank " << status.MPI_SOURCE <<
      " with tag " << status.MPI_TAG
    );
  }
}

void tarch::parallel::Node::triggerDeadlockTimeOut(
  const std::string&  className,
  const std::string&  methodName,
  int                 communicationPartnerRank,
  const std::string&  comment
) {
  _log.warning("triggerDeadlockTimeOut",
      "old function signature triggered - compatibility mode used. A function outputting the tag and the number of expected messages would be available");
  triggerDeadlockTimeOut(
      className,
      methodName,
      communicationPartnerRank,
      0,
      0,
      comment
      );
}

void tarch::parallel::Node::writeTimeOutWarning(
  const std::string&  className,
  const std::string&  methodName,
  int                 communicationPartnerRank
) {
  _log.warning("triggerDeadlockTimeOut",
      "old function signature triggered - compatibility mode used. A function outputting the tag and the number of expected messages would be available");
  writeTimeOutWarning(
      className,
      methodName,
      communicationPartnerRank,
      0,
      0
      );
}

void tarch::parallel::Node::triggerDeadlockTimeOut(
  const std::string&  className,
  const std::string&  methodName,
  int                 communicationPartnerRank,
  int                 tag,
  int                 numberOfExpectedMessages,
  const std::string&  comment
) {
  std::ostringstream out;
  out << "operation " << className << "::" << methodName << " on node "
      << getRank() << " had to wait more than " << _deadlockTimeOut
      << " seconds for " << numberOfExpectedMessages
      << " message(s) from node " << communicationPartnerRank << " with tag " << tag
      << ". Timeout. " << comment;
  _log.error( "triggerDeadlockTimeOut(...)", out.str() );

  plotMessageQueues();

  exit(DEADLOCK_EXIT_CODE);
}


void tarch::parallel::Node::writeTimeOutWarning(
  const std::string&  className,
  const std::string&  methodName,
  int                 communicationPartnerRank,
  int                 tag,
  int                 numberOfExpectedMessages
) {
  std::ostringstream out;
  out << "operation " << className << "::" << methodName << " on node "
      << getRank() << " had to wait more than " << _deadlockTimeOut
      << " seconds for " << numberOfExpectedMessages
      << " message(s) from node " << communicationPartnerRank << " with tag " << tag << ". Application "
      << "will terminate after " << _deadlockTimeOut << " seconds because "
      << "of a deadlock";
  _log.warning( "writeTimeOutWarning(...)", out.str() );
}


clock_t tarch::parallel::Node::getDeadlockWarningTimeStamp() const {
  clock_t result = clock() + _timeOutWarning * CLOCKS_PER_SEC;
  assertion4( result>=0, result, clock(), _timeOutWarning, CLOCKS_PER_SEC);

  return result;
}


clock_t tarch::parallel::Node::getDeadlockTimeOutTimeStamp() const {
  clock_t result = clock() + _deadlockTimeOut * CLOCKS_PER_SEC;
  assertion4( result>=0, result, clock(), _timeOutWarning, CLOCKS_PER_SEC);

  return result;
}


bool tarch::parallel::Node::isTimeOutDeadlockEnabled() const {
  return _deadlockTimeOut > 0;
}


bool tarch::parallel::Node::isTimeOutWarningEnabled() const {
  return _timeOutWarning > 0;
}


std::string tarch::parallel::MPIReturnValueToString( int result ) {
  std::ostringstream out;

  int   resultlen;
  char* string = new char[MPI_MAX_ERROR_STRING];  // (char *)malloc(MPI_MAX_ERROR_STRING * sizeof(char));
  MPI_Error_string(result, string, &resultlen);

  int   errorclass;
  MPI_Error_class(result, &errorclass);

  out << "mpi error class: " << errorclass << "="
      << ", mpi error text: " << string;

  switch ( errorclass ) {
    case MPI_SUCCESS:      out << "MPI_SUCCESS [no error]"; break;
    case MPI_ERR_BUFFER:   out << "MPI_ERR_BUFFER [invalid buffer pointer]"; break;
    case MPI_ERR_COUNT:    out << "MPI_ERR_COUNT [invalid count argument]"; break;
    case MPI_ERR_TYPE:     out << "MPI_ERR_TYPE [invalid datatype]"; break;
    case MPI_ERR_TAG:      out << "MPI_ERR_TAG [invalid tag]"; break;
    case MPI_ERR_COMM:     out << "MPI_ERR_COMM [invalid communicator]"; break;
    case MPI_ERR_RANK:     out << "MPI_ERR_RANK [invalid rank]"; break;
    case MPI_ERR_REQUEST:  out << "MPI_ERR_REQUEST [invalid request handle]"; break;
    case MPI_ERR_ROOT:     out << "MPI_ERR_ROOT [invalid root argument]"; break;
    case MPI_ERR_GROUP:    out << "MPI_ERR_GROUP [invalid group]"; break;
    case MPI_ERR_OP:       out << "MPI_ERR_OP [invalid operation]"; break;
    case MPI_ERR_TOPOLOGY: out << "MPI_ERR_TOPOLOGY [invalid topology]"; break;
    case MPI_ERR_DIMS:     out << "MPI_ERR_DIMS [invalid dimensions]"; break;
    case MPI_ERR_ARG:      out << "MPI_ERR_ARG [invalid argument]"; break;
    case MPI_ERR_UNKNOWN:  out << "MPI_ERR_UNKNOWN [unknown error]"; break;
    case MPI_ERR_TRUNCATE: out << "MPI_ERR_TRUNCATE [message has been truncated by receiver]"; break;
    case MPI_ERR_OTHER:    out << "MPI_ERR_OTHER [other unknown error]"; break;
    case MPI_ERR_INTERN:   out << "MPI_ERR_INTERN [internal mpi error]"; break;
    default: out << "unknown";
  }

  delete[] string;
  return out.str();
}


std::string tarch::parallel::MPIStatusToString( const MPI_Status& status ) {
  std::ostringstream out;
  out << "status flag:"
      << " MPI_ERROR=" << status.MPI_ERROR
      << " (" << MPIReturnValueToString(status.MPI_ERROR)
      << ") ,MPI_SOURCE=" << status.MPI_SOURCE
      << ",MPI_TAG=" << status.MPI_TAG;

  return out.str();
}


tarch::parallel::Node::Node():
  _rank(-1),
  _numberOfProcessors(-1),
  _communicator( MPI_COMM_WORLD),
  _timeOutWarning(0),
  _deadlockTimeOut(0) {
}


tarch::parallel::Node::Node(const parallel::Node& node):
  _rank(-1),
  _numberOfProcessors(-1),
  _communicator( MPI_COMM_WORLD) {
}


tarch::parallel::Node::~Node() {
}


void tarch::parallel::Node::shutdown() {
  assertion( _rank!=-1 );
  MPI_Barrier( _communicator );
  MPI_Finalize();
  _rank         = -1;
  _communicator = MPI_COMM_WORLD;
}


int tarch::parallel::Node::getGlobalMasterRank() {
  return 0;
}


bool tarch::parallel::Node::isGlobalMaster() const {
  assertion(_initIsCalled);
  return getRank() == getGlobalMasterRank();
}


void tarch::parallel::Node::logStatus() const {
  std::ostringstream statusMessage;
  statusMessage << "MPI status:";

  #ifdef CompilerHasUTSName
  utsname* utsdata = new utsname();
  assertion( utsdata!=NULL );
  uname(utsdata);
  statusMessage << " nodename=" << utsdata->nodename;
  delete utsdata;
  #else
  statusMessage << " nodename=undef";
  #endif

  statusMessage << ", rank=" << _rank;
  statusMessage << ", communicator=" << _communicator;
  statusMessage << ", #processors=" << _numberOfProcessors;

  _log.debug( "logStatus()", statusMessage.str() );
}


bool tarch::parallel::Node::init(int* argc, char*** argv) {
  int result = MPI_SUCCESS;

  result = MPI_Init( argc, argv );
  if (result!=MPI_SUCCESS) {
    std::cerr << "init(int*,char***)\t initialisation failed: " + MPIReturnValueToString(result) + " (no logging available yet)";
    return false;
  }

  result = MPI_Comm_size( MPI_COMM_WORLD, &_numberOfProcessors );
  if (result!=MPI_SUCCESS) {
    std::cerr << "init(int*,char***)\t initialisation failed: " + MPIReturnValueToString(result) + " (no logging available yet)";
    return false;
  }

  result = MPI_Comm_rank( MPI_COMM_WORLD, &_rank );
  if (result!=MPI_SUCCESS) {
    std::cerr << "init(int*,char***)\t initialisation failed: " + MPIReturnValueToString(result) + " (no logging available yet)";
    return false;
  }

  _initIsCalled = true;
  return true;
}


int tarch::parallel::Node::getRank() const {
  assertion(_initIsCalled);
  return _rank;
}


tarch::parallel::Node& tarch::parallel::Node::getInstance() {
  static Node singleton;
  return singleton;
}


MPI_Comm tarch::parallel::Node::getCommunicator() const {
  assertion(_initIsCalled);
  return _communicator;
}


int tarch::parallel::Node::getNumberOfNodes() const {
  assertion(_initIsCalled);
  return _numberOfProcessors;
}


void tarch::parallel::Node::setTimeOutWarning( const clock_t & value ) {
  assertion( value>=0 );
  _timeOutWarning = value;
}


void tarch::parallel::Node::setDeadlockTimeOut( const clock_t & value ) {
  assertion( value>=0 );
  _deadlockTimeOut = value;
}


void tarch::parallel::Node::setCommunicator( MPI_Comm communicator ) {
  _communicator = communicator;
}


void tarch::parallel::Node::receiveDanglingMessages() {
  tarch::services::ServiceRepository::getInstance().receiveDanglingMessages();
}

#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/parallel/NodePool.h"
#include "tarch/parallel/Node.h"


#include "tarch/parallel/messages/ActivationMessage.h"
#include "tarch/parallel/messages/JobRequestMessage.h"
#include "tarch/parallel/messages/NodePoolAnswerMessage.h"

#include <sstream>

#include "tarch/services/ServiceFactory.h"
registerService(tarch::parallel::NodePool)


tarch::logging::Log tarch::parallel::NodePool::_log("tarch::parallel::NodePool");


const int tarch::parallel::NodePool::NoFreeNodesMessage = -1;

const int tarch::parallel::NodePool::JobRequestMessageAnswerValues::NewMaster = 0;
const int tarch::parallel::NodePool::JobRequestMessageAnswerValues::Terminate = -1;
const int tarch::parallel::NodePool::JobRequestMessageAnswerValues::HandleLocalProblem = -2;
const int tarch::parallel::NodePool::JobRequestMessageAnswerValues::ExecuteGlobalTask = -3;


tarch::parallel::NodePool::NodePool():
  _masterNode(-1),
  _registrationTag(-1),
  _jobManagementTag(-1),
  _jobServicesTag(-1),
  _isAlive(false),
  _answerJobRequestsWithHandleLocalProblem(false),
  _strategy(0) {
  #ifdef Asserts
  _isInitialised = false;
  #endif
}


void tarch::parallel::NodePool::restart() {
  assertion1( !_isAlive, Node::getInstance().getRank() );
  assertion1( !Node::getInstance().isGlobalMaster() || _strategy!=0 , Node::getInstance().getRank() );
  assertion1( !Node::getInstance().isGlobalMaster() || _strategy->getNumberOfRegisteredNodes()==0, Node::getInstance().getRank() );

  logTraceIn( "restart()" );

  _isAlive = true;

  MPI_Barrier(Node::getInstance().getCommunicator());

  if ( !Node::getInstance().isGlobalMaster() ) {
    logDebug( "restart()", "start to register at node pool" );
    tarch::parallel::messages::RegisterAtNodePoolMessage registerMessage(
      tarch::parallel::StringTools::convert(_log.getMachineInformation())
    );
    registerMessage.send( Node::getGlobalMasterRank(), _registrationTag, true);
    logDebug( "restart()", "register message sent: " << registerMessage.toString() << " on tag " << _registrationTag );
  }
  logTraceOut( "restart()" );
}


int tarch::parallel::NodePool::getNumberOfWorkingNodes() const {
  assertion1( Node::getInstance().isGlobalMaster(), Node::getInstance().getRank() );
  assertion1( _strategy!=0, Node::getInstance().getRank() );

  return _strategy->getNumberOfRegisteredNodes() - _strategy->getNumberOfIdleNodes();
}


void tarch::parallel::NodePool::init() {
  _registrationTag  = Node::getInstance().reserveFreeTag( "tarch::parallel::NodePool[registration]" );
  _jobManagementTag = Node::getInstance().reserveFreeTag( "tarch::parallel::NodePool[job-management]" );
  _jobServicesTag   = Node::getInstance().reserveFreeTag( "tarch::parallel::NodePool[job-services]" );

  tarch::parallel::messages::ActivationMessage::initDatatype();
  tarch::parallel::messages::JobRequestMessage::initDatatype();
  tarch::parallel::messages::NodePoolAnswerMessage::initDatatype();
  tarch::parallel::messages::RegisterAtNodePoolMessage::initDatatype();
  tarch::parallel::messages::WorkerRequestMessage::initDatatype();
}


tarch::parallel::NodePool& tarch::parallel::NodePool::getInstance() {
  static tarch::parallel::NodePool pool;
  return pool;
}


int tarch::parallel::NodePool::getTagForForkMessages() const {
  return _jobManagementTag;
}


void tarch::parallel::NodePool::answerAllJobRequestMessagesWithHandleLocalProblem( bool value ) {
  _answerJobRequestsWithHandleLocalProblem = value;
}

void tarch::parallel::NodePool::sendExecuteGlobalTaskJobMessages() {
  assertion1( Node::getInstance().isGlobalMaster(), Node::getInstance().getRank() );

  _answerJobRequestsWithExecuteGlobalTask = true;

  tarch::parallel::NodePool::getInstance().receiveDanglingMessages();
  while(_strategy->hasIdleNode()) {
    int node = _strategy->reserveNode(Node::getInstance().getRank());

    tarch::parallel::messages::ActivationMessage answerMessage( JobRequestMessageAnswerValues::ExecuteGlobalTask );
    answerMessage.send( node, _jobManagementTag, true );
  }
}


void tarch::parallel::NodePool::setStrategy(NodePoolStrategy* strategy) {
  assertion1( Node::getInstance().isGlobalMaster(), Node::getInstance().getRank() );

  logTraceIn( "setStrategy(...)" );

  if (_strategy!=0) {
    assertion1( _strategy->getNumberOfRegisteredNodes()==0, Node::getInstance().getRank() );

    delete _strategy;
  }

  _strategy = strategy;
  _strategy->setNodePoolTag(_jobServicesTag);

  logTraceOut( "setStrategy(...)" );
}



void tarch::parallel::NodePool::waitForAllNodesToBecomeIdle() {
  assertion1( Node::getInstance().isGlobalMaster(), Node::getInstance().getRank() );
  assertion1( _strategy!=0, Node::getInstance().getRank() );

  while ( _strategy->getNumberOfIdleNodes() < Node::getInstance().getNumberOfNodes()-1) {
    replyToJobRequestMessages();
  }
}


void tarch::parallel::NodePool::shutdown() {
  emptyReceiveBuffers();

  if ( _isAlive ) {
    _log.error("shutdown()", "called destructor for alive node" );
  }

  if (_strategy!=0) {
    delete _strategy;
    _strategy = 0;
  }

  tarch::parallel::messages::ActivationMessage::shutdownDatatype();
  tarch::parallel::messages::JobRequestMessage::shutdownDatatype();
  tarch::parallel::messages::NodePoolAnswerMessage::shutdownDatatype();
  tarch::parallel::messages::RegisterAtNodePoolMessage::shutdownDatatype();
  tarch::parallel::messages::WorkerRequestMessage::shutdownDatatype();
}


tarch::parallel::NodePool::~NodePool() {
  if (_strategy != 0) {
    logError( "~NodePool()", "forgot to call shutdown() on node " << Node::getInstance().getRank() );
  }
}


tarch::parallel::NodePool::JobRequestMessageAnswer tarch::parallel::NodePool::waitForJob() {
  logTraceIn( "waitForJob()" );

  assertion1( _isAlive, Node::getInstance().getRank() );
  assertion( !Node::getInstance().isGlobalMaster() );

  _masterNode = -1;

  tarch::parallel::messages::JobRequestMessage message;
  message.send(Node::getInstance().getGlobalMasterRank(),_jobManagementTag, true);

  MPI_Status   status;
  tarch::parallel::messages::ActivationMessage answer;
  int result = MPI_Recv(
    &answer, 1,
    tarch::parallel::messages::ActivationMessage::Datatype,
    Node::getInstance().getGlobalMasterRank(), _jobManagementTag,
    tarch::parallel::Node::getInstance().getCommunicator(),
    &status
  );
  if ( result != MPI_SUCCESS ) {
    logError(
      "waitForJob()",
      "failed to start to receive tarch::parallel::messages::ActivationMessage from master node: "
      << tarch::parallel::MPIReturnValueToString(result)
    );
  }

  if ( answer.getNewMaster() == JobRequestMessageAnswerValues::Terminate ) {
  	logDebug("waitForJob()", "node received termination signal");
  	_isAlive = false;
    logTraceInWith1Argument( "waitForJob()", false );
  	return JobRequestMessageAnswerValues::Terminate;
  }
  // TODO: check if HandleLocalProblem should be handled here!
  /*else if ( answer.getNewMaster() == JobRequestMessageAnswerValues::HandleLocalProblem ) {
    return JobRequestMessageAnswerValues::HandleLocalProblem;
  }*/
  else if ( answer.getNewMaster() == JobRequestMessageAnswerValues::ExecuteGlobalTask ) {
    return JobRequestMessageAnswerValues::ExecuteGlobalTask;
  }
  else {
    _masterNode = answer.getNewMaster();
    logTraceOutWith2Arguments( "waitForJob()", true, _masterNode );
    assertion1(_masterNode>=0, _masterNode);
    return _masterNode;
  }
}


void tarch::parallel::NodePool::terminate() {
  assertion1( Node::getInstance().isGlobalMaster(), Node::getInstance().getRank() );
  assertion1( _strategy!=0, Node::getInstance().getRank() );
  assertion1( _isAlive, Node::getInstance().getRank() );

  logTraceIn("terminate()" );

  _isAlive = false;

  while ( _strategy->hasIdleNode() ) {
	  int rank = _strategy->removeNextIdleNode();
	  tarch::parallel::messages::ActivationMessage answerMessage( JobRequestMessageAnswerValues::Terminate );
    answerMessage.send( rank, _jobManagementTag, true );
  }

  logDebug( "terminate()", "still working: " << _strategy->getNumberOfRegisteredNodes() << " node(s)" );

  clock_t      timeOutWarning   = Node::getInstance().getDeadlockWarningTimeStamp();
  clock_t      timeOutShutdown  = Node::getInstance().getDeadlockTimeOutTimeStamp();
  bool         triggeredTimeoutWarning = false;

  assertion1( _strategy!=0, Node::getInstance().getRank() );
  while ( _strategy->getNumberOfRegisteredNodes()>0 ) {
    Node::getInstance().receiveDanglingMessages();

    // deadlock aspect
    if ( Node::getInstance().isTimeOutWarningEnabled() && (clock()>timeOutWarning) && (!triggeredTimeoutWarning)) {
      Node::getInstance().writeTimeOutWarning( "tarch::parallel::NodePool", "terminate()", -1, _jobManagementTag, 1);
      triggeredTimeoutWarning = true;
    }
    if ( Node::getInstance().isTimeOutDeadlockEnabled() && (clock()>timeOutShutdown)) {
      Node::getInstance().triggerDeadlockTimeOut( "tarch::parallel::NodePool", "terminate()", -1, _jobManagementTag, 1 );
    }
  }

  logTraceOut( "terminate()" );
}


int tarch::parallel::NodePool::getFreeNode(int forMaster) {
  assertion1( _isAlive, Node::getInstance().getRank() );
  assertion1( _strategy!=0, Node::getInstance().getRank() );

  logTraceInWith1Argument( "getFreeNode(int)", forMaster );

  int result;
  if ( _strategy->hasIdleNode() ) {
    result = _strategy->reserveNode(forMaster);
  }
  else {
    result = NoFreeNodesMessage;
  }

  logTraceOutWith1Argument( "getFreeNode(int)", result );
  return result;
}


int tarch::parallel::NodePool::reserveFreeNodeForServer() {
  assertion1( _isAlive, Node::getInstance().getRank() );

  logTraceIn( "reserveFreeNodeForServer()");

  int activatedNode = getFreeNode( Node::getInstance().getGlobalMasterRank());

  if (activatedNode!=NoFreeNodesMessage) {
    tarch::parallel::messages::ActivationMessage message( Node::getInstance().getGlobalMasterRank() );
    message.send( activatedNode, _jobManagementTag, true );
  }

  logTraceOutWith1Argument( "reserveFreeNodeForServer()", activatedNode);

  return activatedNode;
}


int tarch::parallel::NodePool::reserveFreeNode() {
  if ( Node::getInstance().isGlobalMaster() ) {
    assertion2( _masterNode == -1, _masterNode, Node::getInstance().getRank() );
  	assertion1( _isAlive, Node::getInstance().getRank() );
  	receiveDanglingMessages();
  	return reserveFreeNodeForServer();
  }
  else {
    assertion1( _isAlive, Node::getInstance().getRank() );
    return reserveFreeNodeForClient();
  }
}


int tarch::parallel::NodePool::reserveFreeNodeForClient() {
  assertion1( _isAlive, Node::getInstance().getRank() );
  assertion1( !Node::getInstance().isGlobalMaster(), Node::getInstance().getRank() );

  logTraceIn( "reserveFreeNodeForClient()" );

  tarch::parallel::messages::WorkerRequestMessage queryMessage;
  queryMessage.send(Node::getInstance().getGlobalMasterRank(),_jobServicesTag, true);

  tarch::parallel::messages::NodePoolAnswerMessage answer;
  answer.receive(Node::getInstance().getGlobalMasterRank(),_jobServicesTag, true);

  logTraceOutWith1Argument( "reserveFreeNodeForClient()", answer.getNewWorker() );

  return answer.getNewWorker();
}


tarch::parallel::NodePool::NodePool( const NodePool& pool ) {
  assertionMsg( false, "copy not allowed" );
}


tarch::parallel::NodePool& tarch::parallel::NodePool::operator=( const tarch::parallel::NodePool& pool ) {
  assertionMsg( false, "copy not allowed" );
  return *this;
}


void tarch::parallel::NodePool::receiveDanglingMessages() {
  if ( Node::getInstance().isGlobalMaster() ) {
    replyToRegistrationMessages();
    replyToJobRequestMessages();
    replyToWorkerRequestMessages();
  }
}


int tarch::parallel::NodePool::getMasterRank() const {
  assertion1WithExplanation(
    !Node::getInstance().isGlobalMaster(),
    Node::getInstance().getRank(),
    "You may not call getMasterRank() on the global master (typically rank 0). \nUse isGlobalMaster() to check before whether operation may be called."
  );
  return _masterNode;
}


void tarch::parallel::NodePool::emptyReceiveBuffers() {
  emptyRegisterMessageReceiveBuffer();
  emptyJobRequestMessageBuffer();
  emptyWorkerRequestMessageBuffer();
}


void tarch::parallel::NodePool::replyToRegistrationMessages() {
  assertion1( _strategy!=0, Node::getInstance().getRank() );

  logTraceInWith1Argument( "replyToRegistrationMessages()", tarch::parallel::messages::RegisterAtNodePoolMessage::isMessageInQueue(_registrationTag, true) );

  while ( tarch::parallel::messages::RegisterAtNodePoolMessage::isMessageInQueue(_registrationTag, true) ) {
    tarch::parallel::messages::RegisterAtNodePoolMessage message;
    message.receive( MPI_ANY_SOURCE, _registrationTag, true );
    logDebug(  "replyToRegistrationMessages()", "got registration from rank " << message.getSenderRank() );
    _strategy->addNode( message );
  }

  logTraceOut( "replyToRegistrationMessages()" );
}


void tarch::parallel::NodePool::replyToJobRequestMessages() {
  assertion1( _strategy!=0, Node::getInstance().getRank() );

  while ( tarch::parallel::messages::JobRequestMessage::isMessageInQueue(_jobManagementTag, true) ) {
    tarch::parallel::messages::JobRequestMessage queryMessage;
    queryMessage.receive( MPI_ANY_SOURCE, _jobManagementTag, true );

    assertion1( queryMessage.getSenderRank() !=Node::getInstance().getGlobalMasterRank(), Node::getInstance().getRank() );

    if ( !_strategy->isRegisteredNode(queryMessage.getSenderRank()) ) {
      logDebug(
        "replyToJobRequestMessages()",
        "node pool does not contain entry for rank " << queryMessage.getSenderRank()
         << ". Message from rank " << queryMessage.getSenderRank()
         << " might have overtaken registration message. Waiting for registration"
      );

      while ( !_strategy->isRegisteredNode(queryMessage.getSenderRank()) ) {
        replyToRegistrationMessages();
      }

      logDebug( "replyToJobRequestMessages()", "registration from " << queryMessage.getSenderRank() << " finally arrived" );
    }

    if ( !_isAlive ) {
      _strategy->setNodeIdle( queryMessage.getSenderRank() );
      int rank = _strategy->removeNextIdleNode();
      assertionEquals1( rank, queryMessage.getSenderRank(), Node::getInstance().getRank() );
      tarch::parallel::messages::ActivationMessage answerMessage( JobRequestMessageAnswerValues::Terminate );
      answerMessage.send( rank, _jobManagementTag, true );
    }
    else if (_answerJobRequestsWithExecuteGlobalTask) {
      tarch::parallel::messages::ActivationMessage answerMessage( JobRequestMessageAnswerValues::ExecuteGlobalTask );
      answerMessage.send( queryMessage.getSenderRank(), _jobManagementTag, true );
    }
    else if (_answerJobRequestsWithHandleLocalProblem) {
      tarch::parallel::messages::ActivationMessage answerMessage( JobRequestMessageAnswerValues::HandleLocalProblem );
      answerMessage.send( queryMessage.getSenderRank(), _jobManagementTag, true );
    }
    else {
      _strategy->setNodeIdle( queryMessage.getSenderRank() );
    }
  }

  // maybe ExecuteGlobalTask should be done in a different way, could lead to deadlocks
  if(_answerJobRequestsWithExecuteGlobalTask) {
    _answerJobRequestsWithExecuteGlobalTask = false;
  }
}


void tarch::parallel::NodePool::replyToWorkerRequestMessages() {
  assertion1( _strategy!=0 , Node::getInstance().getRank() );

  logTraceInWith1Argument( "replyToWorkerRequestMessages()", tarch::parallel::messages::WorkerRequestMessage::isMessageInQueue(_jobServicesTag, true) );

  static NodePoolStrategy::RequestQueue queue;

  if (tarch::parallel::messages::WorkerRequestMessage::isMessageInQueue(_jobServicesTag, true)) {
    //NOTE: Take care of recursive calls.
    _strategy->fillWorkerRequestQueue(queue);

    while ( !queue.empty() ) {
      tarch::parallel::messages::WorkerRequestMessage nextRequestToAnswer = _strategy->extractElementFromRequestQueue(queue);
      if ( _isAlive &&  _strategy->hasIdleNode() ) {
        int activatedNode = _strategy->reserveNode(nextRequestToAnswer.getSenderRank());

        tarch::parallel::messages::NodePoolAnswerMessage answerMessage( activatedNode );
        answerMessage.send( nextRequestToAnswer.getSenderRank(), _jobServicesTag, true );

        tarch::parallel::messages::ActivationMessage activationMessage( nextRequestToAnswer.getSenderRank() );
        activationMessage.send( activatedNode, _jobManagementTag, true );
      }
      else {
        tarch::parallel::messages::NodePoolAnswerMessage answerMessage( NoFreeNodesMessage );
        answerMessage.send( nextRequestToAnswer.getSenderRank(), _jobServicesTag, true );
      }
    }
  }

  logTraceOut( "replyToWorkerRequestMessages()" );
}


void tarch::parallel::NodePool::emptyRegisterMessageReceiveBuffer() {
  while ( tarch::parallel::messages::RegisterAtNodePoolMessage::isMessageInQueue(_registrationTag, true) ) {
    tarch::parallel::messages::RegisterAtNodePoolMessage message;
    message.receive( MPI_ANY_SOURCE, _registrationTag, true );
  }
}


void tarch::parallel::NodePool::emptyJobRequestMessageBuffer() {
  while ( tarch::parallel::messages::JobRequestMessage::isMessageInQueue(_jobManagementTag, true) ) {
    tarch::parallel::messages::JobRequestMessage message;
    message.receive( MPI_ANY_SOURCE, _jobManagementTag, true );
  }
}


void tarch::parallel::NodePool::emptyWorkerRequestMessageBuffer() {
  while ( tarch::parallel::messages::WorkerRequestMessage::isMessageInQueue(_jobServicesTag, true) ) {
    tarch::parallel::messages::WorkerRequestMessage message;
    message.receive( MPI_ANY_SOURCE, _jobServicesTag, true );
  }
}

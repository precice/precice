#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/parallel/strategy/FCFSNodePoolStrategy.h"


tarch::logging::Log tarch::parallel::strategy::FCFSNodePoolStrategy::_log( "tarch::parallel::strategy::FCFSNodePoolStrategy" );


tarch::parallel::strategy::FCFSNodePoolStrategy::FCFSNodePoolStrategy():
  NodePoolStrategy(),
  _tag(-1),
  _nodes() {
}


tarch::parallel::strategy::FCFSNodePoolStrategy::~FCFSNodePoolStrategy() {
}


void tarch::parallel::strategy::FCFSNodePoolStrategy::fillWorkerRequestQueue(RequestQueue& queue) {
  assertion( _tag >= 0 );
  while ( tarch::parallel::messages::WorkerRequestMessage::isMessageInQueue(_tag, true) ) {
    tarch::parallel::messages::WorkerRequestMessage message;
    message.receive(MPI_ANY_SOURCE,_tag, true);
    queue.push_back( message );
  }
}


void tarch::parallel::strategy::FCFSNodePoolStrategy::logQueue( const RequestQueue& queue ) const {
  if (queue.empty()) {
	_log.debug( "logQueue()", "queue is empty" );
  }
  else {
    std::ostringstream msg;
    msg << "queue: ";

	for (RequestQueue::const_iterator p = queue.begin(); p != queue.end(); p++ ) {
	   msg << p->getSenderRank() << ",";
	}
	_log.debug( "logQueue()", msg.str() );
  }
}


tarch::parallel::messages::WorkerRequestMessage tarch::parallel::strategy::FCFSNodePoolStrategy::extractElementFromRequestQueue(RequestQueue& queue) {
  assertion( !queue.empty() );
  tarch::parallel::messages::WorkerRequestMessage result = queue.front();
  queue.pop_front();
  return result;
}


void tarch::parallel::strategy::FCFSNodePoolStrategy::addNode(const tarch::parallel::messages::RegisterAtNodePoolMessage& node) {
  assertion( !isRegisteredNode(node.getSenderRank()) );

  logTraceInWith1Argument( "addNode(...)", node.getSenderRank() );
  NodePoolListEntry newEntry(
    node.getSenderRank(),
    tarch::parallel::StringTools::convert(node.getNodeName())
  );
  _nodes.push_back( newEntry ) ;
  _nodes.sort();
  logTraceOutWith1Argument( "addNode(...)", newEntry.toString() );
}


void tarch::parallel::strategy::FCFSNodePoolStrategy::removeNode( int rank ) {
  assertion( isRegisteredNode(rank) );

  for (
    NodeContainer::iterator p = _nodes.begin();
    p != _nodes.end();
    p++
  ) {
    if ( p->getRank() == rank ) {
      #ifdef Debug
      _log.debug( "removeNode(int)", "remove entry " + p->toString() );
      #endif
      _nodes.erase(p);
      _nodes.sort();
      return;
    }
  }
}


bool tarch::parallel::strategy::FCFSNodePoolStrategy::hasIdleNode() const {
  return !_nodes.empty() &&
         _nodes.front().isIdle();
}


int tarch::parallel::strategy::FCFSNodePoolStrategy::removeNextIdleNode() {
  assertion( hasIdleNode() );
  int result = _nodes.front().getRank();
  _nodes.pop_front();
  return result;
}


int tarch::parallel::strategy::FCFSNodePoolStrategy::getNumberOfIdleNodes() const {
  int result = 0;
  NodeContainer::const_iterator p = _nodes.begin();
  while (p != _nodes.end()&& p->isIdle() ) {
	p++;
	result++;
  }
  return result;
}


void tarch::parallel::strategy::FCFSNodePoolStrategy::setNodeIdle( int rank ) {
  for (
    NodeContainer::iterator p = _nodes.begin();
    p != _nodes.end();
    p++
  ) {
    if ( p->getRank() == rank ) {
      p->deActivate();
    }
  }

  _nodes.sort();
}


bool tarch::parallel::strategy::FCFSNodePoolStrategy::isRegisteredNode(int rank) const {
  for (
    NodeContainer::const_iterator p = _nodes.begin();
    p != _nodes.end();
    p++
  ) {
    if ( p->getRank() == rank ) {
      return true;
    }
  }
  return false;
}


bool tarch::parallel::strategy::FCFSNodePoolStrategy::isIdleNode(int rank) const {
  assertion1( isRegisteredNode(rank), rank );
  for (
    NodeContainer::const_iterator p = _nodes.begin();
    p != _nodes.end();
    p++
  ) {
    if ( p->getRank() == rank && p->isIdle() ) {
      return true;
    }
  }
  return false;
}


void tarch::parallel::strategy::FCFSNodePoolStrategy::clearRegisteredNodes() {
  _nodes.clear();
}


int tarch::parallel::strategy::FCFSNodePoolStrategy::getNumberOfRegisteredNodes() const {
  return static_cast<int>( _nodes.size() );
}


std::string tarch::parallel::strategy::FCFSNodePoolStrategy::toString() const {
  std::ostringstream result;
  for (
    NodeContainer::const_iterator p = _nodes.begin();
    p != _nodes.end();
    p++
  ) {
	result << *p;
  }
  return result.str();
}


int tarch::parallel::strategy::FCFSNodePoolStrategy::reserveNode(int forMaster) {
  assertion(hasIdleNode());

  NodePoolListEntry result = _nodes.front();
  _nodes.pop_front();

  logDebug( "getFreeNode(int)", "found free node: " << result );

  result.activate();
  _nodes.push_back(result);
  _nodes.sort();

  return result.getRank();
}


void tarch::parallel::strategy::FCFSNodePoolStrategy::setNodePoolTag(int tag) {
  _tag = tag;
}

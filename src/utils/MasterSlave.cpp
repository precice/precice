//#ifndef PRECICE_NO_MPI

#include "MasterSlave.hpp"

#include "utils/assertion.hpp"
#include "com/Communication.hpp"

namespace precice {
namespace utils {

int MasterSlave::_rank = -1;
int MasterSlave::_size = -1;
bool MasterSlave::_isMaster = false;
bool MasterSlave::_isSlave = false;
com::PtrCommunication MasterSlave::_communication;


logging::Logger MasterSlave:: _log("utils::MasterSlave" );

void MasterSlave:: configure(int rank, int size)
{
  P_TRACE(rank, size);
  P_CHECK(size>=2, "You cannot use a master with a serial participant.");
  _rank = rank;
  _size = size;
  P_ASSERT(_rank != -1 && _size != -1);
  _isMaster = (rank==0);
  _isSlave = (rank!=0);
  P_DEBUG("isSlave: " << _isSlave <<", isMaster: " << _isMaster);
}

int MasterSlave::getRank()
{
  return _rank;
}

int MasterSlave::getSize()
{
  return _size;
}

bool MasterSlave::isMaster()
{
  return _isMaster;
}

bool MasterSlave::isSlave()
{
  return _isSlave;
}


double MasterSlave:: l2norm(const Eigen::VectorXd& vec)
{
  P_TRACE();

  if(not _isMaster && not _isSlave){ //old case
    return vec.norm();
  }

  P_ASSERT(_communication.get() != nullptr);
  P_ASSERT(_communication->isConnected());
  double localSum2 = 0.0;
  double globalSum2 = 0.0;

  for(int i=0; i<vec.size(); i++){
    localSum2 += vec(i)*vec(i);
  }

  // localSum is modified, do not use afterwards
  allreduceSum(localSum2, globalSum2, 1);
   /* old loop over all slaves solution
  if(_isSlave){
    _communication->send(localSum2, 0);
    _communication->receive(globalSum2, 0);
  }
  if(_isMaster){
    globalSum2 += localSum2;
    for(int rankSlave = 1; rankSlave < _size; rankSlave++){
      _communication->receive(localSum2, rankSlave);
      globalSum2 += localSum2;
    }
    for(int rankSlave = 1; rankSlave < _size; rankSlave++){
      _communication->send(globalSum2, rankSlave);
    }
  }
  */
  return sqrt(globalSum2);
}


double MasterSlave:: dot(const Eigen::VectorXd& vec1, const Eigen::VectorXd& vec2)
{
  P_TRACE();

  if(not _isMaster && not _isSlave){ //old case
    return vec1.dot(vec2);
  }

  P_ASSERT(_communication.get() != nullptr);
  P_ASSERT(_communication->isConnected());
  P_ASSERT(vec1.size()==vec2.size(), vec1.size(), vec2.size());
  double localSum = 0.0;
  double globalSum = 0.0;

  for(int i=0; i<vec1.size(); i++){
    localSum += vec1(i)*vec2(i);
  }

  // localSum is modified, do not use afterwards
  allreduceSum(localSum, globalSum, 1);

  // old loop over all slaves solution
  /*
  if(_isSlave){
    _communication->send(localSum, 0);
    _communication->receive(globalSum, 0);
  }
  if(_isMaster){
    globalSum += localSum;
    for(int rankSlave = 1; rankSlave < _size; rankSlave++){
      _communication->receive(localSum, rankSlave);
      globalSum += localSum;
    }
    for(int rankSlave = 1; rankSlave < _size; rankSlave++){
      _communication->send(globalSum, rankSlave);
    }
  }
  */
  return globalSum;
}

void MasterSlave:: reset()
{
  P_TRACE();
  _isMaster = false;
  _isSlave = false;
  _rank = -1;
  _size = -1;
}


void
MasterSlave::reduceSum(double* sendData, double* rcvData, int size) {
  P_TRACE();

  if (not _isMaster && not _isSlave) {
    return;
  }

  P_ASSERT(_communication.get() != nullptr);
  P_ASSERT(_communication->isConnected());

  if (_isSlave) {
    // send local result to master
    _communication->reduceSum(sendData, rcvData, size, 0);
  }

  if (_isMaster) {
    // receive local results from slaves, apply SUM
    _communication->reduceSum(sendData, rcvData, size);
  }
}

void
MasterSlave::reduceSum(int& sendData, int& rcvData, int size) {
  P_TRACE();

  if (not _isMaster && not _isSlave) {
    return;
  }

  P_ASSERT(_communication.get() != nullptr);
  P_ASSERT(_communication->isConnected());

  if (_isSlave) {
    // send local result to master
    _communication->reduceSum(sendData, rcvData, 0);
  }

  if (_isMaster) {
    // receive local results from slaves, apply SUM
    _communication->reduceSum(sendData, rcvData);
  }
}

void
MasterSlave::allreduceSum(double* sendData, double* rcvData, int size) {
  P_TRACE();

  if (not _isMaster && not _isSlave) {
    return;
  }

  P_ASSERT(_communication.get() != nullptr);
  P_ASSERT(_communication->isConnected());

  if (_isSlave) {
    // send local result to master, receive reduced result from master
    _communication->allreduceSum(sendData, rcvData, size, 0);
  }

  if (_isMaster) {
    // receive local results from slaves, apply SUM, send reduced result to slaves
    _communication->allreduceSum(sendData, rcvData, size);
  }
}

void
MasterSlave::allreduceSum(double& sendData, double& rcvData, int size) {
  P_TRACE();

  if (not _isMaster && not _isSlave) {
    return;
  }

  P_ASSERT(_communication.get() != nullptr);
  P_ASSERT(_communication->isConnected());

  if (_isSlave) {
    // send local result to master, receive reduced result from master
    _communication->allreduceSum(sendData, rcvData, 0);
  }

  if (_isMaster) {
    // receive local results from slaves, apply SUM, send reduced result to slaves
    _communication->allreduceSum(sendData, rcvData);
  }
}

void
MasterSlave::allreduceSum(int& sendData, int& rcvData, int size) {
  P_TRACE();

  if (not _isMaster && not _isSlave) {
    return;
  }

  P_ASSERT(_communication.get() != nullptr);
  P_ASSERT(_communication->isConnected());

  if (_isSlave) {
    // send local result to master, receive reduced result from master
    _communication->allreduceSum(sendData, rcvData, 0);
  }

  if (_isMaster) {
    // receive local results from slaves, apply SUM, send reduced result to slaves
    _communication->allreduceSum(sendData, rcvData);
  }
}

void
MasterSlave::broadcast(bool& value) {
  P_TRACE();

  if (not _isMaster && not _isSlave) {
    return;
  }

  P_ASSERT(_communication.get() != nullptr);
  P_ASSERT(_communication->isConnected());

  if (_isMaster) {
    // Broadcast (send) value.
    _communication->broadcast(value);
  }

  if (_isSlave) {
    // Broadcast (receive) value.
    _communication->broadcast(value, 0);
  }
}


void
MasterSlave::broadcast(double& value) {
  P_TRACE();

  if (not _isMaster && not _isSlave) {
    return;
  }

  P_ASSERT(_communication.get() != nullptr);
  P_ASSERT(_communication->isConnected());

  if (_isMaster) {
    // Broadcast (send) value.
    _communication->broadcast(value);
  }

  if (_isSlave) {
    // Broadcast (receive) value.
    _communication->broadcast(value, 0);
  }
}

void
MasterSlave::broadcast(double* values, int size) {
  P_TRACE();

  if (not _isMaster && not _isSlave) {
    return;
  }

  P_ASSERT(_communication.get() != nullptr);
  P_ASSERT(_communication->isConnected());

  if (_isMaster) {
    // Broadcast (send) value.
    _communication->broadcast(values, size);
  }

  if (_isSlave) {
    // Broadcast (receive) value.
    _communication->broadcast(values, size, 0);
  }
}

}} // precice, utils

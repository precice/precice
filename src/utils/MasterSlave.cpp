//#ifndef PRECICE_NO_MPI

#include "MasterSlave.hpp"

#include "EventTimings.hpp"
#include "Globals.hpp"

#include <math.h>

namespace precice {
namespace utils {

int MasterSlave::_rank = -1;
int MasterSlave::_size = -1;
// One day somebody would want to choose master dynamically (e.g. from
// configuration). As a result, `_masterMode' would no longer correspond to 0th
// process. For now we just hardcode `_masterRank' assignment to `0'.
int MasterSlave::_masterRank = 0;
bool MasterSlave::_masterMode = false;
bool MasterSlave::_slaveMode = false;
com::Communication::SharedPointer MasterSlave::_communication;


tarch::logging::Log MasterSlave:: _log ( "precice::utils::MasterSlave" );

void MasterSlave:: configure(int rank, int size)
{
  preciceTrace2("initialize()", rank, size);
  preciceCheck(size>=2, "initialize()", "You cannot use a master with a serial participant.");
  _rank = rank;
  _size = size;
  assertion(_rank != -1 && _size != -1);
  _masterMode = (rank==0);
  _slaveMode = (rank!=0);
  preciceDebug("slaveMode: " << _slaveMode <<", masterMode: " << _masterMode);
}

double MasterSlave:: l2norm(const DynVector& vec)
{
  preciceTrace("l2norm()");

  if(not _masterMode && not _slaveMode){ //old case
    return tarch::la::norm2(vec);
  }

  assertion(_communication.get() != nullptr);
  assertion(_communication->isConnected());
  double localSum2 = 0.0;
  double globalSum2 = 0.0;

  for(int i=0; i<vec.size(); i++){
    localSum2 += vec[i]*vec[i];
  }

  // localSum is modified, do not use afterwards
  allreduceSum(localSum2, globalSum2, 1);

  /* old loop over all slaves solution
  if(_slaveMode){
    _communication->send(localSum2, 0);
    _communication->receive(globalSum2, 0);
  }
  if(_masterMode){
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

double MasterSlave:: l2norm(const EigenVector& vec)
{
  preciceTrace("l2norm()");

  if(not _masterMode && not _slaveMode){ //old case
    return vec.norm();
  }

  assertion(_communication.get() != nullptr);
  assertion(_communication->isConnected());
  double localSum2 = 0.0;
  double globalSum2 = 0.0;

  for(int i=0; i<vec.size(); i++){
    localSum2 += vec(i)*vec(i);
  }

  // localSum is modified, do not use afterwards
  allreduceSum(localSum2, globalSum2, 1);
   /* old loop over all slaves solution
  if(_slaveMode){
    _communication->send(localSum2, 0);
    _communication->receive(globalSum2, 0);
  }
  if(_masterMode){
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


double MasterSlave:: dot(const DynVector& vec1, const DynVector& vec2)
{
  preciceTrace("dot()");

  if(not _masterMode && not _slaveMode){ //old case
    return tarch::la::dot(vec1, vec2);
  }

  assertion(_communication.get() != nullptr);
  assertion(_communication->isConnected());
  assertion(vec1.size()==vec2.size(), vec1.size(), vec2.size());
  double localSum = 0.0;
  double globalSum = 0.0;

  for(int i=0; i<vec1.size(); i++){
    localSum += vec1[i]*vec2[i];
  }

  // localSum is modified, do not use afterwards
  allreduceSum(localSum, globalSum, 1);

  /* old loop over all slaves solution
  if(_slaveMode){
    _communication->send(localSum, 0);
    _communication->receive(globalSum, 0);
  }
  if(_masterMode){
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

double MasterSlave:: dot(const EigenVector& vec1, const EigenVector& vec2)
{
  preciceTrace("dot()");

  if(not _masterMode && not _slaveMode){ //old case
    return vec1.dot(vec2);
  }

  assertion(_communication.get() != nullptr);
  assertion(_communication->isConnected());
  assertion(vec1.size()==vec2.size(), vec1.size(), vec2.size());
  double localSum = 0.0;
  double globalSum = 0.0;

  for(int i=0; i<vec1.size(); i++){
    localSum += vec1(i)*vec2(i);
  }

  // localSum is modified, do not use afterwards
  allreduceSum(localSum, globalSum, 1);

  // old loop over all slaves solution
  /*
  if(_slaveMode){
    _communication->send(localSum, 0);
    _communication->receive(globalSum, 0);
  }
  if(_masterMode){
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
  preciceTrace("reset()");
  _masterMode = false;
  _slaveMode = false;
  _rank = -1;
  _size = -1;
}


void
MasterSlave::reduceSum(double* sendData, double* rcvData, int size) {
  preciceTrace("reduceSum(double*)");

  if (not _masterMode && not _slaveMode) {
    return;
  }

  assertion(_communication.get() != nullptr);
  assertion(_communication->isConnected());

  //Event e("MasterSlave::allreduce_sum");

  if (_slaveMode) {
    // send local result to master
    _communication->reduceSum(sendData, rcvData, size, 0);
  }

  if (_masterMode) {
    // receive local results from slaves, apply SUM
    _communication->reduceSum(sendData, rcvData, size);
  }
}

void
MasterSlave::reduceSum(int& sendData, int& rcvData, int size) {
  preciceTrace("reduceSum(int)");

  if (not _masterMode && not _slaveMode) {
    return;
  }

  assertion(_communication.get() != nullptr);
  assertion(_communication->isConnected());

  //Event e("MasterSlave::allreduce_sum");

  if (_slaveMode) {
    // send local result to master
    _communication->reduceSum(sendData, rcvData, 0);
  }

  if (_masterMode) {
    // receive local results from slaves, apply SUM
    _communication->reduceSum(sendData, rcvData);
  }
}

void
MasterSlave::allreduceSum(double* sendData, double* rcvData, int size) {
  preciceTrace("allreduceSum(double*)");

  if (not _masterMode && not _slaveMode) {
    return;
  }

  assertion(_communication.get() != nullptr);
  assertion(_communication->isConnected());

  //Event e("MasterSlave::allreduce_sum");

  if (_slaveMode) {
    // send local result to master, receive reduced result from master
    _communication->allreduceSum(sendData, rcvData, size, 0);
  }

  if (_masterMode) {
    // receive local results from slaves, apply SUM, send reduced result to slaves
    _communication->allreduceSum(sendData, rcvData, size);
  }
}

void
MasterSlave::allreduceSum(double& sendData, double& rcvData, int size) {
  preciceTrace("allreduceSum(double)");

  if (not _masterMode && not _slaveMode) {
    return;
  }

  assertion(_communication.get() != nullptr);
  assertion(_communication->isConnected());

  //Event e("MasterSlave::allreduce_sum");

  if (_slaveMode) {
    // send local result to master, receive reduced result from master
    _communication->allreduceSum(sendData, rcvData, 0);
  }

  if (_masterMode) {
    // receive local results from slaves, apply SUM, send reduced result to slaves
    _communication->allreduceSum(sendData, rcvData);
  }
}

void
MasterSlave::allreduceSum(int& sendData, int& rcvData, int size) {
  preciceTrace("allreduceSum(double)");

  if (not _masterMode && not _slaveMode) {
    return;
  }

  assertion(_communication.get() != nullptr);
  assertion(_communication->isConnected());

  //Event e("MasterSlave::allreduce_sum");

  if (_slaveMode) {
    // send local result to master, receive reduced result from master
    _communication->allreduceSum(sendData, rcvData, 0);
  }

  if (_masterMode) {
    // receive local results from slaves, apply SUM, send reduced result to slaves
    _communication->allreduceSum(sendData, rcvData);
  }
}

void
MasterSlave::broadcast(bool& value) {
  preciceTrace("broadcast(bool&)");

  if (not _masterMode && not _slaveMode) {
    return;
  }

  assertion(_communication.get() != nullptr);
  assertion(_communication->isConnected());

  //Event e("MasterSlave::broadcast");

  if (_masterMode) {
    // Broadcast (send) value.
    _communication->broadcast(value);
  }

  if (_slaveMode) {
    // Broadcast (receive) value.
    _communication->broadcast(value, 0);
  }
}


void
MasterSlave::broadcast(double& value) {
  preciceTrace("broadcast(double&)");

  if (not _masterMode && not _slaveMode) {
    return;
  }

  assertion(_communication.get() != nullptr);
  assertion(_communication->isConnected());

  //Event e("MasterSlave::broadcast");

  if (_masterMode) {
    // Broadcast (send) value.
    _communication->broadcast(value);
  }

  if (_slaveMode) {
    // Broadcast (receive) value.
    _communication->broadcast(value, 0);
  }
}

void
MasterSlave::broadcast(double* values, int size) {
  preciceTrace("broadcast(double*)");

  if (not _masterMode && not _slaveMode) {
    return;
  }

  assertion(_communication.get() != nullptr);
  assertion(_communication->isConnected());

  //Event e("MasterSlave::broadcast");

  if (_masterMode) {
    // Broadcast (send) value.
    _communication->broadcast(values, size);
  }

  if (_slaveMode) {
    // Broadcast (receive) value.
    _communication->broadcast(values, size, 0);
  }
}

}} // precice, utils

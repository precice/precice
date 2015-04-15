// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
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

  assertion(_communication.get() != NULL);
  assertion(_communication->isConnected());
  double localSum2 = 0.0;
  double globalSum2 = 0.0;

  for(int i=0; i<vec.size(); i++){
    localSum2 += vec[i]*vec[i];
  }

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
  return sqrt(globalSum2);
}


double MasterSlave:: dot(const DynVector& vec1, const DynVector& vec2)
{
  preciceTrace("dot()");

  if(not _masterMode && not _slaveMode){ //old case
    return tarch::la::dot(vec1, vec2);
  }

  assertion(_communication.get() != NULL);
  assertion(_communication->isConnected());
  assertion2(vec1.size()==vec2.size(), vec1.size(), vec2.size());
  double localSum = 0.0;
  double globalSum = 0.0;

  for(int i=0; i<vec1.size(); i++){
    localSum += vec1[i]*vec2[i];
  }

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
MasterSlave::broadcast(bool& value) {
  preciceTrace("broadcast(bool&)");

  if (not _masterMode && not _slaveMode) {
    return;
  }

  assertion(_communication.get() != NULL);
  assertion(_communication->isConnected());

  Event e("MasterSlave::broadcast");

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

  assertion(_communication.get() != NULL);
  assertion(_communication->isConnected());

  Event e("MasterSlave::broadcast");

  if (_masterMode) {
    // Broadcast (send) value.
    _communication->broadcast(value);
  }

  if (_slaveMode) {
    // Broadcast (receive) value.
    _communication->broadcast(value, 0);
  }
}

}} // precice, utils

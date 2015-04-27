// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License

#include "M2N.hpp"

#include "DistributedCommunication.hpp"
#include "DistributedComFactory.hpp"
#include "GatherScatterCommunication.hpp"
#include "com/Communication.hpp"
#include "utils/EventTimings.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/Publisher.hpp"
#include "mesh/Mesh.hpp"

using precice::utils::Event;
using precice::utils::Publisher;

namespace precice {
namespace m2n {

tarch::logging::Log M2N::_log("precice::m2n::M2N");

M2N:: M2N(  com::Communication::SharedPointer masterCom, DistributedComFactory::SharedPointer distrFactory )
:
  _distComs(),
  _masterCom(masterCom),
  _distrFactory(distrFactory),
  _isMasterConnected(false),
  _areSlavesConnected(false)
{}

M2N:: ~M2N()
{
  if (isConnected()){
    closeConnection();
  }
}

bool M2N:: isConnected()
{
  return _isMasterConnected;
}

void M2N:: acceptMasterConnection (
  const std::string& nameAcceptor,
  const std::string& nameRequester)
{
  preciceTrace2("acceptMasterConnection()", nameAcceptor, nameRequester);

  Event e("M2N::acceptMasterConnection");

  if(not utils::MasterSlave::_slaveMode){
    assertion(_masterCom.use_count()>0);
    _masterCom->acceptConnection(nameAcceptor, nameRequester, 0, 1);
    _isMasterConnected = _masterCom->isConnected();
  }

  utils::MasterSlave::broadcast(_isMasterConnected);
}

void M2N:: requestMasterConnection (
  const std::string& nameAcceptor,
  const std::string& nameRequester)
{
  preciceTrace2("requestMasterConnection()", nameAcceptor, nameRequester);

  Event e("M2N::requestMasterConnection");

  if(not utils::MasterSlave::_slaveMode){
    assertion(_masterCom.use_count()>0);

    Publisher::ScopedSetEventNamePrefix ssenp(
        "M2N::requestMasterConnection"
        "/");

    _masterCom->requestConnection(nameAcceptor, nameRequester, 0, 1);
    _isMasterConnected = _masterCom->isConnected();
  }

  utils::MasterSlave::broadcast(_isMasterConnected);
}

void M2N:: acceptSlavesConnection (
  const std::string& nameAcceptor,
  const std::string& nameRequester)
{
  preciceTrace2("acceptSlavesConnection()", nameAcceptor, nameRequester);
  _areSlavesConnected = true;
  for( const auto& pair : _distComs){
    pair.second->acceptConnection(nameAcceptor, nameRequester);
    _areSlavesConnected = _areSlavesConnected && pair.second->isConnected();
  }
  assertion(_areSlavesConnected);
}

void M2N:: requestSlavesConnection (
  const std::string& nameAcceptor,
  const std::string& nameRequester)
{
  preciceTrace2("requestSlavesConnection()", nameAcceptor, nameRequester);
  _areSlavesConnected = true;
  for( const auto& pair : _distComs){
    pair.second->requestConnection(nameAcceptor, nameRequester);
    _areSlavesConnected = _areSlavesConnected && pair.second->isConnected();
  }
  assertion(_areSlavesConnected);
}

void M2N:: closeConnection()
{
  preciceTrace("closeConnection()");
  if(not utils::MasterSlave::_slaveMode && _masterCom->isConnected()){
    _masterCom->closeConnection();
    _isMasterConnected = false;
  }

  utils::MasterSlave::broadcast(_isMasterConnected);

  if(utils::MasterSlave::_slaveMode || utils::MasterSlave::_masterMode){
    _areSlavesConnected = false;
    for( const auto& pair : _distComs){
      pair.second->closeConnection();
      _areSlavesConnected = _areSlavesConnected || pair.second->isConnected();
    }
    assertion(not _areSlavesConnected);
  }
}

com::Communication::SharedPointer M2N:: getMasterCommunication()
{
  assertion(not utils::MasterSlave::_slaveMode);
  return _masterCom; //TODO maybe it would be a nicer design to not offer this
}

void M2N:: createDistributedCommunication(mesh::PtrMesh mesh){
  DistributedCommunication::SharedPointer distCom =
      _distrFactory->newDistributedCommunication(mesh);
  _distComs[mesh->getID()] = distCom;
}

void M2N:: startSendPackage ( int rankReceiver )
{
  if(not utils::MasterSlave::_slaveMode){
    _masterCom->startSendPackage(rankReceiver);
  }
}

void M2N:: finishSendPackage()
{
  if(not utils::MasterSlave::_slaveMode){
    _masterCom->finishSendPackage();
  }
}

int M2N:: startReceivePackage ( int rankSender )
{
  if(not utils::MasterSlave::_slaveMode){
    return _masterCom->startReceivePackage(rankSender);
  }
  return -1;
}

void M2N:: finishReceivePackage()
{
  if(not utils::MasterSlave::_slaveMode){
    _masterCom->finishReceivePackage();
  }
}

void M2N:: send (
  double* itemsToSend,
  int     size,
  int     meshID,
  int     valueDimension )
{
  if(utils::MasterSlave::_slaveMode || utils::MasterSlave::_masterMode){
    assertion(_areSlavesConnected);
    assertion(_distComs.find(meshID) != _distComs.end());
    assertion(_distComs[meshID].get() != NULL);
    _distComs[meshID]->send(itemsToSend,size,valueDimension);
  }
  else{//coupling mode
    assertion(_isMasterConnected);
    _masterCom->send(itemsToSend, size, 0);
  }
}

void M2N:: send (
  bool   itemToSend)
{
  preciceTrace1("send(bool)", utils::MasterSlave::_rank);
  if(not utils::MasterSlave::_slaveMode){
    _masterCom->send(itemToSend, 0);
  }
}

void M2N:: send (
  double itemToSend)
{
  preciceTrace1("send(double)", utils::MasterSlave::_rank);
  if(not utils::MasterSlave::_slaveMode){
    _masterCom->send(itemToSend, 0);
  }
}

void M2N:: receive (
  double* itemsToReceive,
  int     size,
  int     meshID,
  int     valueDimension )
{
  if(utils::MasterSlave::_slaveMode || utils::MasterSlave::_masterMode){
    assertion(_areSlavesConnected);
    assertion(_distComs.find(meshID) != _distComs.end());
    assertion(_distComs[meshID].get() != NULL);
    _distComs[meshID]->receive(itemsToReceive,size,valueDimension);
  }
  else{//coupling mode
    assertion(_isMasterConnected);
    _masterCom->receive(itemsToReceive, size, 0);
  }
}

void M2N:: receive (
  bool&  itemToReceive )
{
  preciceTrace1("receive(bool)", utils::MasterSlave::_rank);
  if(not utils::MasterSlave::_slaveMode){
    _masterCom->receive(itemToReceive, 0);
  }

  utils::MasterSlave::broadcast(itemToReceive);

  preciceDebug("receive(bool): " << itemToReceive);
}

void M2N:: receive (
  double&  itemToReceive)
{
  preciceTrace1("receive(double)", utils::MasterSlave::_rank);
  if(not utils::MasterSlave::_slaveMode){ //coupling mode
    _masterCom->receive(itemToReceive, 0);
  }

  utils::MasterSlave::broadcast(itemToReceive);

  preciceDebug("receive(double): " << itemToReceive);
}


}} // namespace precice, m2n

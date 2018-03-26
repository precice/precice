
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

extern bool testMode;

namespace m2n {

logging::Logger M2N::_log("m2n::M2N");

M2N:: M2N(  com::PtrCommunication masterCom, DistributedComFactory::SharedPointer distrFactory )
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
  TRACE(nameAcceptor, nameRequester);

  //Event e("M2N::acceptMasterConnection");

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
  TRACE(nameAcceptor, nameRequester);

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
  TRACE(nameAcceptor, nameRequester);
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
  TRACE(nameAcceptor, nameRequester);
  _areSlavesConnected = true;
  for( const auto& pair : _distComs){
    pair.second->requestConnection(nameAcceptor, nameRequester);
    _areSlavesConnected = _areSlavesConnected && pair.second->isConnected();
  }
  assertion(_areSlavesConnected);
}

void M2N:: closeConnection()
{
  TRACE();
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

com::PtrCommunication M2N:: getMasterCommunication()
{
  assertion(not utils::MasterSlave::_slaveMode);
  return _masterCom; //TODO maybe it would be a nicer design to not offer this
}

void M2N:: createDistributedCommunication(mesh::PtrMesh mesh){
  DistributedCommunication::SharedPointer distCom =
      _distrFactory->newDistributedCommunication(mesh);
  _distComs[mesh->getID()] = distCom;
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
    assertion(_distComs[meshID].get() != nullptr);

#ifdef M2N_PRE_SYNCHRONIZE
    if(not precice::testMode){
//      Event e("M2N::send/synchronize", true);

      if(not utils::MasterSlave::_slaveMode){
        bool ack;

        _masterCom->send(ack, 0);
        _masterCom->receive(ack, 0);
        _masterCom->send(ack, 0);
      }
    }
#endif

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
  TRACE(utils::MasterSlave::_rank);
  if(not utils::MasterSlave::_slaveMode){
    _masterCom->send(itemToSend, 0);
  }
}


void M2N:: send (
  double itemToSend)
{
  TRACE(utils::MasterSlave::_rank);
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
    assertion(_distComs[meshID].get() != nullptr);

#ifdef M2N_PRE_SYNCHRONIZE
    if(not precice::testMode){
//      Event e("M2N::receive/synchronize", true);

      if(not utils::MasterSlave::_slaveMode){
        bool ack;

        _masterCom->receive(ack, 0);
        _masterCom->send(ack, 0);
        _masterCom->receive(ack, 0);
      }
    }
#endif

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
  TRACE(utils::MasterSlave::_rank);
  if(not utils::MasterSlave::_slaveMode){
    _masterCom->receive(itemToReceive, 0);
  }

  utils::MasterSlave::broadcast(itemToReceive);

  DEBUG("receive(bool): " << itemToReceive);
}


void M2N:: receive (
  double&  itemToReceive)
{
  TRACE(utils::MasterSlave::_rank);
  if(not utils::MasterSlave::_slaveMode){ //coupling mode
    _masterCom->receive(itemToReceive, 0);
  }

  utils::MasterSlave::broadcast(itemToReceive);

  DEBUG("receive(double): " << itemToReceive);
}


}} // namespace precice, m2n

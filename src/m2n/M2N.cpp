#include "M2N.hpp"
#include "DistributedComFactory.hpp"
#include "DistributedCommunication.hpp"
#include "com/Communication.hpp"
#include "mesh/Mesh.hpp"
#include "utils/EventTimings.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/Publisher.hpp"

using precice::utils::Event;
using precice::utils::Publisher;

namespace precice
{
extern bool testMode;
extern bool syncMode;

namespace m2n
{

M2N::M2N(com::PtrCommunication masterCom, DistributedComFactory::SharedPointer distrFactory)
    : _masterCom(masterCom),
      _distrFactory(distrFactory)
{
}

M2N::~M2N()
{
  if (isConnected()) {
    closeConnection();
  }
}

bool M2N::isConnected()
{
  return _isMasterConnected;
}

void M2N::acceptMasterConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  TRACE(acceptorName, requesterName);

  Event e("m2n.acceptMasterConnection");

  if (not utils::MasterSlave::_slaveMode) {
    assertion(_masterCom.use_count() > 0);
    _masterCom->acceptConnection(acceptorName, requesterName, utils::MasterSlave::_rank);
    _isMasterConnected = _masterCom->isConnected();
  }

  utils::MasterSlave::broadcast(_isMasterConnected);
}

void M2N::requestMasterConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  TRACE(acceptorName, requesterName);

  Event e("m2n.requestMasterConnection");

  if (not utils::MasterSlave::_slaveMode) {
    assertion(_masterCom.use_count() > 0);

    utils::ScopedEventPrefix sep("M2N::requestMasterConnection/");

    _masterCom->requestConnection(acceptorName, requesterName, 0, 1);
    _isMasterConnected = _masterCom->isConnected();
  }

  utils::MasterSlave::broadcast(_isMasterConnected);
}

void M2N::acceptSlavesConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  TRACE(acceptorName, requesterName);
  Event e("m2n.acceptSlavesConnection");

  _areSlavesConnected = true;
  for (const auto &pair : _distComs) {
    pair.second->acceptConnection(acceptorName, requesterName);
    _areSlavesConnected = _areSlavesConnected && pair.second->isConnected();
  }
  assertion(_areSlavesConnected);
}

void M2N::acceptSlavesPreConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  TRACE(acceptorName, requesterName);
  _areSlavesConnected = true;
  for (const auto &pair : _distComs) {
    pair.second->acceptPreConnection(acceptorName, requesterName);
    _areSlavesConnected = _areSlavesConnected && pair.second->isConnected();
    }
  assertion(_areSlavesConnected);
}

void M2N::requestSlavesConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  TRACE(acceptorName, requesterName);
  Event e("m2n.requestSlavesConnection");

  _areSlavesConnected = true;
  for (const auto &pair : _distComs) {
    pair.second->requestConnection(acceptorName, requesterName);
    _areSlavesConnected = _areSlavesConnected && pair.second->isConnected();
  }
  assertion(_areSlavesConnected);
}

void M2N::requestSlavesPreConnection(const std::string &acceptorName,
                                     const std::string &requesterName)
{
  TRACE(acceptorName, requesterName);
  _areSlavesConnected = true;
  for (const auto &pair : _distComs) {
    pair.second->requestPreConnection(acceptorName, requesterName);
    _areSlavesConnected = _areSlavesConnected && pair.second->isConnected();
    }
  assertion(_areSlavesConnected);
}

void M2N::updateAcceptorCommunicationMap()
{
  for (const auto &pair : _distComs) {
    pair.second->updateAcceptorCommunicationMap();
  }
}

void M2N::updateRequesterCommunicationMap()
{ 
  for (const auto &pair : _distComs) {
    pair.second->updateRequesterCommunicationMap();
  }
}


void M2N::closeConnection()
{
  TRACE();
  if (not utils::MasterSlave::_slaveMode && _masterCom->isConnected()) {
    _masterCom->closeConnection();
    _isMasterConnected = false;
  }

  utils::MasterSlave::broadcast(_isMasterConnected);

  if (utils::MasterSlave::_slaveMode || utils::MasterSlave::_masterMode) {
    _areSlavesConnected = false;
    for (const auto &pair : _distComs) {
      pair.second->closeConnection();
      _areSlavesConnected = _areSlavesConnected || pair.second->isConnected();
    }
    assertion(not _areSlavesConnected);
  }
}

com::PtrCommunication M2N::getMasterCommunication()
{
  assertion(not utils::MasterSlave::_slaveMode);
  return _masterCom; /// @todo maybe it would be a nicer design to not offer this
}

void M2N::createDistributedCommunication(mesh::PtrMesh mesh)
{
  DistributedCommunication::SharedPointer distCom = _distrFactory->newDistributedCommunication(mesh);
  _distComs[mesh->getID()]                        = distCom;
}

void M2N::send(
    double *itemsToSend,
    int     size,
    int     meshID,
    int     valueDimension)
{
  if (utils::MasterSlave::_slaveMode || utils::MasterSlave::_masterMode) {
    assertion(_areSlavesConnected);
    assertion(_distComs.find(meshID) != _distComs.end());
    assertion(_distComs[meshID].get() != nullptr);

    if (precice::syncMode and not precice::testMode) {
      if (not utils::MasterSlave::_slaveMode) {
        bool ack = true;
        _masterCom->send(ack, 0);
        _masterCom->receive(ack, 0);
        _masterCom->send(ack, 0);
      }
    }

    _distComs[meshID]->send(itemsToSend, size, valueDimension);
  } else { //coupling mode
    assertion(_isMasterConnected);
    _masterCom->send(itemsToSend, size, 0);
  }
}

void M2N::send(bool itemToSend)
{
  TRACE(utils::MasterSlave::_rank);
  if (not utils::MasterSlave::_slaveMode) {
    _masterCom->send(itemToSend, 0);
  }
}

void M2N::send(double itemToSend)
{
  TRACE(utils::MasterSlave::_rank);
  if (not utils::MasterSlave::_slaveMode) {
    _masterCom->send(itemToSend, 0);
  }
}

void M2N::broadcastSendLocalMesh(mesh::Mesh &mesh)
{
  int meshID = mesh.getID();
  if (utils::MasterSlave::_slaveMode || utils::MasterSlave::_masterMode) {
    assertion(_areSlavesConnected);
    assertion(_distComs.find(meshID) != _distComs.end());
    assertion(_distComs[meshID].get() != nullptr);    
    _distComs[meshID]->sendMesh(mesh);
  } else { //coupling mode
  }
}

void M2N::sendCommunicationMap(mesh::Mesh::FeedbackMap &localCommunicationMap,
                               mesh::Mesh &mesh)
{ 
  if (utils::MasterSlave::_slaveMode || utils::MasterSlave::_masterMode) {
    int meshID = mesh.getID();
    assertion(_areSlavesConnected);
    _distComs[meshID]->sendCommunicationMap(localCommunicationMap);
  } else { //coupling mode
  }
}

void M2N::receive(double *itemsToReceive,
                  int     size,
                  int     meshID,
                  int     valueDimension)
{
  if (utils::MasterSlave::_slaveMode || utils::MasterSlave::_masterMode) {
    assertion(_areSlavesConnected);
    assertion(_distComs.find(meshID) != _distComs.end());
    assertion(_distComs[meshID].get() != nullptr);

    if (precice::syncMode and not precice::testMode) {
      if (not utils::MasterSlave::_slaveMode) {
        bool ack;

        _masterCom->receive(ack, 0);
        _masterCom->send(ack, 0);
        _masterCom->receive(ack, 0);
      }
    }

    _distComs[meshID]->receive(itemsToReceive, size, valueDimension);
  } else { //coupling mode
    assertion(_isMasterConnected);
    _masterCom->receive(itemsToReceive, size, 0);
  }
}

void M2N::receive(bool &itemToReceive)
{
  TRACE(utils::MasterSlave::_rank);
  if (not utils::MasterSlave::_slaveMode) {
    _masterCom->receive(itemToReceive, 0);
  }

  utils::MasterSlave::broadcast(itemToReceive);

  DEBUG("receive(bool): " << itemToReceive);
}

void M2N::receive(double &itemToReceive)
{
  TRACE(utils::MasterSlave::_rank);
  if (not utils::MasterSlave::_slaveMode) { //coupling mode
    _masterCom->receive(itemToReceive, 0);
  }

  utils::MasterSlave::broadcast(itemToReceive);

  DEBUG("receive(double): " << itemToReceive);
}

void M2N::broadcastReceiveLocalMesh(mesh::Mesh &mesh)
{
  int meshID = mesh.getID();
  if (utils::MasterSlave::_slaveMode || utils::MasterSlave::_masterMode) {
    assertion(_areSlavesConnected);
    assertion(_distComs.find(meshID) != _distComs.end());
    assertion(_distComs[meshID].get() != nullptr);
    _distComs[meshID]->receiveMesh(mesh);
  } else { //coupling mode
  }
}

void M2N::receiveCommunicationMap(mesh::Mesh::FeedbackMap &localCommunicationMap,
                                  mesh::Mesh &mesh)
{
  if (utils::MasterSlave::_slaveMode || utils::MasterSlave::_masterMode) {
    int meshID = mesh.getID();
    assertion(_areSlavesConnected);
    _distComs[meshID]->receiveCommunicationMap(localCommunicationMap);
  } else { //coupling mode
  }
}

} // namespace m2n
} // namespace precice

#include "M2N.hpp"
#include "DistributedComFactory.hpp"
#include "DistributedCommunication.hpp"
#include "com/Communication.hpp"
#include "mesh/Mesh.hpp"
#include "utils/Event.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/Publisher.hpp"

using precice::utils::Event;
using precice::utils::Publisher;

namespace precice
{
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

  Event e("m2n.acceptMasterConnection", precice::syncMode);

  if (not utils::MasterSlave::isSlave()) {
    assertion(_masterCom);
    _masterCom->acceptConnection(acceptorName, requesterName, utils::MasterSlave::getRank());
    _isMasterConnected = _masterCom->isConnected();
  }

  utils::MasterSlave::broadcast(_isMasterConnected);
}

void M2N::requestMasterConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  TRACE(acceptorName, requesterName);

  Event e("m2n.requestMasterConnection", precice::syncMode);

  if (not utils::MasterSlave::isSlave()) {
    assertion(_masterCom);

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
  Event e("m2n.acceptSlavesConnection", precice::syncMode);

  _areSlavesConnected = true;
  for (const auto &pair : _distComs) {
    pair.second->acceptConnection(acceptorName, requesterName);
    _areSlavesConnected = _areSlavesConnected && pair.second->isConnected();
  }
  assertion(_areSlavesConnected);
}

void M2N::requestSlavesConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  TRACE(acceptorName, requesterName);
  Event e("m2n.requestSlavesConnection", precice::syncMode);

  _areSlavesConnected = true;
  for (const auto &pair : _distComs) {
    pair.second->requestConnection(acceptorName, requesterName);
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

void M2N::requestSlavesPreConnection(
  const std::string &acceptorName,
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

void M2N::closeConnection()
{
  TRACE();
  if (not utils::MasterSlave::isSlave() && _masterCom->isConnected()) {
    _masterCom->closeConnection();
    _isMasterConnected = false;
  }

  utils::MasterSlave::broadcast(_isMasterConnected);

  if (utils::MasterSlave::isSlave() || utils::MasterSlave::isMaster()) {
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
  assertion(not utils::MasterSlave::isSlave());
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
  if (utils::MasterSlave::isSlave() || utils::MasterSlave::isMaster()) {
    assertion(_areSlavesConnected);
    assertion(_distComs.find(meshID) != _distComs.end());
    assertion(_distComs[meshID].get() != nullptr);

    if (precice::syncMode) {
      if (not utils::MasterSlave::isSlave()) {
        bool ack = true;
        _masterCom->send(ack, 0);
        _masterCom->receive(ack, 0);
        _masterCom->send(ack, 0);
      }
    }
    Event e("m2n.sendData", precice::syncMode);
    _distComs[meshID]->send(itemsToSend, size, valueDimension);
  } else { //coupling mode
    assertion(_isMasterConnected);
    _masterCom->send(itemsToSend, size, 0);
  }
}

void M2N::send(bool itemToSend)
{
  TRACE(utils::MasterSlave::getRank());
  if (not utils::MasterSlave::isSlave()) {
    _masterCom->send(itemToSend, 0);
  }
}

void M2N::send(double itemToSend)
{
  TRACE(utils::MasterSlave::getRank());
  if (not utils::MasterSlave::isSlave()) {
    _masterCom->send(itemToSend, 0);
  }
}

void M2N::broadcastSendLocalMesh(mesh::Mesh &mesh)
{
  int meshID = mesh.getID();
  if (utils::MasterSlave::isSlave() || utils::MasterSlave::isMaster()) {
    assertion(_areSlavesConnected);
    assertion(_distComs.find(meshID) != _distComs.end());
    assertion(_distComs[meshID].get() != nullptr);    
    _distComs[meshID]->broadcastSendMesh(mesh);
  } else { //coupling mode
  }
}

void M2N::broadcastSendLCM(std::map<int, std::vector<int>> &localCommunicationMap,
                        mesh::Mesh &mesh)
{
  if (utils::MasterSlave::isSlave() || utils::MasterSlave::isMaster()) {
    int meshID = mesh.getID();
    assertion(_areSlavesConnected);
    _distComs[meshID]->broadcastSendLCM(localCommunicationMap);
  } else { //coupling mode
  }  
}

void M2N::receive(double *itemsToReceive,
                  int     size,
                  int     meshID,
                  int     valueDimension)
{
  if (utils::MasterSlave::isSlave() || utils::MasterSlave::isMaster()) {
    assertion(_areSlavesConnected);
    assertion(_distComs.find(meshID) != _distComs.end());
    assertion(_distComs[meshID].get() != nullptr);

    if (precice::syncMode) {
      if (not utils::MasterSlave::isSlave()) {
        bool ack;

        _masterCom->receive(ack, 0);
        _masterCom->send(ack, 0);
        _masterCom->receive(ack, 0);
      }
    }
    Event e("m2n.receiveData", precice::syncMode);
    _distComs[meshID]->receive(itemsToReceive, size, valueDimension);
  } else { //coupling mode
    assertion(_isMasterConnected);
    _masterCom->receive(itemsToReceive, size, 0);
  }
}

void M2N::receive(bool &itemToReceive)
{
  TRACE(utils::MasterSlave::getRank());
  if (not utils::MasterSlave::isSlave()) {
    _masterCom->receive(itemToReceive, 0);
  }

  utils::MasterSlave::broadcast(itemToReceive);

  DEBUG("receive(bool): " << itemToReceive);
}

void M2N::receive(double &itemToReceive)
{
  TRACE(utils::MasterSlave::getRank());
  if (not utils::MasterSlave::isSlave()) { //coupling mode
    _masterCom->receive(itemToReceive, 0);
  }

  utils::MasterSlave::broadcast(itemToReceive);

  DEBUG("receive(double): " << itemToReceive);
}

void M2N::broadcastReceiveLocalMesh(mesh::Mesh &mesh)
{
  int meshID = mesh.getID();
  if (utils::MasterSlave::isSlave() || utils::MasterSlave::isMaster()) {
    assertion(_areSlavesConnected);
    assertion(_distComs.find(meshID) != _distComs.end());
    assertion(_distComs[meshID].get() != nullptr);
    _distComs[meshID]->broadcastReceiveMesh(mesh);
  } else { //coupling mode
  }
}

void M2N::broadcastReceiveLCM(std::map<int, std::vector<int>> &localCommunicationMap, mesh::Mesh &mesh)
{
  if (utils::MasterSlave::isSlave() || utils::MasterSlave::isMaster()) {
    int meshID = mesh.getID();
    assertion(_areSlavesConnected);
    _distComs[meshID]->broadcastReceiveLCM(localCommunicationMap);
  } else { //coupling mode
  }  
}

} // namespace m2n
} // namespace precice

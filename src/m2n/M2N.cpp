#include <utility>

#include "DistributedComFactory.hpp"
#include "DistributedCommunication.hpp"
#include "M2N.hpp"
#include "com/Communication.hpp"
#include "logging/LogMacros.hpp"
#include "mesh/Mesh.hpp"
#include "precice/types.hpp"
#include "utils/Event.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/assertion.hpp"

using precice::utils::Event;

namespace precice {
extern bool syncMode;

namespace m2n {

M2N::M2N(com::PtrCommunication primaryCom, DistributedComFactory::SharedPointer distrFactory, bool useOnlyMasterCom, bool useTwoLevelInit)
    : _primaryCom(std::move(primaryCom)),
      _distrFactory(std::move(distrFactory)),
      _useOnlyMasterCom(useOnlyMasterCom),
      _useTwoLevelInit(useTwoLevelInit)
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
  PRECICE_TRACE(acceptorName, requesterName);

  Event e("m2n.acceptMasterConnection", precice::syncMode);

  if (not utils::MasterSlave::isSlave()) {
    PRECICE_DEBUG("Accept primary-primary connection");
    PRECICE_ASSERT(_primaryCom);
    _primaryCom->acceptConnection(acceptorName, requesterName, "MASTERCOM", utils::MasterSlave::getRank());
    _isMasterConnected = _primaryCom->isConnected();
  }

  utils::MasterSlave::broadcast(_isMasterConnected);
}

void M2N::requestMasterConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);

  Event e("m2n.requestMasterConnection", precice::syncMode);

  if (not utils::MasterSlave::isSlave()) {
    PRECICE_ASSERT(_primaryCom);
    PRECICE_DEBUG("Request primary-primary connection");
    _primaryCom->requestConnection(acceptorName, requesterName, "MASTERCOM", 0, 1);
    _isMasterConnected = _primaryCom->isConnected();
  }

  utils::MasterSlave::broadcast(_isMasterConnected);
}

void M2N::acceptSlavesConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not _useOnlyMasterCom);
  Event e("m2n.acceptSlavesConnection", precice::syncMode);

  _areSlavesConnected = true;
  for (const auto &pair : _distComs) {
    PRECICE_DEBUG("Accept secondaries-secondaries connections");
    pair.second->acceptConnection(acceptorName, requesterName);
    _areSlavesConnected = _areSlavesConnected && pair.second->isConnected();
  }
  PRECICE_ASSERT(_areSlavesConnected);
}

void M2N::requestSlavesConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not _useOnlyMasterCom);
  Event e("m2n.requestSlavesConnection", precice::syncMode);

  _areSlavesConnected = true;
  for (const auto &pair : _distComs) {
    PRECICE_DEBUG("Request secondaries connections");
    pair.second->requestConnection(acceptorName, requesterName);
    _areSlavesConnected = _areSlavesConnected && pair.second->isConnected();
  }
  PRECICE_ASSERT(_areSlavesConnected);
}

void M2N::prepareEstablishment(const std::string &acceptorName,
                               const std::string &requesterName)
{
  PRECICE_TRACE();
  _primaryCom->prepareEstablishment(acceptorName, requesterName);
}

void M2N::cleanupEstablishment(const std::string &acceptorName,
                               const std::string &requesterName)
{
  PRECICE_TRACE();
  _primaryCom->cleanupEstablishment(acceptorName, requesterName);
}

void M2N::acceptSlavesPreConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not _useOnlyMasterCom);
  _areSlavesConnected = true;
  for (const auto &pair : _distComs) {
    pair.second->acceptPreConnection(acceptorName, requesterName);
    _areSlavesConnected = _areSlavesConnected && pair.second->isConnected();
  }
  PRECICE_ASSERT(_areSlavesConnected);
}

void M2N::requestSlavesPreConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not _useOnlyMasterCom);
  _areSlavesConnected = true;
  for (const auto &pair : _distComs) {
    pair.second->requestPreConnection(acceptorName, requesterName);
    _areSlavesConnected = _areSlavesConnected && pair.second->isConnected();
  }
  PRECICE_ASSERT(_areSlavesConnected);
}

void M2N::completeSlavesConnection()
{
  PRECICE_ASSERT(not _useOnlyMasterCom);
  for (const auto &pair : _distComs) {
    pair.second->completeSlavesConnection();
  }
}

void M2N::closeConnection()
{
  PRECICE_TRACE();
  closeMasterConnection();
  closeDistributedConnections();
}

void M2N::closeMasterConnection()
{
  PRECICE_TRACE();
  if (not utils::MasterSlave::isSlave() && _primaryCom->isConnected()) {
    _primaryCom->closeConnection();
    _isMasterConnected = false;
  }

  utils::MasterSlave::broadcast(_isMasterConnected);
  PRECICE_ASSERT(not _isMasterConnected);
}

void M2N::closeDistributedConnections()
{
  PRECICE_TRACE();
  if (_useOnlyMasterCom) {
    return;
  }

  _areSlavesConnected = false;
  for (const auto &pair : _distComs) {
    pair.second->closeConnection();
    _areSlavesConnected |= pair.second->isConnected();
  }
  PRECICE_ASSERT(not _areSlavesConnected);
}

com::PtrCommunication M2N::getMasterCommunication()
{
  PRECICE_ASSERT(not utils::MasterSlave::isSlave());
  return _primaryCom; /// @todo maybe it would be a nicer design to not offer this
}

void M2N::createDistributedCommunication(const mesh::PtrMesh &mesh)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(not _useOnlyMasterCom);
  DistributedCommunication::SharedPointer distCom = _distrFactory->newDistributedCommunication(mesh);
  _distComs[mesh->getID()]                        = distCom;
}

void M2N::send(
    precice::span<double const> itemsToSend,
    int                         meshID,
    int                         valueDimension)
{
  if (not _useOnlyMasterCom) {
    PRECICE_ASSERT(_areSlavesConnected);
    PRECICE_ASSERT(_distComs.find(meshID) != _distComs.end());
    PRECICE_ASSERT(_distComs[meshID].get() != nullptr);

    if (precice::syncMode && not utils::MasterSlave::isSlave()) {
      bool ack = true;
      _primaryCom->send(ack, 0);
      _primaryCom->receive(ack, 0);
      _primaryCom->send(ack, 0);
    }

    Event e("m2n.sendData", precice::syncMode);

    _distComs[meshID]->send(itemsToSend, valueDimension);
  } else {
    PRECICE_ASSERT(_isMasterConnected);
    _primaryCom->send(itemsToSend, 0);
  }
}

void M2N::send(bool itemToSend)
{
  PRECICE_TRACE(utils::MasterSlave::getRank());
  if (not utils::MasterSlave::isSlave()) {
    _primaryCom->send(itemToSend, 0);
  }
}

void M2N::send(double itemToSend)
{
  PRECICE_TRACE(utils::MasterSlave::getRank());
  if (not utils::MasterSlave::isSlave()) {
    _primaryCom->send(itemToSend, 0);
  }
}

void M2N::broadcastSendMesh(mesh::Mesh &mesh)
{
  MeshID meshID = mesh.getID();
  PRECICE_ASSERT(utils::MasterSlave::isParallel(),
                 "This method can only be used for parallel participants");
  PRECICE_ASSERT(_areSlavesConnected);
  PRECICE_ASSERT(_distComs.find(meshID) != _distComs.end());
  PRECICE_ASSERT(_distComs[meshID].get() != nullptr);
  _distComs[meshID]->broadcastSendMesh();
}

void M2N::scatterAllCommunicationMap(std::map<int, std::vector<int>> &localCommunicationMap,
                                     mesh::Mesh &                     mesh)
{
  PRECICE_ASSERT(utils::MasterSlave::isParallel(),
                 "This method can only be used for parallel participants");
  MeshID meshID = mesh.getID();
  PRECICE_ASSERT(_areSlavesConnected);
  _distComs[meshID]->scatterAllCommunicationMap(localCommunicationMap);
}

void M2N::broadcastSend(int &itemToSend, mesh::Mesh &mesh)
{
  PRECICE_ASSERT(utils::MasterSlave::isParallel(),
                 "This method can only be used for parallel participants");
  MeshID meshID = mesh.getID();
  PRECICE_ASSERT(_areSlavesConnected);
  _distComs[meshID]->broadcastSend(itemToSend);
}

void M2N::receive(precice::span<double> itemsToReceive,
                  int                   meshID,
                  int                   valueDimension)
{
  if (not _useOnlyMasterCom) {
    PRECICE_ASSERT(_areSlavesConnected);
    PRECICE_ASSERT(_distComs.find(meshID) != _distComs.end());
    PRECICE_ASSERT(_distComs[meshID].get() != nullptr);

    if (precice::syncMode) {
      if (not utils::MasterSlave::isSlave()) {
        bool ack;

        _primaryCom->receive(ack, 0);
        _primaryCom->send(ack, 0);
        _primaryCom->receive(ack, 0);
      }
    }

    Event e("m2n.receiveData", precice::syncMode);

    _distComs[meshID]->receive(itemsToReceive, valueDimension);
  } else {
    PRECICE_ASSERT(_isMasterConnected);
    _primaryCom->receive(itemsToReceive, 0);
  }
}

void M2N::receive(bool &itemToReceive)
{
  PRECICE_TRACE(utils::MasterSlave::getRank());
  if (not utils::MasterSlave::isSlave()) {
    _primaryCom->receive(itemToReceive, 0);
  }

  utils::MasterSlave::broadcast(itemToReceive);

  PRECICE_DEBUG("receive(bool): {}", itemToReceive);
}

void M2N::receive(double &itemToReceive)
{
  PRECICE_TRACE(utils::MasterSlave::getRank());
  if (not utils::MasterSlave::isSlave()) { //coupling mode
    _primaryCom->receive(itemToReceive, 0);
  }

  utils::MasterSlave::broadcast(itemToReceive);

  PRECICE_DEBUG("receive(double): {}", itemToReceive);
}

void M2N::broadcastReceiveAll(std::vector<int> &itemToReceive, mesh::Mesh &mesh)
{
  PRECICE_ASSERT(utils::MasterSlave::isParallel(),
                 "This method can only be used for parallel participants");
  MeshID meshID = mesh.getID();
  PRECICE_ASSERT(_areSlavesConnected);
  _distComs[meshID]->broadcastReceiveAll(itemToReceive);
}

void M2N::broadcastReceiveAllMesh(mesh::Mesh &mesh)
{
  PRECICE_ASSERT(utils::MasterSlave::isParallel(),
                 "This method can only be used for parallel participants");
  MeshID meshID = mesh.getID();
  PRECICE_ASSERT(_areSlavesConnected);
  PRECICE_ASSERT(_distComs.find(meshID) != _distComs.end());
  PRECICE_ASSERT(_distComs[meshID].get() != nullptr);
  _distComs[meshID]->broadcastReceiveAllMesh();
}

void M2N::gatherAllCommunicationMap(std::map<int, std::vector<int>> &localCommunicationMap, mesh::Mesh &mesh)
{
  PRECICE_ASSERT(utils::MasterSlave::isParallel(),
                 "This method can only be used for parallel participants");
  MeshID meshID = mesh.getID();
  PRECICE_ASSERT(_areSlavesConnected);
  _distComs[meshID]->gatherAllCommunicationMap(localCommunicationMap);
}

} // namespace m2n
} // namespace precice

#include <utility>

#include "DistributedComFactory.hpp"
#include "DistributedCommunication.hpp"
#include "M2N.hpp"
#include "com/Communication.hpp"
#include "logging/LogMacros.hpp"
#include "mesh/Mesh.hpp"
#include "precice/types.hpp"
#include "utils/Event.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

using precice::utils::Event;

namespace precice {
extern bool syncMode;

namespace m2n {

M2N::M2N(com::PtrCommunication masterCom, DistributedComFactory::SharedPointer distrFactory, bool useOnlyPrimaryCom, bool useTwoLevelInit)
    : _masterCom(std::move(masterCom)),
      _distrFactory(std::move(distrFactory)),
      _useOnlyPrimaryCom(useOnlyPrimaryCom),
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
  return _isPrimaryConnected;
}

void M2N::acceptPrimaryConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);

  Event e("m2n.acceptPrimaryConnection", precice::syncMode);

  if (not utils::IntraComm::isSecondary()) {
    PRECICE_DEBUG("Accept master-master connection");
    PRECICE_ASSERT(_masterCom);
    _masterCom->acceptConnection(acceptorName, requesterName, "MASTERCOM", utils::IntraComm::getRank());
    _isPrimaryConnected = _masterCom->isConnected();
  }

  utils::IntraComm::broadcast(_isPrimaryConnected);
}

void M2N::requestPrimaryConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);

  Event e("m2n.requestPrimaryConnection", precice::syncMode);

  if (not utils::IntraComm::isSecondary()) {
    PRECICE_ASSERT(_masterCom);
    PRECICE_DEBUG("Request master-master connection");
    _masterCom->requestConnection(acceptorName, requesterName, "MASTERCOM", 0, 1);
    _isPrimaryConnected = _masterCom->isConnected();
  }

  utils::IntraComm::broadcast(_isPrimaryConnected);
}

void M2N::acceptSecondariesConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not _useOnlyPrimaryCom);
  Event e("m2n.acceptSecondariesConnection", precice::syncMode);

  _areSecondariesConnected = true;
  for (const auto &pair : _distComs) {
    PRECICE_DEBUG("Accept slaves-slaves connections");
    pair.second->acceptConnection(acceptorName, requesterName);
    _areSecondariesConnected = _areSecondariesConnected && pair.second->isConnected();
  }
  PRECICE_ASSERT(_areSecondariesConnected);
}

void M2N::requestSecondariesConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not _useOnlyPrimaryCom);
  Event e("m2n.requestSecondariesConnection", precice::syncMode);

  _areSecondariesConnected = true;
  for (const auto &pair : _distComs) {
    PRECICE_DEBUG("Request slaves connections");
    pair.second->requestConnection(acceptorName, requesterName);
    _areSecondariesConnected = _areSecondariesConnected && pair.second->isConnected();
  }
  PRECICE_ASSERT(_areSecondariesConnected);
}

void M2N::prepareEstablishment(const std::string &acceptorName,
                               const std::string &requesterName)
{
  PRECICE_TRACE();
  _masterCom->prepareEstablishment(acceptorName, requesterName);
}

void M2N::cleanupEstablishment(const std::string &acceptorName,
                               const std::string &requesterName)
{
  PRECICE_TRACE();
  _masterCom->cleanupEstablishment(acceptorName, requesterName);
}

void M2N::acceptSecondariesPreConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not _useOnlyPrimaryCom);
  _areSecondariesConnected = true;
  for (const auto &pair : _distComs) {
    pair.second->acceptPreConnection(acceptorName, requesterName);
    _areSecondariesConnected = _areSecondariesConnected && pair.second->isConnected();
  }
  PRECICE_ASSERT(_areSecondariesConnected);
}

void M2N::requestSecondariesPreConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not _useOnlyPrimaryCom);
  _areSecondariesConnected = true;
  for (const auto &pair : _distComs) {
    pair.second->requestPreConnection(acceptorName, requesterName);
    _areSecondariesConnected = _areSecondariesConnected && pair.second->isConnected();
  }
  PRECICE_ASSERT(_areSecondariesConnected);
}

void M2N::completeSecondariesConnection()
{
  PRECICE_ASSERT(not _useOnlyPrimaryCom);
  for (const auto &pair : _distComs) {
    pair.second->completeSecondariesConnection();
  }
}

void M2N::closeConnection()
{
  PRECICE_TRACE();
  closePrimaryConnection();
  closeDistributedConnections();
}

void M2N::closePrimaryConnection()
{
  PRECICE_TRACE();
  if (not utils::IntraComm::isSecondary() && _masterCom->isConnected()) {
    _masterCom->closeConnection();
    _isPrimaryConnected = false;
  }

  utils::IntraComm::broadcast(_isPrimaryConnected);
  PRECICE_ASSERT(not _isPrimaryConnected);
}

void M2N::closeDistributedConnections()
{
  PRECICE_TRACE();
  if (_useOnlyPrimaryCom) {
    return;
  }

  _areSecondariesConnected = false;
  for (const auto &pair : _distComs) {
    pair.second->closeConnection();
    _areSecondariesConnected |= pair.second->isConnected();
  }
  PRECICE_ASSERT(not _areSecondariesConnected);
}

com::PtrCommunication M2N::getPrimaryCommunication()
{
  PRECICE_ASSERT(not utils::IntraComm::isSecondary());
  return _masterCom; /// @todo maybe it would be a nicer design to not offer this
}

void M2N::createDistributedCommunication(const mesh::PtrMesh &mesh)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(not _useOnlyPrimaryCom);
  DistributedCommunication::SharedPointer distCom = _distrFactory->newDistributedCommunication(mesh);
  _distComs[mesh->getID()]                        = distCom;
}

void M2N::send(
    precice::span<double const> itemsToSend,
    int                         meshID,
    int                         valueDimension)
{
  if (not _useOnlyPrimaryCom) {
    PRECICE_ASSERT(_areSecondariesConnected);
    PRECICE_ASSERT(_distComs.find(meshID) != _distComs.end());
    PRECICE_ASSERT(_distComs[meshID].get() != nullptr);

    if (precice::syncMode && not utils::IntraComm::isSecondary()) {
      bool ack = true;
      _masterCom->send(ack, 0);
      _masterCom->receive(ack, 0);
      _masterCom->send(ack, 0);
    }

    Event e("m2n.sendData", precice::syncMode);

    _distComs[meshID]->send(itemsToSend, valueDimension);
  } else {
    PRECICE_ASSERT(_isPrimaryConnected);
    _masterCom->send(itemsToSend, 0);
  }
}

void M2N::send(bool itemToSend)
{
  PRECICE_TRACE(utils::IntraComm::getRank());
  if (not utils::IntraComm::isSecondary()) {
    _masterCom->send(itemToSend, 0);
  }
}

void M2N::send(double itemToSend)
{
  PRECICE_TRACE(utils::IntraComm::getRank());
  if (not utils::IntraComm::isSecondary()) {
    _masterCom->send(itemToSend, 0);
  }
}

void M2N::broadcastSendMesh(mesh::Mesh &mesh)
{
  MeshID meshID = mesh.getID();
  PRECICE_ASSERT(utils::IntraComm::isParallel(),
                 "This method can only be used for parallel participants");
  PRECICE_ASSERT(_areSecondariesConnected);
  PRECICE_ASSERT(_distComs.find(meshID) != _distComs.end());
  PRECICE_ASSERT(_distComs[meshID].get() != nullptr);
  _distComs[meshID]->broadcastSendMesh();
}

void M2N::scatterAllCommunicationMap(std::map<int, std::vector<int>> &localCommunicationMap,
                                     mesh::Mesh &                     mesh)
{
  PRECICE_ASSERT(utils::IntraComm::isParallel(),
                 "This method can only be used for parallel participants");
  MeshID meshID = mesh.getID();
  PRECICE_ASSERT(_areSecondariesConnected);
  _distComs[meshID]->scatterAllCommunicationMap(localCommunicationMap);
}

void M2N::broadcastSend(int &itemToSend, mesh::Mesh &mesh)
{
  PRECICE_ASSERT(utils::IntraComm::isParallel(),
                 "This method can only be used for parallel participants");
  MeshID meshID = mesh.getID();
  PRECICE_ASSERT(_areSecondariesConnected);
  _distComs[meshID]->broadcastSend(itemToSend);
}

void M2N::receive(precice::span<double> itemsToReceive,
                  int                   meshID,
                  int                   valueDimension)
{
  if (not _useOnlyPrimaryCom) {
    PRECICE_ASSERT(_areSecondariesConnected);
    PRECICE_ASSERT(_distComs.find(meshID) != _distComs.end());
    PRECICE_ASSERT(_distComs[meshID].get() != nullptr);

    if (precice::syncMode) {
      if (not utils::IntraComm::isSecondary()) {
        bool ack;

        _masterCom->receive(ack, 0);
        _masterCom->send(ack, 0);
        _masterCom->receive(ack, 0);
      }
    }

    Event e("m2n.receiveData", precice::syncMode);

    _distComs[meshID]->receive(itemsToReceive, valueDimension);
  } else {
    PRECICE_ASSERT(_isPrimaryConnected);
    _masterCom->receive(itemsToReceive, 0);
  }
}

void M2N::receive(bool &itemToReceive)
{
  PRECICE_TRACE(utils::IntraComm::getRank());
  if (not utils::IntraComm::isSecondary()) {
    _masterCom->receive(itemToReceive, 0);
  }

  utils::IntraComm::broadcast(itemToReceive);

  PRECICE_DEBUG("receive(bool): {}", itemToReceive);
}

void M2N::receive(double &itemToReceive)
{
  PRECICE_TRACE(utils::IntraComm::getRank());
  if (not utils::IntraComm::isSecondary()) { //coupling mode
    _masterCom->receive(itemToReceive, 0);
  }

  utils::IntraComm::broadcast(itemToReceive);

  PRECICE_DEBUG("receive(double): {}", itemToReceive);
}

void M2N::broadcastReceiveAll(std::vector<int> &itemToReceive, mesh::Mesh &mesh)
{
  PRECICE_ASSERT(utils::IntraComm::isParallel(),
                 "This method can only be used for parallel participants");
  MeshID meshID = mesh.getID();
  PRECICE_ASSERT(_areSecondariesConnected);
  _distComs[meshID]->broadcastReceiveAll(itemToReceive);
}

void M2N::broadcastReceiveAllMesh(mesh::Mesh &mesh)
{
  PRECICE_ASSERT(utils::IntraComm::isParallel(),
                 "This method can only be used for parallel participants");
  MeshID meshID = mesh.getID();
  PRECICE_ASSERT(_areSecondariesConnected);
  PRECICE_ASSERT(_distComs.find(meshID) != _distComs.end());
  PRECICE_ASSERT(_distComs[meshID].get() != nullptr);
  _distComs[meshID]->broadcastReceiveAllMesh();
}

void M2N::gatherAllCommunicationMap(std::map<int, std::vector<int>> &localCommunicationMap, mesh::Mesh &mesh)
{
  PRECICE_ASSERT(utils::IntraComm::isParallel(),
                 "This method can only be used for parallel participants");
  MeshID meshID = mesh.getID();
  PRECICE_ASSERT(_areSecondariesConnected);
  _distComs[meshID]->gatherAllCommunicationMap(localCommunicationMap);
}

} // namespace m2n
} // namespace precice

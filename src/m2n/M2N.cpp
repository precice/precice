#include <utility>

#include "DistributedComFactory.hpp"
#include "DistributedCommunication.hpp"
#include "M2N.hpp"
#include "com/Communication.hpp"
#include "logging/LogConfiguration.hpp"
#include "logging/LogMacros.hpp"
#include "mesh/Mesh.hpp"
#include "precice/impl/Types.hpp"
#include "precice/impl/versions.hpp"
#include "profiling/Event.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

using precice::profiling::Event;

namespace precice {
extern bool syncMode;

namespace m2n {

M2N::M2N(com::PtrCommunication interComm, DistributedComFactory::SharedPointer distrFactory, bool useOnlyPrimaryCom, bool useTwoLevelInit)
    : _interComm(std::move(interComm)),
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
  return _isPrimaryRankConnected;
}

void M2N::acceptPrimaryRankConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);

  Event e("m2n.acceptPrimaryRankConnection." + requesterName, profiling::Fundamental, profiling::Synchronize);

  if (not utils::IntraComm::isSecondary()) {
    PRECICE_DEBUG("Accept primary connection");
    PRECICE_ASSERT(_interComm);
    _interComm->acceptConnection(acceptorName, requesterName, "PRIMARYCOM", utils::IntraComm::getRank());
    _isPrimaryRankConnected      = _interComm->isConnected();
    std::string acceptorVersion  = PRECICE_VERSION;
    auto        requesterVersion = acceptorVersion;
    _interComm->send(acceptorVersion, 0);
    _interComm->receive(requesterVersion, 0);
    PRECICE_CHECK(requesterVersion != acceptorVersion, "This participant {} uses preCICE version {} but the requester participant {} uses preCICE version {}.", acceptorName, acceptorVersion, requesterName, requesterVersion);
  }

  utils::IntraComm::broadcast(_isPrimaryRankConnected);
}

void M2N::requestPrimaryRankConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);

  Event e("m2n.requestPrimaryRankConnection." + acceptorName, profiling::Fundamental, profiling::Synchronize);

  if (not utils::IntraComm::isSecondary()) {
    PRECICE_ASSERT(_interComm);
    PRECICE_DEBUG("Request primary connection");
    _interComm->requestConnection(acceptorName, requesterName, "PRIMARYCOM", 0, 1);
    _isPrimaryRankConnected = _interComm->isConnected();

    //check that local and remote have the same precice version
    std::string requesterVersion = PRECICE_VERSION;
    auto        acceptorVersion  = requesterVersion;
    _interComm->receive(acceptorVersion, 0);
    _interComm->send(requesterVersion, 0);
    PRECICE_CHECK(requesterVersion != acceptorVersion, "This participant {} uses preCICE version {} but the acceptor participant {} uses preCICE version {}.", requesterName, requesterVersion, acceptorName, acceptorVersion);
  }

  utils::IntraComm::broadcast(_isPrimaryRankConnected);
}

void M2N::acceptSecondaryRanksConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not _useOnlyPrimaryCom);
  Event e("m2n.acceptSecondaryRanksConnection", profiling::Synchronize);

  _areSecondaryRanksConnected = true;
  for (const auto &pair : _distComs) {
    PRECICE_DEBUG("Accept secondary connections");
    pair.second->acceptConnection(acceptorName, requesterName);
    _areSecondaryRanksConnected = _areSecondaryRanksConnected && pair.second->isConnected();
  }
  PRECICE_ASSERT(_areSecondaryRanksConnected);
}

void M2N::requestSecondaryRanksConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not _useOnlyPrimaryCom);
  Event e("m2n.requestSecondaryRanksConnection", profiling::Synchronize);

  _areSecondaryRanksConnected = true;
  for (const auto &pair : _distComs) {
    PRECICE_DEBUG("Request secondary connections");
    pair.second->requestConnection(acceptorName, requesterName);
    _areSecondaryRanksConnected = _areSecondaryRanksConnected && pair.second->isConnected();
  }
  PRECICE_ASSERT(_areSecondaryRanksConnected);
}

void M2N::prepareEstablishment(const std::string &acceptorName,
                               const std::string &requesterName)
{
  PRECICE_TRACE();
  _interComm->prepareEstablishment(acceptorName, requesterName);
}

void M2N::cleanupEstablishment(const std::string &acceptorName,
                               const std::string &requesterName)
{
  PRECICE_TRACE();
  _interComm->cleanupEstablishment(acceptorName, requesterName);
}

void M2N::acceptSecondaryRanksPreConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not _useOnlyPrimaryCom);
  _areSecondaryRanksConnected = true;
  for (const auto &pair : _distComs) {
    pair.second->acceptPreConnection(acceptorName, requesterName);
    _areSecondaryRanksConnected = _areSecondaryRanksConnected && pair.second->isConnected();
  }
  PRECICE_ASSERT(_areSecondaryRanksConnected);
}

void M2N::requestSecondaryRanksPreConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not _useOnlyPrimaryCom);
  _areSecondaryRanksConnected = true;
  for (const auto &pair : _distComs) {
    pair.second->requestPreConnection(acceptorName, requesterName);
    _areSecondaryRanksConnected = _areSecondaryRanksConnected && pair.second->isConnected();
  }
  PRECICE_ASSERT(_areSecondaryRanksConnected);
}

void M2N::completeSecondaryRanksConnection()
{
  PRECICE_ASSERT(not _useOnlyPrimaryCom);
  for (const auto &pair : _distComs) {
    pair.second->completeSecondaryRanksConnection();
  }
}

void M2N::closeConnection()
{
  PRECICE_TRACE();
  closePrimaryRankConnection();
  closeDistributedConnections();
}

void M2N::closePrimaryRankConnection()
{
  PRECICE_TRACE();
  if (not utils::IntraComm::isSecondary() && _interComm->isConnected()) {
    _interComm->closeConnection();
    _isPrimaryRankConnected = false;
  }

  utils::IntraComm::broadcast(_isPrimaryRankConnected);
  PRECICE_ASSERT(not _isPrimaryRankConnected);
}

void M2N::closeDistributedConnections()
{
  PRECICE_TRACE();
  if (_useOnlyPrimaryCom) {
    return;
  }

  _areSecondaryRanksConnected = false;
  for (const auto &pair : _distComs) {
    pair.second->closeConnection();
    _areSecondaryRanksConnected |= pair.second->isConnected();
  }
  PRECICE_ASSERT(not _areSecondaryRanksConnected);
}

com::PtrCommunication M2N::getPrimaryRankCommunication()
{
  PRECICE_ASSERT(not utils::IntraComm::isSecondary());
  return _interComm; /// @todo maybe it would be a nicer design to not offer this
}

void M2N::createDistributedCommunication(const mesh::PtrMesh &mesh)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(not _useOnlyPrimaryCom);
  _distComs[mesh->getID()] = _distrFactory->newDistributedCommunication(mesh);
}

void M2N::send(
    precice::span<double const> itemsToSend,
    int                         meshID,
    int                         valueDimension)
{
  if (not _useOnlyPrimaryCom) {
    PRECICE_ASSERT(_areSecondaryRanksConnected);
    PRECICE_ASSERT(_distComs.find(meshID) != _distComs.end());
    PRECICE_ASSERT(_distComs[meshID].get() != nullptr);

    if (precice::syncMode && not utils::IntraComm::isSecondary()) {
      Event es("mn2.sendData.interSync");
      bool  ack = true;
      _interComm->send(ack, 0);
      _interComm->receive(ack, 0);
      _interComm->send(ack, 0);
    }

    Event e("m2n.sendData", profiling::Synchronize);

    _distComs[meshID]->send(itemsToSend, valueDimension);
  } else {
    PRECICE_ASSERT(_isPrimaryRankConnected);
    _interComm->send(itemsToSend, 0);
  }
}

void M2N::send(bool itemToSend)
{
  PRECICE_TRACE(utils::IntraComm::getRank());
  if (not utils::IntraComm::isSecondary()) {
    _interComm->send(itemToSend, 0);
  }
}

void M2N::send(double itemToSend)
{
  PRECICE_TRACE(utils::IntraComm::getRank());
  if (not utils::IntraComm::isSecondary()) {
    _interComm->send(itemToSend, 0);
  }
}

void M2N::send(precice::span<double const> itemsToSend)
{
  PRECICE_TRACE(utils::IntraComm::getRank());
  if (not utils::IntraComm::isSecondary()) {
    _interComm->send(itemsToSend, 0);
  }
}

void M2N::send(int itemToSend)
{
  PRECICE_TRACE(utils::IntraComm::getRank());
  if (not utils::IntraComm::isSecondary()) {
    _interComm->send(itemToSend, 0);
  }
}

void M2N::broadcastSendMesh(mesh::Mesh &mesh)
{
  MeshID meshID = mesh.getID();
  PRECICE_ASSERT(utils::IntraComm::isParallel(),
                 "This method can only be used for parallel participants");
  PRECICE_ASSERT(_areSecondaryRanksConnected);
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
  PRECICE_ASSERT(_areSecondaryRanksConnected);
  _distComs[meshID]->scatterAllCommunicationMap(localCommunicationMap);
}

void M2N::broadcastSend(int itemToSend, mesh::Mesh &mesh)
{
  PRECICE_ASSERT(utils::IntraComm::isParallel(),
                 "This method can only be used for parallel participants");
  MeshID meshID = mesh.getID();
  PRECICE_ASSERT(_areSecondaryRanksConnected);
  _distComs[meshID]->broadcastSend(itemToSend);
}

void M2N::receive(precice::span<double> itemsToReceive,
                  int                   meshID,
                  int                   valueDimension)
{
  if (not _useOnlyPrimaryCom) {
    PRECICE_ASSERT(_areSecondaryRanksConnected);
    PRECICE_ASSERT(_distComs.find(meshID) != _distComs.end());
    PRECICE_ASSERT(_distComs[meshID].get() != nullptr);

    if (precice::syncMode && not utils::IntraComm::isSecondary()) {
      Event es("m2n.receiveData.interSync");
      bool  ack;
      _interComm->receive(ack, 0);
      _interComm->send(ack, 0);
      _interComm->receive(ack, 0);
    }

    Event e("m2n.receiveData", profiling::Synchronize);

    _distComs[meshID]->receive(itemsToReceive, valueDimension);
  } else {
    PRECICE_ASSERT(_isPrimaryRankConnected);
    _interComm->receive(itemsToReceive, 0);
  }
}

void M2N::receive(bool &itemToReceive)
{
  PRECICE_TRACE(utils::IntraComm::getRank());
  if (not utils::IntraComm::isSecondary()) {
    _interComm->receive(itemToReceive, 0);
  }

  utils::IntraComm::broadcast(itemToReceive);

  PRECICE_DEBUG("receive(bool): {}", itemToReceive);
}

void M2N::receive(double &itemToReceive)
{
  PRECICE_TRACE(utils::IntraComm::getRank());
  if (not utils::IntraComm::isSecondary()) {
    _interComm->receive(itemToReceive, 0);
  }

  utils::IntraComm::broadcast(itemToReceive);

  PRECICE_DEBUG("receive(double): {}", itemToReceive);
}

void M2N::receive(precice::span<double> itemsToReceive)
{
  PRECICE_TRACE(utils::IntraComm::getRank());
  if (not utils::IntraComm::isSecondary()) {
    _interComm->receive(itemsToReceive, 0);
  }

  utils::IntraComm::broadcast(itemsToReceive);

  PRECICE_DEBUG("receive(span<double>) .size() = {}", itemsToReceive.size());
}

void M2N::receive(int &itemToReceive)
{
  PRECICE_TRACE(utils::IntraComm::getRank());
  if (not utils::IntraComm::isSecondary()) {
    _interComm->receive(itemToReceive, 0);
  }

  utils::IntraComm::broadcast(itemToReceive);
}

void M2N::broadcastReceiveAll(std::vector<int> &itemToReceive, mesh::Mesh &mesh)
{
  PRECICE_ASSERT(utils::IntraComm::isParallel(),
                 "This method can only be used for parallel participants");
  MeshID meshID = mesh.getID();
  PRECICE_ASSERT(_areSecondaryRanksConnected);
  _distComs[meshID]->broadcastReceiveAll(itemToReceive);
}

void M2N::broadcastReceiveAllMesh(mesh::Mesh &mesh)
{
  PRECICE_ASSERT(utils::IntraComm::isParallel(),
                 "This method can only be used for parallel participants");
  MeshID meshID = mesh.getID();
  PRECICE_ASSERT(_areSecondaryRanksConnected);
  PRECICE_ASSERT(_distComs.find(meshID) != _distComs.end());
  PRECICE_ASSERT(_distComs[meshID].get() != nullptr);
  _distComs[meshID]->broadcastReceiveAllMesh();
}

void M2N::gatherAllCommunicationMap(std::map<int, std::vector<int>> &localCommunicationMap, mesh::Mesh &mesh)
{
  PRECICE_ASSERT(utils::IntraComm::isParallel(),
                 "This method can only be used for parallel participants");
  MeshID meshID = mesh.getID();
  PRECICE_ASSERT(_areSecondaryRanksConnected);
  _distComs[meshID]->gatherAllCommunicationMap(localCommunicationMap);
}

} // namespace m2n
} // namespace precice

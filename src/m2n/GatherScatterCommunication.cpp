#include <algorithm>
#include <map>
#include <memory>
#include <ostream>
#include <utility>

#include "GatherScatterCommunication.hpp"
#include "com/Communication.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/DistributedCommunication.hpp"
#include "mesh/Mesh.hpp"
#include "precice/types.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace m2n {
GatherScatterCommunication::GatherScatterCommunication(
    com::PtrCommunication com,
    mesh::PtrMesh         mesh)
    : DistributedCommunication(std::move(mesh)),
      _com(std::move(com)),
      _isConnected(false)
{
}

GatherScatterCommunication::~GatherScatterCommunication()
{
  if (isConnected()) {
    closeConnection();
  }
}

bool GatherScatterCommunication::isConnected() const
{
  return _isConnected;
}

void GatherScatterCommunication::acceptConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(utils::MasterSlave::isSlave() || _com->isConnected());
  _isConnected = true;
}

void GatherScatterCommunication::requestConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(utils::MasterSlave::isSlave() || _com->isConnected());
  _isConnected = true;
}

void GatherScatterCommunication::closeConnection()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(utils::MasterSlave::isSlave() || not _com->isConnected());
  _isConnected = false;
}

void GatherScatterCommunication::send(precice::span<double const> itemsToSend, int valueDimension)
{
  PRECICE_TRACE(itemsToSend.size());

  // Gather data
  if (utils::MasterSlave::isSlave()) { // Slave
    if (!itemsToSend.empty()) {
      utils::MasterSlave::getCommunication()->send(itemsToSend, 0);
    }
  } else { // Master or coupling mode
    PRECICE_ASSERT(utils::MasterSlave::getRank() == 0);
    mesh::Mesh::VertexDistribution &vertexDistribution = _mesh->getVertexDistribution();
    int                             globalSize         = _mesh->getGlobalNumberOfVertices() * valueDimension;
    PRECICE_DEBUG("Global Size = {}", globalSize);
    std::vector<double> globalItemsToSend(globalSize);

    // Master data
    for (size_t i = 0; i < vertexDistribution[0].size(); i++) {
      for (int j = 0; j < valueDimension; j++) {
        globalItemsToSend[vertexDistribution[0][i] * valueDimension + j] += itemsToSend[i * valueDimension + j];
      }
    }

    // Slaves data
    for (Rank rankSlave : utils::MasterSlave::allSlaves()) {
      PRECICE_ASSERT(utils::MasterSlave::getCommunication() != nullptr);
      PRECICE_ASSERT(utils::MasterSlave::getCommunication()->isConnected());

      int slaveSize = vertexDistribution[rankSlave].size() * valueDimension;
      PRECICE_DEBUG("Slave Size = {}", slaveSize);
      if (slaveSize > 0) {
        std::vector<double> valuesSlave(slaveSize);
        utils::MasterSlave::getCommunication()->receive(valuesSlave, rankSlave);
        for (size_t i = 0; i < vertexDistribution[rankSlave].size(); i++) {
          for (int j = 0; j < valueDimension; j++) {
            globalItemsToSend[vertexDistribution[rankSlave][i] * valueDimension + j] += valuesSlave[i * valueDimension + j];
          }
        }
      }
    }

    // Send data to other master
    _com->send(globalItemsToSend, 0);
  }
}

void GatherScatterCommunication::receive(precice::span<double> itemsToReceive, int valueDimension)
{
  PRECICE_TRACE(itemsToReceive.size());

  std::vector<double> globalItemsToReceive;

  // Receive data at master
  if (not utils::MasterSlave::isSlave()) {
    int globalSize = _mesh->getGlobalNumberOfVertices() * valueDimension;
    PRECICE_DEBUG("Global Size = {}", globalSize);
    globalItemsToReceive.resize(globalSize);
    _com->receive(globalItemsToReceive, 0);
  }

  // Scatter data
  if (utils::MasterSlave::isSlave()) { // Slave
    if (!itemsToReceive.empty()) {
      PRECICE_DEBUG("itemsToRec[0] = {}", itemsToReceive[0]);
      utils::MasterSlave::getCommunication()->receive(itemsToReceive, 0);
      PRECICE_DEBUG("itemsToRec[0] = {}", itemsToReceive[0]);
    }
  } else { // Master or coupling mode
    PRECICE_ASSERT(utils::MasterSlave::getRank() == 0);
    mesh::Mesh::VertexDistribution &vertexDistribution = _mesh->getVertexDistribution();

    // Master data
    for (size_t i = 0; i < vertexDistribution[0].size(); i++) {
      for (int j = 0; j < valueDimension; j++) {
        itemsToReceive[i * valueDimension + j] = globalItemsToReceive[vertexDistribution[0][i] * valueDimension + j];
      }
    }

    // Slaves data
    for (Rank rankSlave : utils::MasterSlave::allSlaves()) {
      PRECICE_ASSERT(utils::MasterSlave::getCommunication() != nullptr);
      PRECICE_ASSERT(utils::MasterSlave::getCommunication()->isConnected());

      int slaveSize = vertexDistribution[rankSlave].size() * valueDimension;
      PRECICE_DEBUG("Slave Size = {}", slaveSize);
      if (slaveSize > 0) {
        std::vector<double> valuesSlave(slaveSize);
        for (size_t i = 0; i < vertexDistribution[rankSlave].size(); i++) {
          for (int j = 0; j < valueDimension; j++) {
            valuesSlave[i * valueDimension + j] = globalItemsToReceive[vertexDistribution[rankSlave][i] * valueDimension + j];
          }
        }
        utils::MasterSlave::getCommunication()->send(valuesSlave, rankSlave);
        PRECICE_DEBUG("valuesSlave[0] = {}", valuesSlave[0]);
      }
    }
  } // Master
}

void GatherScatterCommunication::acceptPreConnection(
    std::string const &acceptorName,
    std::string const &requesterName)
{
  PRECICE_ASSERT(false, "Not available for GatherScatterCommunication.");
}

void GatherScatterCommunication::requestPreConnection(
    std::string const &acceptorName,
    std::string const &requesterName)
{
  PRECICE_ASSERT(false, "Not available for GatherScatterCommunication.");
}

void GatherScatterCommunication::broadcastSend(const int &itemToSend)
{
  PRECICE_ASSERT(false, "Not available for GatherScatterCommunication.");
}

void GatherScatterCommunication::broadcastReceiveAll(std::vector<int> &itemToReceive)
{
  PRECICE_ASSERT(false, "Not available for GatherScatterCommunication.");
}

void GatherScatterCommunication::broadcastSendMesh()
{
  PRECICE_ASSERT(false, "Not available for GatherScatterCommunication.");
}

void GatherScatterCommunication::broadcastReceiveAllMesh()
{
  PRECICE_ASSERT(false, "Not available for GatherScatterCommunication.");
}

void GatherScatterCommunication::scatterAllCommunicationMap(CommunicationMap &localCommunicationMap)
{
  PRECICE_ASSERT(false, "Not available for GatherScatterCommunication.");
}

void GatherScatterCommunication::gatherAllCommunicationMap(CommunicationMap &localCommunicationMap)
{
  PRECICE_ASSERT(false, "Not available for GatherScatterCommunication.");
}

void GatherScatterCommunication::completeSlavesConnection()
{
  PRECICE_ASSERT(false, "Not available for GatherScatterCommunication.");
}

} // namespace m2n
} // namespace precice

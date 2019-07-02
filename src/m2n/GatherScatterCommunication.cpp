#include "GatherScatterCommunication.hpp"
#include "com/Communication.hpp"
#include "mesh/Mesh.hpp"
#include "utils/MasterSlave.hpp"

namespace precice
{
namespace m2n
{
GatherScatterCommunication::GatherScatterCommunication(
    com::PtrCommunication com,
    mesh::PtrMesh         mesh)
    : DistributedCommunication(mesh),
      _com(com),
      _isConnected(false)
{
}

GatherScatterCommunication::~GatherScatterCommunication()
{
  if (isConnected()) {
    closeConnection();
  }
}

bool GatherScatterCommunication::isConnected()
{
  return _isConnected;
}

void GatherScatterCommunication::acceptConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  TRACE(acceptorName, requesterName);
  assertion(utils::MasterSlave::isSlave() || _com->isConnected());
  _isConnected = true;
}

void GatherScatterCommunication::requestConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  TRACE(acceptorName, requesterName);
  assertion(utils::MasterSlave::isSlave() || _com->isConnected());
  _isConnected = true;
}

void GatherScatterCommunication::closeConnection()
{
  TRACE();
  assertion(utils::MasterSlave::isSlave() || not _com->isConnected());
  _isConnected = false;
}

void GatherScatterCommunication::send(
    double *itemsToSend,
    size_t  size,
    int     valueDimension)
{
  TRACE(size);
  assertion(utils::MasterSlave::isSlave() || utils::MasterSlave::isMaster());
  assertion(utils::MasterSlave::_communication.get() != nullptr);
  assertion(utils::MasterSlave::_communication->isConnected());
  assertion(utils::MasterSlave::getSize() > 1);
  assertion(utils::MasterSlave::getRank() != -1);

  // Gather data
  if (utils::MasterSlave::isSlave()) { // Slave
    if (size > 0) {
      utils::MasterSlave::_communication->send(itemsToSend, size, 0);
    }
  } else { // Master
    assertion(utils::MasterSlave::getRank() == 0);
    mesh::Mesh::VertexDistribution  &vertexDistribution = _mesh->getVertexDistribution();
    int                              globalSize         = _mesh->getGlobalNumberOfVertices() * valueDimension;
    DEBUG("Global Size = " << globalSize);
    std::vector<double> globalItemsToSend(globalSize);

    // Master data
    for (size_t i = 0; i < vertexDistribution[0].size(); i++) {
      for (int j = 0; j < valueDimension; j++) {
        globalItemsToSend[vertexDistribution[0][i] * valueDimension + j] += itemsToSend[i * valueDimension + j];
      }
    }

    // Slaves data
    for (int rankSlave = 1; rankSlave < utils::MasterSlave::getSize(); rankSlave++) {
      int slaveSize = vertexDistribution[rankSlave].size() * valueDimension;
      DEBUG("Slave Size = " << slaveSize);
      if (slaveSize > 0) {
        std::vector<double> valuesSlave(slaveSize);
        utils::MasterSlave::_communication->receive(valuesSlave.data(), slaveSize, rankSlave);
        for (size_t i = 0; i < vertexDistribution[rankSlave].size(); i++) {
          for (int j = 0; j < valueDimension; j++) {
            globalItemsToSend[vertexDistribution[rankSlave][i] * valueDimension + j] += valuesSlave[i * valueDimension + j];
          }
        }
      }
    }

    // Send data to other master
    _com->send(globalItemsToSend.data(), globalSize, 0);
  } // Master
}

void GatherScatterCommunication::receive(
    double *itemsToReceive,
    size_t  size,
    int     valueDimension)
{
  TRACE(size);
  assertion(utils::MasterSlave::isSlave() || utils::MasterSlave::isMaster());
  assertion(utils::MasterSlave::_communication.get() != nullptr);
  assertion(utils::MasterSlave::_communication->isConnected());
  assertion(utils::MasterSlave::getSize() > 1);
  assertion(utils::MasterSlave::getRank() != -1);

  std::vector<double> globalItemsToReceive;

  // Receive data at master
  if (utils::MasterSlave::isMaster()) {
    int globalSize = _mesh->getGlobalNumberOfVertices() * valueDimension;
    DEBUG("Global Size = " << globalSize);
    globalItemsToReceive.resize(globalSize);
    _com->receive(globalItemsToReceive.data(), globalSize, 0);
  }

  // Scatter data
  if (utils::MasterSlave::isSlave()) { // Slave
    if (size > 0) {
      DEBUG("itemsToRec[0] = " << itemsToReceive[0]);
      utils::MasterSlave::_communication->receive(itemsToReceive, size, 0);
      DEBUG("itemsToRec[0] = " << itemsToReceive[0]);
    }
  } else { // Master
    assertion(utils::MasterSlave::getRank() == 0);
    mesh::Mesh::VertexDistribution &vertexDistribution = _mesh->getVertexDistribution();

    // Master data
    for (size_t i = 0; i < vertexDistribution[0].size(); i++) {
      for (int j = 0; j < valueDimension; j++) {
        itemsToReceive[i * valueDimension + j] = globalItemsToReceive[vertexDistribution[0][i] * valueDimension + j];
      }
    }

    // Slaves data
    for (int rankSlave = 1; rankSlave < utils::MasterSlave::getSize(); rankSlave++) {
      int slaveSize = vertexDistribution[rankSlave].size() * valueDimension;
      DEBUG("Slave Size = " << slaveSize);
      if (slaveSize > 0) {
        std::vector<double> valuesSlave(slaveSize);
        for (size_t i = 0; i < vertexDistribution[rankSlave].size(); i++) {
          for (int j = 0; j < valueDimension; j++) {
            valuesSlave[i * valueDimension + j] = globalItemsToReceive[vertexDistribution[rankSlave][i] * valueDimension + j];
          }
        }
        utils::MasterSlave::_communication->send(valuesSlave.data(), slaveSize, rankSlave);
        DEBUG("valuesSlave[0] = " << valuesSlave[0]);
      }
    }
  } // Master
}

} // namespace m2n
} // namespace precice

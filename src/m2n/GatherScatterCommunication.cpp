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

void GatherScatterCommunication::send(
    double const *itemsToSend,
    size_t  size,
    int     valueDimension)
{
  PRECICE_TRACE(size);
  PRECICE_ASSERT(utils::MasterSlave::isSlave() || utils::MasterSlave::isMaster());
  PRECICE_ASSERT(utils::MasterSlave::_communication.get() != nullptr);
  PRECICE_ASSERT(utils::MasterSlave::_communication->isConnected());
  PRECICE_ASSERT(utils::MasterSlave::getSize() > 1);
  PRECICE_ASSERT(utils::MasterSlave::getRank() != -1);

  // Gather data
  if (utils::MasterSlave::isSlave()) { // Slave
    if (size > 0) {
      utils::MasterSlave::_communication->send(itemsToSend, size, 0);
    }
  } else { // Master
    PRECICE_ASSERT(utils::MasterSlave::getRank() == 0);
    mesh::Mesh::VertexDistribution  &vertexDistribution = _mesh->getVertexDistribution();
    int                              globalSize         = _mesh->getGlobalNumberOfVertices() * valueDimension;
    PRECICE_DEBUG("Global Size = " << globalSize);
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
      PRECICE_DEBUG("Slave Size = " << slaveSize);
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
  PRECICE_TRACE(size);
  PRECICE_ASSERT(utils::MasterSlave::isSlave() || utils::MasterSlave::isMaster());
  PRECICE_ASSERT(utils::MasterSlave::_communication.get() != nullptr);
  PRECICE_ASSERT(utils::MasterSlave::_communication->isConnected());
  PRECICE_ASSERT(utils::MasterSlave::getSize() > 1);
  PRECICE_ASSERT(utils::MasterSlave::getRank() != -1);

  std::vector<double> globalItemsToReceive;

  // Receive data at master
  if (utils::MasterSlave::isMaster()) {
    int globalSize = _mesh->getGlobalNumberOfVertices() * valueDimension;
    PRECICE_DEBUG("Global Size = " << globalSize);
    globalItemsToReceive.resize(globalSize);
    _com->receive(globalItemsToReceive.data(), globalSize, 0);
  }

  // Scatter data
  if (utils::MasterSlave::isSlave()) { // Slave
    if (size > 0) {
      PRECICE_DEBUG("itemsToRec[0] = " << itemsToReceive[0]);
      utils::MasterSlave::_communication->receive(itemsToReceive, size, 0);
      PRECICE_DEBUG("itemsToRec[0] = " << itemsToReceive[0]);
    }
  } else { // Master
    PRECICE_ASSERT(utils::MasterSlave::getRank() == 0);
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
      PRECICE_DEBUG("Slave Size = " << slaveSize);
      if (slaveSize > 0) {
        std::vector<double> valuesSlave(slaveSize);
        for (size_t i = 0; i < vertexDistribution[rankSlave].size(); i++) {
          for (int j = 0; j < valueDimension; j++) {
            valuesSlave[i * valueDimension + j] = globalItemsToReceive[vertexDistribution[rankSlave][i] * valueDimension + j];
          }
        }
        utils::MasterSlave::_communication->send(valuesSlave.data(), slaveSize, rankSlave);
        PRECICE_DEBUG("valuesSlave[0] = " << valuesSlave[0]);
      }
    }
  } // Master
}

void GatherScatterCommunication::acceptPreConnection(
  std::string const &acceptorName,
  std::string const &requesterName)
{  
  PRECICE_ASSERT(false, "This method can only be used with the point to point communication scheme");
}
 
void GatherScatterCommunication::requestPreConnection(
  std::string const &acceptorName,
  std::string const &requesterName)
{
  PRECICE_ASSERT(false, "This method can only be used with the point to point communication scheme");
}

void GatherScatterCommunication::broadcastSend(const int &itemToSend)
{
  PRECICE_ASSERT(false, "This method can only be used with the point to point communication scheme");
}

void GatherScatterCommunication::broadcastReceiveAll(std::vector<int> &itemToReceive)
{
  PRECICE_ASSERT(false, "This method can only be used with the point to point communication scheme");
}

void GatherScatterCommunication::broadcastSendMesh()
{
  PRECICE_ASSERT(false, "This method can only be used with the point to point communication scheme");
}

void GatherScatterCommunication::broadcastReceiveMesh()
{
  PRECICE_ASSERT(false, "This method can only be used with the point to point communication scheme");
}

void GatherScatterCommunication::broadcastSendLCM(CommunicationMap &localCommunicationMap)
{
  PRECICE_ASSERT(false, "This method can only be used with the point to point communication scheme");
}

void GatherScatterCommunication::broadcastReceiveLCM(CommunicationMap &localCommunicationMap)
{
  PRECICE_ASSERT(false, "This method can only be used with the point to point communication scheme");
}

void GatherScatterCommunication::updateVertexList()
{
  PRECICE_ASSERT(false, "This method can only be used with the point to point communication scheme");
}

} // namespace m2n
} // namespace precice

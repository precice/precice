#include "partition/ProvidedPartition.hpp"
#include "com/CommunicateMesh.hpp"
#include "com/Communication.hpp"
#include "m2n/M2N.hpp"
#include "utils/EventTimings.hpp"
#include "utils/MasterSlave.hpp"

using precice::utils::Event;

namespace precice
{
namespace partition
{

ProvidedPartition::ProvidedPartition(
    mesh::PtrMesh mesh,
    bool          hasToSend)
    : Partition(mesh),
      _hasToSend(hasToSend)
{
}

void ProvidedPartition::communicate()
{
  TRACE();

  //TODO communication to more than one participant

  if (_hasToSend) {
    Event e1("gather mesh");

    // Temporary globalMesh such that the master also keeps his local mesh
    mesh::Mesh globalMesh(_mesh->getName(), _mesh->getDimensions(), _mesh->isFlipNormals());

    if (not utils::MasterSlave::_slaveMode) {
      globalMesh.addMesh(*_mesh); // Add local master mesh to global mesh
    }

    // Gather Mesh
    INFO("Gather mesh " << _mesh->getName());
    if (utils::MasterSlave::_slaveMode or utils::MasterSlave::_masterMode) {
      if (utils::MasterSlave::_slaveMode) {
        com::CommunicateMesh(utils::MasterSlave::_communication).sendMesh(*_mesh, 0);
      } else { // Master
        assertion(utils::MasterSlave::_rank == 0);
        assertion(utils::MasterSlave::_size > 1);

        for (int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++) {
          com::CommunicateMesh(utils::MasterSlave::_communication).receiveMesh(globalMesh, rankSlave);
          DEBUG("Received sub-mesh, from slave: " << rankSlave << ", global vertexCount: " << globalMesh.vertices().size());
        }
      }
    }

    // Set global index
    if (not utils::MasterSlave::_slaveMode) {
      int globalIndex = 0;
      for (mesh::Vertex &v : globalMesh.vertices()) {
        v.setGlobalIndex(globalIndex);
        globalIndex++;
      }
    }

    e1.stop();

    // Send (global) Mesh
    INFO("Send global mesh " << _mesh->getName());
    Event e2("send global mesh");
    if (not utils::MasterSlave::_slaveMode) {
      CHECK(globalMesh.vertices().size() > 0, "The provided mesh " << globalMesh.getName() << " is invalid (possibly empty).");
      com::CommunicateMesh(_m2n->getMasterCommunication()).sendMesh(globalMesh, 0);
    }
    e2.stop();

  } //_hasToSend
}

void ProvidedPartition::compute()
{
  TRACE();
  INFO("Compute partition for mesh " << _mesh->getName());
  Event e6("feedback mesh");

  int numberOfVertices = _mesh->vertices().size();

  // Set global indices at every slave and vertexDistribution at master
  if (utils::MasterSlave::_slaveMode) {
    int globalVertexCounter = -1;
    utils::MasterSlave::_communication->send(numberOfVertices, 0);
    utils::MasterSlave::_communication->receive(globalVertexCounter, 0);
    for (int i = 0; i < numberOfVertices; i++) {
      _mesh->vertices()[i].setGlobalIndex(globalVertexCounter + i);
    }
    int globalNumberOfVertices = -1;
    utils::MasterSlave::_communication->broadcast(globalNumberOfVertices, 0);
    assertion(globalNumberOfVertices != -1);
    _mesh->setGlobalNumberOfVertices(globalNumberOfVertices);
  } else if (utils::MasterSlave::_masterMode) {
    assertion(utils::MasterSlave::_size > 1);
    int vertexCounter = 0;

    // Add master vertices
    for (int i = 0; i < numberOfVertices; i++) {
      _mesh->getVertexDistribution()[0].push_back(vertexCounter);
      _mesh->vertices()[i].setGlobalIndex(vertexCounter);
      vertexCounter++;
    }

    for (int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++) {
      utils::MasterSlave::_communication->receive(numberOfVertices, rankSlave);
      utils::MasterSlave::_communication->send(vertexCounter, rankSlave);

      for (int i = 0; i < numberOfVertices; i++) {
        _mesh->getVertexDistribution()[rankSlave].push_back(vertexCounter);
        vertexCounter++;
      }
    }
    _mesh->setGlobalNumberOfVertices(vertexCounter);
    utils::MasterSlave::_communication->broadcast(vertexCounter);
  } else { // Coupling mode
    for (int i = 0; i < numberOfVertices; i++) {
      _mesh->vertices()[i].setGlobalIndex(i);
    }
    _mesh->setGlobalNumberOfVertices(numberOfVertices);
  }

  createOwnerInformation();

  computeVertexOffsets();
}

void ProvidedPartition::createOwnerInformation()
{
  TRACE();
  for (mesh::Vertex &v : _mesh->vertices()) {
    v.setOwner(true);
  }
}

} // namespace partition
} // namespace precice

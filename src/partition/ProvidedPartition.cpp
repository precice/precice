#include "partition/ProvidedPartition.hpp"
#include "com/CommunicateMesh.hpp"
#include "com/Communication.hpp"
#include "m2n/M2N.hpp"
#include "utils/Event.hpp"
#include "utils/MasterSlave.hpp"

using precice::utils::Event;

namespace precice {
extern bool syncMode;

namespace partition {

ProvidedPartition::ProvidedPartition(
    mesh::PtrMesh mesh)
    : Partition(mesh)
{
}

void ProvidedPartition::communicate()
{
  PRECICE_TRACE();

  if(_m2ns.size()>0){ // if there is no connected participant we also don't need to gather the mesh
    Event e1("partition.gatherMesh." + _mesh->getName(), precice::syncMode);

    // Temporary globalMesh such that the master also keeps his local mesh
    mesh::Mesh globalMesh(_mesh->getName(), _mesh->getDimensions(), _mesh->isFlipNormals());

    if (not utils::MasterSlave::isSlave()) {
      globalMesh.addMesh(*_mesh); // Add local master mesh to global mesh
    }

    // Gather Mesh
    PRECICE_INFO("Gather mesh " + _mesh->getName());
    if (utils::MasterSlave::isSlave() ) {
        com::CommunicateMesh(utils::MasterSlave::_communication).sendMesh(*_mesh, 0);
    }
    if (utils::MasterSlave::isMaster())  {
      PRECICE_ASSERT(utils::MasterSlave::getRank() == 0);
      PRECICE_ASSERT(utils::MasterSlave::getSize() > 1);

      for (int rankSlave = 1; rankSlave < utils::MasterSlave::getSize(); rankSlave++) {
        com::CommunicateMesh(utils::MasterSlave::_communication).receiveMesh(globalMesh, rankSlave);
        PRECICE_DEBUG("Received sub-mesh, from slave: " << rankSlave << ", global vertexCount: " << globalMesh.vertices().size());
      }
    }

    // Set global index
    if (not utils::MasterSlave::isSlave()) {
      int globalIndex = 0;
      for (mesh::Vertex &v : globalMesh.vertices()) {
        v.setGlobalIndex(globalIndex);
        globalIndex++;
      }
    }

    e1.stop();

    // Send (global) Mesh
    PRECICE_INFO("Send global mesh " << _mesh->getName());
    Event e2("partition.sendGlobalMesh." + _mesh->getName(), precice::syncMode);

    for(auto m2n : _m2ns) {
      if (not utils::MasterSlave::isSlave()) {
        PRECICE_CHECK(globalMesh.vertices().size() > 0, "The provided mesh " << globalMesh.getName() << " is invalid (possibly empty).");
        com::CommunicateMesh(m2n->getMasterCommunication()).sendMesh(globalMesh, 0);
      }
    }
    e2.stop();

  }
}

void ProvidedPartition::compute()
{
  PRECICE_TRACE();
  PRECICE_INFO("Compute partition for mesh " << _mesh->getName());
  Event e6("partition.feedbackMesh." + _mesh->getName(), precice::syncMode);

  int numberOfVertices = _mesh->vertices().size();

  // Set global indices at every slave and vertexDistribution at master
  if (utils::MasterSlave::isSlave()) {
    int globalVertexCounter = -1;
    PRECICE_DEBUG("Send number of vertices: " << numberOfVertices);
    utils::MasterSlave::_communication->send(numberOfVertices, 0);
    utils::MasterSlave::_communication->receive(globalVertexCounter, 0);
    PRECICE_DEBUG("Set global vertex indices");
    for (int i = 0; i < numberOfVertices; i++) {
      _mesh->vertices()[i].setGlobalIndex(globalVertexCounter + i);
    }
    int globalNumberOfVertices = -1;
    utils::MasterSlave::_communication->broadcast(globalNumberOfVertices, 0);
    PRECICE_ASSERT(globalNumberOfVertices != -1);
    _mesh->setGlobalNumberOfVertices(globalNumberOfVertices);
  } else if (utils::MasterSlave::isMaster()) {
    PRECICE_ASSERT(utils::MasterSlave::getSize() > 1);
    int vertexCounter = 0;

    // Add master vertices
    PRECICE_DEBUG("Add master vertices to vertex distribution");
    for (int i = 0; i < numberOfVertices; i++) {
      _mesh->getVertexDistribution()[0].push_back(vertexCounter);
      _mesh->vertices()[i].setGlobalIndex(vertexCounter);
      vertexCounter++;
    }

    for (int rankSlave = 1; rankSlave < utils::MasterSlave::getSize(); rankSlave++) {
      utils::MasterSlave::_communication->receive(numberOfVertices, rankSlave);
      utils::MasterSlave::_communication->send(vertexCounter, rankSlave);

      for (int i = 0; i < numberOfVertices; i++) {
        _mesh->getVertexDistribution()[rankSlave].push_back(vertexCounter);
        vertexCounter++;
      }
    }
    _mesh->setGlobalNumberOfVertices(vertexCounter);
    PRECICE_DEBUG("broadcast global number of vertices: " << vertexCounter);
    utils::MasterSlave::_communication->broadcast(vertexCounter);
  } else { // Coupling mode
    for (int i = 0; i < numberOfVertices; i++) {
      _mesh->vertices()[i].setGlobalIndex(i);
    }
    _mesh->setGlobalNumberOfVertices(numberOfVertices);
  }

  PRECICE_DEBUG("Create owner information");
  createOwnerInformation();

  computeVertexOffsets();
}

void ProvidedPartition::createOwnerInformation()
{
  PRECICE_TRACE();
  for (mesh::Vertex &v : _mesh->vertices()) {
    v.setOwner(true);
  }
}

} // namespace partition
} // namespace precice

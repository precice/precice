#include "partition/ProvidedPartition.hpp"
#include <algorithm>
#include <map>
#include <memory>
#include <ostream>
#include <utility>
#include <vector>
#include "com/CommunicateBoundingBox.hpp"
#include "com/CommunicateMesh.hpp"
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/M2N.hpp"
#include "m2n/SharedPointer.hpp"
#include "mesh/BoundingBox.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "partition/Partition.hpp"
#include "utils/Event.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/assertion.hpp"

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

  prepare();

  if (_m2ns.empty())
    return;

  // Temporary globalMesh such that the master also keeps his local mesh
  mesh::Mesh globalMesh(_mesh->getName(), _mesh->getDimensions(), _mesh->isFlipNormals(), mesh::Mesh::MESH_ID_UNDEFINED);
  bool       hasMeshBeenGathered = false;

  bool twoLevelInitAlreadyUsed = false;

  for (auto &m2n : _m2ns) {
    if (m2n->usesTwoLevelInitialization()) {

      PRECICE_CHECK(not twoLevelInitAlreadyUsed, "Two-level initialization does not yet support multiple receivers of a provided mesh. "
                                                 "Please either switch two-level initialization off in your m2n definition, or "
                                                 "adapt your mesh setup such that each provided mesh is only received by maximum one "
                                                 "participant.");
      twoLevelInitAlreadyUsed = true;

      Event e("partition.broadcastMeshPartitions." + _mesh->getName(), precice::syncMode);

      // communicate the total number of vertices to the other participants master
      if (utils::MasterSlave::isMaster()) {
        _m2ns[0]->getMasterCommunication()->send(_mesh->getGlobalNumberOfVertices(), 0);
      }

      // the min and max of global vertex IDs of this rank's partition
      int minGlobalVertexID = _mesh->getVertexOffsets()[utils::MasterSlave::getRank()] - _mesh->vertices().size();
      int maxGlobalVertexID = _mesh->getVertexOffsets()[utils::MasterSlave::getRank()] - 1;

      // each rank sends its min/max global vertex index to connected remote ranks
      _m2ns[0]->broadcastSend(minGlobalVertexID, *_mesh);
      _m2ns[0]->broadcastSend(maxGlobalVertexID, *_mesh);

      // each rank sends its mesh partition to connected remote ranks
      _m2ns[0]->broadcastSendMesh(*_mesh);

    } else {

      if (not hasMeshBeenGathered) {
        //Gather mesh
        Event e("partition.gatherMesh." + _mesh->getName(), precice::syncMode);
        if (not utils::MasterSlave::isSlave()) {
          globalMesh.addMesh(*_mesh); // Add local master mesh to global mesh
        }
        PRECICE_INFO("Gather mesh " + _mesh->getName());
        if (utils::MasterSlave::isMaster()) {
          PRECICE_ASSERT(utils::MasterSlave::getRank() == 0);
          PRECICE_ASSERT(utils::MasterSlave::getSize() > 1);

          for (int rankSlave = 1; rankSlave < utils::MasterSlave::getSize(); rankSlave++) {
            com::CommunicateMesh(utils::MasterSlave::_communication).receiveMesh(globalMesh, rankSlave);
            PRECICE_DEBUG("Received sub-mesh, from slave: " << rankSlave << ", global vertexCount: " << globalMesh.vertices().size());
          }
        }
        if (utils::MasterSlave::isSlave()) {
          com::CommunicateMesh(utils::MasterSlave::_communication).sendMesh(*_mesh, 0);
        }
        hasMeshBeenGathered = true;
      }

      // Send (global) Mesh
      PRECICE_INFO("Send global mesh " << _mesh->getName());
      Event e("partition.sendGlobalMesh." + _mesh->getName(), precice::syncMode);

      if (not utils::MasterSlave::isSlave()) {
        PRECICE_CHECK(globalMesh.vertices().size() > 0,
                      "The provided mesh \"" << globalMesh.getName() << "\" is empty. Please set the mesh using setMeshXXX() prior to calling initialize().");
        com::CommunicateMesh(m2n->getMasterCommunication()).sendMesh(globalMesh, 0);
      }
    }
  }
}

void ProvidedPartition::prepare()
{
  PRECICE_TRACE();
  PRECICE_INFO("Prepare partition for mesh " << _mesh->getName());
  Event e("partition.prepareMesh." + _mesh->getName(), precice::syncMode);

  int numberOfVertices = _mesh->vertices().size();

  if (utils::MasterSlave::isMaster()) {
    PRECICE_ASSERT(utils::MasterSlave::getSize() > 1);

    // set globals IDs on master
    for (int i = 0; i < numberOfVertices; i++) {
      _mesh->vertices()[i].setGlobalIndex(i);
    }

    _mesh->getVertexOffsets().resize(utils::MasterSlave::getSize());
    _mesh->getVertexOffsets()[0] = numberOfVertices;
    int globalNumberOfVertices   = numberOfVertices;

    // receive number of slave vertices and fill vertex offsets
    for (int rankSlave = 1; rankSlave < utils::MasterSlave::getSize(); rankSlave++) {
      int numberOfSlaveVertices = -1;
      utils::MasterSlave::_communication->receive(numberOfSlaveVertices, rankSlave);
      _mesh->getVertexOffsets()[rankSlave] = numberOfSlaveVertices + _mesh->getVertexOffsets()[rankSlave - 1];
      utils::MasterSlave::_communication->send(globalNumberOfVertices, rankSlave);
      globalNumberOfVertices += numberOfSlaveVertices;
    }

    // set and broadcast global number of vertices
    _mesh->setGlobalNumberOfVertices(globalNumberOfVertices);
    PRECICE_DEBUG("Broadcast global number of vertices: " << globalNumberOfVertices);
    utils::MasterSlave::_communication->broadcast(globalNumberOfVertices);

    // broadcast vertex offsets
    PRECICE_DEBUG("My vertex offsets: " << _mesh->getVertexOffsets());
    utils::MasterSlave::_communication->broadcast(_mesh->getVertexOffsets());

    // fill vertex distribution
    if (std::any_of(_m2ns.begin(), _m2ns.end(), [](const m2n::PtrM2N &m2n) { return not m2n->usesTwoLevelInitialization(); })) {
      if (utils::MasterSlave::isMaster()) {
        PRECICE_DEBUG("Fill vertex distribution");
        auto &localIds = _mesh->getVertexDistribution()[0];
        for (int i = 0; i < _mesh->getVertexOffsets()[0]; i++) {
          localIds.push_back(i);
        }
        for (int rankSlave = 1; rankSlave < utils::MasterSlave::getSize(); rankSlave++) {
          // This always creates an entry for each slave
          auto &slaveIds = _mesh->getVertexDistribution()[rankSlave];
          for (int i = _mesh->getVertexOffsets()[rankSlave - 1]; i < _mesh->getVertexOffsets()[rankSlave]; i++) {
            slaveIds.push_back(i);
          }
        }
        PRECICE_ASSERT(_mesh->getVertexDistribution().size() == static_cast<decltype(_mesh->getVertexDistribution().size())>(utils::MasterSlave::getSize()));
      }
    }
  } else if (utils::MasterSlave::isSlave()) {

    // send number of own vertices
    PRECICE_DEBUG("Send number of vertices: " << numberOfVertices);
    utils::MasterSlave::_communication->send(numberOfVertices, 0);

    // set global IDs
    int globalVertexCounter = -1;
    utils::MasterSlave::_communication->receive(globalVertexCounter, 0);
    PRECICE_DEBUG("Set global vertex indices");
    for (int i = 0; i < numberOfVertices; i++) {
      _mesh->vertices()[i].setGlobalIndex(globalVertexCounter + i);
    }

    // set global number of vertices
    int globalNumberOfVertices = -1;
    utils::MasterSlave::_communication->broadcast(globalNumberOfVertices, 0);
    PRECICE_ASSERT(globalNumberOfVertices != -1);
    _mesh->setGlobalNumberOfVertices(globalNumberOfVertices);

    // set vertex offsets
    utils::MasterSlave::_communication->broadcast(_mesh->getVertexOffsets(), 0);
    PRECICE_DEBUG("My vertex offsets: " << _mesh->getVertexOffsets());

  } else { // Coupling mode

    for (int i = 0; i < numberOfVertices; i++) {
      _mesh->getVertexDistribution()[0].push_back(i);
      _mesh->vertices()[i].setGlobalIndex(i);
    }
    _mesh->getVertexOffsets().push_back(numberOfVertices);
    _mesh->setGlobalNumberOfVertices(numberOfVertices);
  }

  PRECICE_DEBUG("Set owner information");
  for (mesh::Vertex &v : _mesh->vertices()) {
    v.setOwner(true);
  }
}

void ProvidedPartition::compute()
{
  PRECICE_TRACE();
  for (auto m2n : _m2ns) {
    if (m2n->usesTwoLevelInitialization()) {
      // @todo this will probably not work for more than one m2n
      PRECICE_ASSERT(_m2ns.size() <= 1);
      // receive communication map from all remote connected ranks
      m2n->gatherAllCommunicationMap(_mesh->getCommunicationMap(), *_mesh);
    }
  }
}

void ProvidedPartition::compareBoundingBoxes()
{
  PRECICE_TRACE();
  if (_m2ns.empty())
    return;

  //@todo coupling mode

  //@todo treatment of multiple m2ns
  if (not _m2ns[0]->usesTwoLevelInitialization())
    return;

  // each rank sends its bb to master
  if (utils::MasterSlave::isSlave()) { //slave
    PRECICE_ASSERT(_mesh->getBoundingBox().getDimension() == _mesh->getDimensions(), "The boundingbox of the local mesh is invalid!");
    com::CommunicateBoundingBox(utils::MasterSlave::_communication).sendBoundingBox(_mesh->getBoundingBox(), 0);
  } else { // Master

    PRECICE_ASSERT(utils::MasterSlave::getRank() == 0);
    PRECICE_ASSERT(utils::MasterSlave::getSize() > 1);

    // to store the collection of bounding boxes
    mesh::Mesh::BoundingBoxMap bbm;
    mesh::BoundingBox          bb(_mesh->getDimensions());
    bbm.emplace(0, _mesh->getBoundingBox());
    PRECICE_ASSERT(!bbm.empty(), "The bounding box of the local mesh is invalid!");

    // master receives bbs from slaves and stores them in bbm
    for (int rankSlave = 1; rankSlave < utils::MasterSlave::getSize(); rankSlave++) {
      // initialize bbm
      bbm.emplace(rankSlave, bb);
      com::CommunicateBoundingBox(utils::MasterSlave::_communication).receiveBoundingBox(bbm.at(rankSlave), rankSlave);
    }

    // master sends number of ranks and bbm to the other master
    _m2ns[0]->getMasterCommunication()->send(utils::MasterSlave::getSize(), 0);
    com::CommunicateBoundingBox(_m2ns[0]->getMasterCommunication()).sendBoundingBoxMap(bbm, 0);
  }

  // size of the feedbackmap
  int              remoteConnectionMapSize = 0;
  std::vector<int> connectedRanksList;

  std::map<int, std::vector<int>> remoteConnectionMap;

  if (utils::MasterSlave::isMaster()) {

    // master receives feedback map (map of other participant ranks -> connected ranks at this participant)
    // from other participants master
    _m2ns[0]->getMasterCommunication()->receive(connectedRanksList, 0);
    remoteConnectionMapSize = connectedRanksList.size();

    for (auto &rank : connectedRanksList) {
      remoteConnectionMap[rank] = {-1};
    }
    if (remoteConnectionMapSize != 0) {
      com::CommunicateBoundingBox(_m2ns[0]->getMasterCommunication()).receiveConnectionMap(remoteConnectionMap, 0);
    }

    // broadcast the received feedbackMap
    utils::MasterSlave::_communication->broadcast(connectedRanksList);
    if (remoteConnectionMapSize != 0) {
      com::CommunicateBoundingBox(utils::MasterSlave::_communication).broadcastSendConnectionMap(remoteConnectionMap);
    }

    // master checks which ranks are connected to it
    for (auto &remoteRank : remoteConnectionMap) {
      for (auto &includedRank : remoteRank.second) {
        if (utils::MasterSlave::getRank() == includedRank) {
          _mesh->getConnectedRanks().push_back(remoteRank.first);
        }
      }
    }

  } else { // Slave

    utils::MasterSlave::_communication->broadcast(connectedRanksList, 0);

    if (!connectedRanksList.empty()) {
      for (int rank : connectedRanksList) {
        remoteConnectionMap[rank] = {-1};
      }
      com::CommunicateBoundingBox(utils::MasterSlave::_communication).broadcastReceiveConnectionMap(remoteConnectionMap);
    }

    for (const auto &remoteRank : remoteConnectionMap) {
      for (int includedRanks : remoteRank.second) {
        if (utils::MasterSlave::getRank() == includedRanks) {
          _mesh->getConnectedRanks().push_back(remoteRank.first);
        }
      }
    }
  }
}

} // namespace partition
} // namespace precice

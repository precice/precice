#include <algorithm>
#include <map>
#include <memory>
#include <numeric>
#include <ostream>
#include <utility>
#include <vector>

#include "com/Communication.hpp"
#include "com/Extra.hpp"
#include "com/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/M2N.hpp"
#include "m2n/SharedPointer.hpp"
#include "mesh/BoundingBox.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "partition/Partition.hpp"
#include "partition/ProvidedPartition.hpp"
#include "precice/types.hpp"
#include "profiling/Event.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

using precice::profiling::Event;

namespace precice::partition {

ProvidedPartition::ProvidedPartition(
    mesh::PtrMesh mesh)
    : Partition(std::move(mesh))
{
}

void ProvidedPartition::communicate()
{
  PRECICE_TRACE();

  prepare();

  if (_m2ns.empty())
    return;

  // Temporary globalMesh such that the primary rank also keeps his local mesh
  mesh::Mesh globalMesh(_mesh->getName(), _mesh->getDimensions(), mesh::Mesh::MESH_ID_UNDEFINED);
  bool       hasMeshBeenGathered = false;

  bool twoLevelInitAlreadyUsed = false;

  for (auto &m2n : _m2ns) {
    if (m2n->usesTwoLevelInitialization()) {

      PRECICE_CHECK(not twoLevelInitAlreadyUsed, "Two-level initialization does not yet support multiple receivers of a provided mesh. "
                                                 "Please either switch two-level initialization off in your m2n definition, or "
                                                 "adapt your mesh setup such that each provided mesh is only received by maximum one "
                                                 "participant.");
      twoLevelInitAlreadyUsed = true;

      Event e("partition.broadcastMeshPartitions." + _mesh->getName(), profiling::Synchronize);

      // communicate the total number of vertices to the other participants primary rank
      if (utils::IntraComm::isPrimary()) {
        _m2ns[0]->getPrimaryRankCommunication()->send(_mesh->getGlobalNumberOfVertices(), 0);
      }

      // the min and max of global vertex IDs of this rank's partition
      PRECICE_ASSERT(_mesh->getVertexOffsets().size() == static_cast<decltype(_mesh->getVertexOffsets().size())>(utils::IntraComm::getSize()));
      const int vertexOffset      = _mesh->getVertexOffsets()[utils::IntraComm::getRank()];
      const int minGlobalVertexID = vertexOffset - _mesh->vertices().size();
      const int maxGlobalVertexID = vertexOffset - 1;

      // each rank sends its min/max global vertex index to connected remote ranks
      _m2ns[0]->broadcastSend(minGlobalVertexID, *_mesh);
      _m2ns[0]->broadcastSend(maxGlobalVertexID, *_mesh);

      // each rank sends its mesh partition to connected remote ranks
      _m2ns[0]->broadcastSendMesh(*_mesh);

    } else {

      if (not hasMeshBeenGathered) {
        // Gather mesh
        Event e("partition.gatherMesh." + _mesh->getName(), profiling::Synchronize);
        if (not utils::IntraComm::isSecondary()) {
          globalMesh.addMesh(*_mesh); // Add local primary mesh to global mesh
        }
        PRECICE_INFO("Gather mesh {}", _mesh->getName());
        if (utils::IntraComm::isPrimary()) {
          PRECICE_ASSERT(utils::IntraComm::getRank() == 0);
          PRECICE_ASSERT(utils::IntraComm::getSize() > 1);

          for (Rank secondaryRank : utils::IntraComm::allSecondaryRanks()) {
            com::receiveMesh(*utils::IntraComm::getCommunication(), secondaryRank, globalMesh);
            PRECICE_DEBUG("Received sub-mesh, from secondary rank: {}, global vertexCount: {}", secondaryRank, globalMesh.vertices().size());
          }
        }
        if (utils::IntraComm::isSecondary()) {
          com::sendMesh(*utils::IntraComm::getCommunication(), 0, *_mesh);
        }
        hasMeshBeenGathered = true;
      }

      // Send (global) Mesh
      PRECICE_INFO("Send global mesh {}", _mesh->getName());
      Event e("partition.sendGlobalMesh." + _mesh->getName(), profiling::Synchronize);

      if (not utils::IntraComm::isSecondary()) {
        PRECICE_CHECK(globalMesh.vertices().size() > 0,
                      "The provided mesh \"{}\" is empty. Please set the mesh using setMeshXXX() prior to calling initialize().",
                      globalMesh.getName());
        com::sendMesh(*m2n->getPrimaryRankCommunication(), 0, globalMesh);
      }
    }
  }
}

void ProvidedPartition::prepare()
{
  PRECICE_TRACE();
  PRECICE_INFO("Prepare partition for mesh {}", _mesh->getName());
  Event e("partition.prepareMesh." + _mesh->getName(), profiling::Synchronize);

  int numberOfVertices = _mesh->vertices().size();

  if (utils::IntraComm::isPrimary()) {
    PRECICE_ASSERT(utils::IntraComm::getSize() > 1);

    // set globals IDs on primary rank
    for (int i = 0; i < numberOfVertices; i++) {
      _mesh->vertices()[i].setGlobalIndex(i);
    }

    mesh::Mesh::VertexOffsets vertexOffsets(utils::IntraComm::getSize());
    vertexOffsets[0]           = numberOfVertices;
    int globalNumberOfVertices = numberOfVertices;

    // receive number of secondary vertices and fill vertex offsets
    for (Rank secondaryRank : utils::IntraComm::allSecondaryRanks()) {
      int numberOfSecondaryRankVertices = -1;
      utils::IntraComm::getCommunication()->receive(numberOfSecondaryRankVertices, secondaryRank);
      vertexOffsets[secondaryRank] = numberOfSecondaryRankVertices + vertexOffsets[secondaryRank - 1];
      utils::IntraComm::getCommunication()->send(globalNumberOfVertices, secondaryRank);
      globalNumberOfVertices += numberOfSecondaryRankVertices;
    }
    PRECICE_ASSERT(std::all_of(vertexOffsets.begin(), vertexOffsets.end(), [](auto i) { return i >= 0; }));
    PRECICE_ASSERT(_mesh->getVertexOffsets().empty());
    _mesh->setVertexOffsets(vertexOffsets);

    // set and broadcast global number of vertices
    _mesh->setGlobalNumberOfVertices(globalNumberOfVertices);
    PRECICE_DEBUG("Broadcast global number of vertices: {}", globalNumberOfVertices);
    utils::IntraComm::getCommunication()->broadcast(globalNumberOfVertices);

    // broadcast vertex offsets to secondary ranks
    PRECICE_DEBUG("My vertex offsets: {}", vertexOffsets);
    utils::IntraComm::getCommunication()->broadcast(vertexOffsets);

    // fill vertex distribution
    if (std::any_of(_m2ns.begin(), _m2ns.end(), [](const m2n::PtrM2N &m2n) { return not m2n->usesTwoLevelInitialization(); }) && utils::IntraComm::isPrimary()) {
      PRECICE_DEBUG("Fill vertex distribution");
      PRECICE_ASSERT(_mesh->getVertexDistribution().empty());
      /// @TODO are these distributions allowed to contain verices already?
      mesh::Mesh::VertexDistribution vertexDistribution;
      auto &                         localIds = vertexDistribution[0];
      localIds.resize(vertexOffsets[0]);
      std::iota(localIds.begin(), localIds.end(), 0);

      for (Rank secondaryRank : utils::IntraComm::allSecondaryRanks()) {
        // This always creates an entry for each secondary rank
        auto &secondaryIds = vertexDistribution[secondaryRank];
        for (int i = vertexOffsets[secondaryRank - 1]; i < vertexOffsets[secondaryRank]; i++) {
          secondaryIds.push_back(i);
        }
      }
      PRECICE_ASSERT(vertexDistribution.size() == static_cast<mesh::Mesh::VertexDistribution::size_type>(utils::IntraComm::getSize()));
      _mesh->setVertexDistribution(std::move(vertexDistribution));
    }
  } else if (utils::IntraComm::isSecondary()) {

    // send number of own vertices
    PRECICE_DEBUG("Send number of vertices: {}", numberOfVertices);
    utils::IntraComm::getCommunication()->send(numberOfVertices, 0);

    // set global IDs
    int globalVertexCounter = -1;
    utils::IntraComm::getCommunication()->receive(globalVertexCounter, 0);
    PRECICE_DEBUG("Set global vertex indices");
    for (int i = 0; i < numberOfVertices; i++) {
      _mesh->vertices()[i].setGlobalIndex(globalVertexCounter + i);
    }

    // set global number of vertices
    int globalNumberOfVertices = -1;
    utils::IntraComm::getCommunication()->broadcast(globalNumberOfVertices, 0);
    PRECICE_ASSERT(globalNumberOfVertices != -1);
    _mesh->setGlobalNumberOfVertices(globalNumberOfVertices);

    // receive set vertex offsets
    mesh::Mesh::VertexOffsets vertexOffsets;
    utils::IntraComm::getCommunication()->broadcast(vertexOffsets, 0);
    PRECICE_DEBUG("My vertex offsets: {}", vertexOffsets);
    PRECICE_ASSERT(_mesh->getVertexOffsets().empty());
    _mesh->setVertexOffsets(std::move(vertexOffsets));

  } else {

    // The only rank of the participant contains all vertices
    PRECICE_ASSERT(_mesh->getVertexDistribution().empty());
    _mesh->setVertexDistribution([&] {
      mesh::Mesh::VertexDistribution vertexDistribution;
      for (int i = 0; i < numberOfVertices; i++) {
        vertexDistribution[0].push_back(i);
        _mesh->vertices()[i].setGlobalIndex(i);
      }
      return vertexDistribution;
    }());
    PRECICE_ASSERT(_mesh->getVertexOffsets().empty());
    _mesh->setVertexOffsets({numberOfVertices});
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
  for (const auto &m2n : _m2ns) {
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

  _mesh->clearPartitioning();

  //@todo coupling mode

  //@todo treatment of multiple m2ns
  if (not _m2ns[0]->usesTwoLevelInitialization())
    return;

  // each secondary rank sends its bb to the primary rank
  if (utils::IntraComm::isSecondary()) { // secondary
    PRECICE_ASSERT(_mesh->getBoundingBox().getDimension() == _mesh->getDimensions(), "The boundingbox of the local mesh is invalid!");
    com::sendBoundingBox(*utils::IntraComm::getCommunication(), 0, _mesh->getBoundingBox());
  } else { // Primary

    PRECICE_ASSERT(utils::IntraComm::getRank() == 0);
    PRECICE_ASSERT(utils::IntraComm::getSize() > 1);

    // to store the collection of bounding boxes
    mesh::Mesh::BoundingBoxMap bbm;
    mesh::BoundingBox          bb(_mesh->getDimensions());
    bbm.emplace(0, _mesh->getBoundingBox());
    PRECICE_ASSERT(!bbm.empty(), "The bounding box of the local mesh is invalid!");

    // primary rank receives bbs from secondary ranks and stores them in bbm
    for (Rank secondaryRank : utils::IntraComm::allSecondaryRanks()) {
      // initialize bbm
      bbm.emplace(secondaryRank, bb);
      com::receiveBoundingBox(*utils::IntraComm::getCommunication(), secondaryRank, bbm.at(secondaryRank));
    }

    // primary rank sends number of ranks and bbm to the other primary rank
    _m2ns[0]->getPrimaryRankCommunication()->send(utils::IntraComm::getSize(), 0);
    com::sendBoundingBoxMap(*_m2ns[0]->getPrimaryRankCommunication(), 0, bbm);
  }

  // size of the feedbackmap
  int remoteConnectionMapSize = 0;

  if (utils::IntraComm::isPrimary()) {

    // primary rank receives feedback map (map of other participant ranks -> connected ranks at this participant)
    // from other participants primary rank
    std::vector<Rank> connectedRanksList = _m2ns[0]->getPrimaryRankCommunication()->receiveRange(0, com::AsVectorTag<Rank>{});
    remoteConnectionMapSize              = connectedRanksList.size();

    mesh::Mesh::CommunicationMap remoteConnectionMap;
    for (auto &rank : connectedRanksList) {
      remoteConnectionMap[rank] = {-1};
    }
    if (remoteConnectionMapSize != 0) {
      com::receiveConnectionMap(*_m2ns[0]->getPrimaryRankCommunication(), 0, remoteConnectionMap);
    }

    // broadcast the received feedbackMap
    utils::IntraComm::getCommunication()->broadcast(connectedRanksList);
    if (remoteConnectionMapSize != 0) {
      com::broadcastSendConnectionMap(*utils::IntraComm::getCommunication(), remoteConnectionMap);
    }

    // primary rank checks which ranks are connected to it
    PRECICE_ASSERT(_mesh->getConnectedRanks().empty());
    _mesh->setConnectedRanks([&] {
      std::vector<Rank> ranks;
      for (const auto &remoteRank : remoteConnectionMap) {
        for (const auto &includedRank : remoteRank.second) {
          if (utils::IntraComm::getRank() == includedRank) {
            ranks.push_back(remoteRank.first);
          }
        }
      }
      return ranks;
    }());

  } else { // Secondary rank
    std::vector<Rank> connectedRanksList;
    utils::IntraComm::getCommunication()->broadcast(connectedRanksList, 0);

    mesh::Mesh::CommunicationMap remoteConnectionMap;
    if (!connectedRanksList.empty()) {
      for (Rank rank : connectedRanksList) {
        remoteConnectionMap[rank] = {-1};
      }
      com::broadcastReceiveConnectionMap(*utils::IntraComm::getCommunication(), remoteConnectionMap);
    }

    PRECICE_ASSERT(_mesh->getConnectedRanks().empty());
    _mesh->setConnectedRanks([&] {
      std::vector<Rank> ranks;
      for (const auto &remoteRank : remoteConnectionMap) {
        for (int includedRanks : remoteRank.second) {
          if (utils::IntraComm::getRank() == includedRanks) {
            ranks.push_back(remoteRank.first);
          }
        }
      }
      return ranks;
    }());
  }
}

} // namespace precice::partition

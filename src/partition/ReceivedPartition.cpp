#include "partition/ReceivedPartition.hpp"
#include "com/CommunicateMesh.hpp"
#include "com/Communication.hpp"
#include "m2n/M2N.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"
#include "utils/Event.hpp"
#include "utils/Helpers.hpp"
#include "utils/MasterSlave.hpp"

using precice::utils::Event;

namespace precice {
extern bool syncMode;

namespace partition {

ReceivedPartition::ReceivedPartition(
    mesh::PtrMesh mesh, GeometricFilter geometricFilter, double safetyFactor)
    : Partition(mesh),
      _geometricFilter(geometricFilter),
      _bb(mesh->getDimensions(), std::make_pair(std::numeric_limits<double>::max(),
                                                std::numeric_limits<double>::lowest())),
      _dimensions(mesh->getDimensions()),
      _safetyFactor(safetyFactor)
{
}

void ReceivedPartition::communicate()
{
  PRECICE_TRACE();
  PRECICE_INFO("Receive global mesh " << _mesh->getName());
  Event e("partition.receiveGlobalMesh." + _mesh->getName(), precice::syncMode);
  if (not utils::MasterSlave::isSlave()) {
    PRECICE_ASSERT(_mesh->vertices().empty());
    // a ReceivedPartition can only have one communication, @todo nicer design
    com::CommunicateMesh(_m2ns[0]->getMasterCommunication()).receiveMesh(*_mesh, 0);
  }
}

void ReceivedPartition::compute()
{
  PRECICE_TRACE(_geometricFilter);

  // handle coupling mode first (i.e. serial participant)
  if (not utils::MasterSlave::isSlave() && not utils::MasterSlave::isMaster()) { //coupling mode
    PRECICE_DEBUG("Handle partition data structures for serial participant");
    _mesh->setGlobalNumberOfVertices(_mesh->vertices().size());
    computeVertexOffsets();
    for (mesh::Vertex &v : _mesh->vertices()) {
      v.setOwner(true);
    }
    return;
  }

  // check to prevent false configuration
  if (not utils::MasterSlave::isSlave()) {
    PRECICE_CHECK(_fromMapping || _toMapping,
          "The received mesh " << _mesh->getName()
          << " needs a mapping, either from it, to it, or both. Maybe you don't want to receive this mesh at all?")
  }


  // To understand the following steps, it is recommended to look at BU's thesis, especially Figure 69 on page 89
  // for RBF-based filtering. https://mediatum.ub.tum.de/doc/1320661/document.pdf


  // (0) set global number of vertices before filtering
  if (utils::MasterSlave::isMaster()) {
    PRECICE_DEBUG("Set global number of vertices");
    _mesh->setGlobalNumberOfVertices(_mesh->vertices().size());
  }

  // (1) Bounding-Box-Filter

  if (_geometricFilter == FILTER_FIRST) { //pre-filter-post-filter

    PRECICE_INFO("Pre-filter mesh " << _mesh->getName() << " by bounding-box");
    Event e("partition.preFilterMesh." + _mesh->getName(), precice::syncMode);

    if (utils::MasterSlave::isSlave()) {
      PRECICE_DEBUG("Send bounding box to master");
      prepareBoundingBox();
      com::CommunicateMesh(utils::MasterSlave::_communication).sendBoundingBox(_bb, 0);
      PRECICE_DEBUG("Receive filtered mesh");
      com::CommunicateMesh(utils::MasterSlave::_communication).receiveMesh(*_mesh, 0);

      if(areProvidedMeshesEmpty()) {
        std::string msg = "The re-partitioning completely filtered out the mesh " + _mesh->getName() +
          " received on this rank at the coupling interface. "
          "Most probably, the coupling interfaces of your coupled participants do not match geometry-wise. "
          "Please check your geometry setup again. Small overlaps or gaps are no problem. "
          "If your geometry setup is correct and if you have very different mesh resolutions on both sides, increasing the safety-factor "
          "of the decomposition strategy might be necessary.";
        PRECICE_CHECK(not _mesh->vertices().empty(), msg);
      }

    } else { // Master
      PRECICE_ASSERT(utils::MasterSlave::getRank() == 0);
      PRECICE_ASSERT(utils::MasterSlave::getSize() > 1);

      for (int rankSlave = 1; rankSlave < utils::MasterSlave::getSize(); rankSlave++) {
        com::CommunicateMesh(utils::MasterSlave::_communication).receiveBoundingBox(_bb, rankSlave);

        PRECICE_DEBUG("From slave " << rankSlave << ", bounding mesh: " << _bb[0].first
              << ", " << _bb[0].second << " and " << _bb[1].first << ", " << _bb[1].second);
        mesh::Mesh slaveMesh("SlaveMesh", _dimensions, _mesh->isFlipNormals());
        filterMesh(slaveMesh, true);
        PRECICE_DEBUG("Send filtered mesh to slave: " << rankSlave);
        com::CommunicateMesh(utils::MasterSlave::_communication).sendMesh(slaveMesh, rankSlave);
      }

      // Now also filter the remaining master mesh
      prepareBoundingBox();
      mesh::Mesh filteredMesh("FilteredMesh", _dimensions, _mesh->isFlipNormals());
      filterMesh(filteredMesh, true);
      _mesh->clear();
      _mesh->addMesh(filteredMesh);
      _mesh->computeState();
      PRECICE_DEBUG("Master mesh after filtering, #vertices " << _mesh->vertices().size());

      if(areProvidedMeshesEmpty()) {
        std::string msg = "The re-partitioning completely filtered out the mesh " + _mesh->getName() + " received on this rank at the coupling interface. "
          "Most probably, the coupling interfaces of your coupled participants do not match geometry-wise. "
          "Please check your geometry setup again. Small overlaps or gaps are no problem. "
          "If your geometry setup is correct and if you have very different mesh resolutions on both sides, increasing the safety-factor "
          "of the decomposition strategy might be necessary.";
        PRECICE_CHECK(not _mesh->vertices().empty(), msg);
      }
    }
  } else {
    PRECICE_INFO("Broadcast mesh " << _mesh->getName());
    Event e1("partition.broadcastMesh." + _mesh->getName(), precice::syncMode);

    if (utils::MasterSlave::isSlave()) {
      com::CommunicateMesh(utils::MasterSlave::_communication).broadcastReceiveMesh(*_mesh);
    } else { // Master
      PRECICE_ASSERT(utils::MasterSlave::getRank() == 0);
      PRECICE_ASSERT(utils::MasterSlave::getSize() > 1);
      com::CommunicateMesh(utils::MasterSlave::_communication).broadcastSendMesh(*_mesh);
    }

    e1.stop();

    if (_geometricFilter == BROADCAST_FILTER) {

      PRECICE_INFO("Filter mesh " << _mesh->getName() << " by bounding-box");
      Event e2("partition.filterMeshBB." + _mesh->getName(), precice::syncMode);

      prepareBoundingBox();
      mesh::Mesh filteredMesh("FilteredMesh", _dimensions, _mesh->isFlipNormals());
      filterMesh(filteredMesh, true);

      if(areProvidedMeshesEmpty()) {
        std::string msg = "The re-partitioning completely filtered out the mesh " + _mesh->getName() + " received on this rank at the coupling interface. "
          "Most probably, the coupling interfaces of your coupled participants do not match geometry-wise. "
          "Please check your geometry setup again. Small overlaps or gaps are no problem. "
          "If your geometry setup is correct and if you have very different mesh resolutions on both sides, increasing the safety-factor "
          "of the decomposition strategy might be necessary.";
        PRECICE_CHECK(not filteredMesh.vertices().empty(), msg);
      }

      PRECICE_DEBUG("Bounding box filter, filtered from " << _mesh->vertices().size() << " vertices to " << filteredMesh.vertices().size() << " vertices.");
      _mesh->clear();
      _mesh->addMesh(filteredMesh);
      _mesh->computeState();
      e2.stop();
    } else {
      PRECICE_ASSERT(_geometricFilter == NO_FILTER);
    }
  }

  // (2) Tag vertices 1st round (i.e. who could be owned by this rank)
  PRECICE_DEBUG("Tag vertices for filtering: 1st round.");
  // go to both meshes, vertex is tagged if already one mesh tags him
  if (_fromMapping)
    _fromMapping->tagMeshFirstRound();
  if (_toMapping)
    _toMapping->tagMeshFirstRound();

  // (3) Define which vertices are owned by this rank
  PRECICE_DEBUG("Create owner information.");
  createOwnerInformation();

  // (4) Tag vertices 2nd round (what should be filtered out)
  PRECICE_DEBUG("Tag vertices for filtering: 2nd round.");
  if (_fromMapping)
    _fromMapping->tagMeshSecondRound();
  if (_toMapping)
    _toMapping->tagMeshSecondRound();

  // (5) Filter mesh according to tag
  PRECICE_INFO("Filter mesh " << _mesh->getName() << " by mappings");
  Event e5("partition.filterMeshMappings" + _mesh->getName(), precice::syncMode);
  mesh::Mesh filteredMesh("FilteredMesh", _dimensions, _mesh->isFlipNormals());
  filterMesh(filteredMesh, false);
  PRECICE_DEBUG("Mapping filter, filtered from " << _mesh->vertices().size() << " vertices to " << filteredMesh.vertices().size() << " vertices.");
  _mesh->clear();
  _mesh->addMesh(filteredMesh);
  _mesh->computeState();
  e5.stop();

  // (6) Compute distribution
  PRECICE_INFO("Feedback distribution for mesh " << _mesh->getName());
  Event e6("partition.feedbackMesh." + _mesh->getName(), precice::syncMode);
  if (utils::MasterSlave::isSlave()) {
    int numberOfVertices = _mesh->vertices().size();
    utils::MasterSlave::_communication->send(numberOfVertices, 0);
    if (numberOfVertices != 0) {
      std::vector<int> vertexIDs(numberOfVertices, -1);
      for (int i = 0; i < numberOfVertices; i++) {
        vertexIDs[i] = _mesh->vertices()[i].getGlobalIndex();
      }
      PRECICE_DEBUG("Send partition feedback to master");
      utils::MasterSlave::_communication->send(vertexIDs, 0);
    }
    int globalNumberOfVertices = -1;
    utils::MasterSlave::_communication->broadcast(globalNumberOfVertices, 0);
    PRECICE_ASSERT(globalNumberOfVertices >= 0);
    _mesh->setGlobalNumberOfVertices(globalNumberOfVertices);
  } else { // Master
    int              numberOfVertices = _mesh->vertices().size();
    std::vector<int> vertexIDs(numberOfVertices, -1);
    for (int i = 0; i < numberOfVertices; i++) {
      vertexIDs[i] = _mesh->vertices()[i].getGlobalIndex();
    }
    _mesh->getVertexDistribution()[0] = std::move(vertexIDs);

    for (int rankSlave = 1; rankSlave < utils::MasterSlave::getSize(); rankSlave++) {
      int numberOfSlaveVertices = -1;
      utils::MasterSlave::_communication->receive(numberOfSlaveVertices, rankSlave);
      PRECICE_ASSERT(numberOfSlaveVertices >= 0);
      std::vector<int> slaveVertexIDs(numberOfSlaveVertices, -1);
      if (numberOfSlaveVertices != 0) {
        PRECICE_DEBUG("Receive partition feedback from slave rank " << rankSlave);
        utils::MasterSlave::_communication->receive(slaveVertexIDs, rankSlave);
      }
      _mesh->getVertexDistribution()[rankSlave] = std::move(slaveVertexIDs);
    }
    utils::MasterSlave::_communication->broadcast(_mesh->getGlobalNumberOfVertices());
  }
  e6.stop();

  computeVertexOffsets();
}

void ReceivedPartition::filterMesh(mesh::Mesh &filteredMesh, const bool filterByBB)
{
  PRECICE_TRACE(filterByBB);

  PRECICE_DEBUG("Filter mesh " << _mesh->getName());
  PRECICE_DEBUG("Bounding mesh. #vertices: " << _mesh->vertices().size()
        << ", #edges: " << _mesh->edges().size()
        << ", #triangles: " << _mesh->triangles().size()
        << ", rank: " << utils::MasterSlave::getRank());

  std::map<int, mesh::Vertex *> vertexMap;
  std::map<int, mesh::Edge *>   edgeMap;
  int                           vertexCounter = 0;

  PRECICE_DEBUG("Filter vertices");
  for (const mesh::Vertex &vertex : _mesh->vertices()) {
    if ((filterByBB && isVertexInBB(vertex)) || (not filterByBB && vertex.isTagged())) {
      mesh::Vertex &v = filteredMesh.createVertex(vertex.getCoords());
      v.setGlobalIndex(vertex.getGlobalIndex());
      if (vertex.isTagged())
        v.tag();
      v.setOwner(vertex.isOwner());
      vertexMap[vertex.getID()] = &v;
    }
    vertexCounter++;
  }

  // Add all edges formed by the contributing vertices
  PRECICE_DEBUG("Add edges to filtered mesh");
  for (mesh::Edge &edge : _mesh->edges()) {
    int vertexIndex1 = edge.vertex(0).getID();
    int vertexIndex2 = edge.vertex(1).getID();
    if (utils::contained(vertexIndex1, vertexMap) &&
        utils::contained(vertexIndex2, vertexMap)) {
      mesh::Edge &e         = filteredMesh.createEdge(*vertexMap[vertexIndex1], *vertexMap[vertexIndex2]);
      edgeMap[edge.getID()] = &e;
    }
  }

  // Add all triangles formed by the contributing edges
  if (_dimensions == 3) {
    PRECICE_DEBUG("Add triangles to filtered mesh");
    for (mesh::Triangle &triangle : _mesh->triangles()) {
      int edgeIndex1 = triangle.edge(0).getID();
      int edgeIndex2 = triangle.edge(1).getID();
      int edgeIndex3 = triangle.edge(2).getID();
      if (utils::contained(edgeIndex1, edgeMap) &&
          utils::contained(edgeIndex2, edgeMap) &&
          utils::contained(edgeIndex3, edgeMap)) {
        filteredMesh.createTriangle(*edgeMap[edgeIndex1], *edgeMap[edgeIndex2], *edgeMap[edgeIndex3]);
      }
    }
  }

  PRECICE_DEBUG("Filtered mesh. #vertices: " << filteredMesh.vertices().size()
        << ", #edges: " << filteredMesh.edges().size()
        << ", #triangles: " << filteredMesh.triangles().size()
        << ", rank: " << utils::MasterSlave::getRank());
}

void ReceivedPartition::prepareBoundingBox()
{
  PRECICE_TRACE(_safetyFactor);

  PRECICE_DEBUG("Merge bounding boxes and increase by safety factor");

  _bb.resize(_dimensions,
             std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()));

  // Create BB around both "other" meshes
  if (_fromMapping) {
    auto other_bb = _fromMapping->getOutputMesh()->getBoundingBox();
    for (int d = 0; d < _dimensions; d++) {
      _bb[d].first = std::min(_bb[d].first, other_bb[d].first);
      _bb[d].second = std::max(_bb[d].second, other_bb[d].second);
    }
  }
  if (_toMapping) {
    auto other_bb = _toMapping->getInputMesh()->getBoundingBox();
      for (int d = 0; d < _dimensions; d++) {
        _bb[d].first = std::min(_bb[d].first, other_bb[d].first);
        _bb[d].second = std::max(_bb[d].second, other_bb[d].second);
    }
  }

  // Enlarge BB
  PRECICE_ASSERT(_safetyFactor >= 0.0);

  double maxSideLength = 1e-6; // we need some minimum > 0 here

  for (int d = 0; d < _dimensions; d++) {
    if(_bb[d].second > _bb[d].first)
      maxSideLength = std::max(maxSideLength, _bb[d].second - _bb[d].first);
  }
  for (int d = 0; d < _dimensions; d++) {
    _bb[d].second += _safetyFactor * maxSideLength;
    _bb[d].first -= _safetyFactor * maxSideLength;
    PRECICE_DEBUG("Merged BoundingBox, dim: " << d << ", first: " << _bb[d].first << ", second: " << _bb[d].second);
  }
}

bool ReceivedPartition::isVertexInBB(const mesh::Vertex &vertex)
{
  for (int d = 0; d < _dimensions; d++) {
    if (vertex.getCoords()[d] < _bb[d].first or vertex.getCoords()[d] > _bb[d].second) {
      return false;
    }
  }
  return true;
}

void ReceivedPartition::createOwnerInformation()
{
  PRECICE_TRACE();
  Event e("partition.createOwnerInformation." + _mesh->getName(), precice::syncMode);

  if (utils::MasterSlave::isSlave()) {
    int numberOfVertices = _mesh->vertices().size();
    utils::MasterSlave::_communication->send(numberOfVertices, 0);

    if (numberOfVertices != 0) {
      PRECICE_DEBUG("Tag vertices, number of vertices " << numberOfVertices);
      std::vector<int> tags(numberOfVertices, -1);
      std::vector<int> globalIDs(numberOfVertices, -1);
      bool             atInterface = false;
      for (int i = 0; i < numberOfVertices; i++) {
        globalIDs[i] = _mesh->vertices()[i].getGlobalIndex();
        if (_mesh->vertices()[i].isTagged()) {
          tags[i]     = 1;
          atInterface = true;
        } else {
          tags[i] = 0;
        }
      }
      PRECICE_DEBUG("My tags: " << tags);
      PRECICE_DEBUG("My global IDs: " << globalIDs);
      PRECICE_DEBUG("Send tags and global IDs");
      utils::MasterSlave::_communication->send(tags, 0);
      utils::MasterSlave::_communication->send(globalIDs, 0);
      utils::MasterSlave::_communication->send(atInterface, 0);

      PRECICE_DEBUG("Receive owner information");
      std::vector<int> ownerVec(numberOfVertices, -1);
      utils::MasterSlave::_communication->receive(ownerVec, 0);
      PRECICE_DEBUG("My owner information: " << ownerVec);
      setOwnerInformation(ownerVec);
    }
  }

  else if (utils::MasterSlave::isMaster()) {
    // To temporary store which vertices already have an owner
    std::vector<int> globalOwnerVec(_mesh->getGlobalNumberOfVertices(), 0);
    // The same per rank
    std::vector<std::vector<int>> slaveOwnerVecs(utils::MasterSlave::getSize());
    // Global IDs per rank
    std::vector<std::vector<int>> slaveGlobalIDs(utils::MasterSlave::getSize());
    // Tag information per rank
    std::vector<std::vector<int>> slaveTags(utils::MasterSlave::getSize());

    // Fill master data
    PRECICE_DEBUG("Tag master vertices");
    bool masterAtInterface = false;
    slaveOwnerVecs[0].resize(_mesh->vertices().size());
    slaveGlobalIDs[0].resize(_mesh->vertices().size());
    slaveTags[0].resize(_mesh->vertices().size());
    for (size_t i = 0; i < _mesh->vertices().size(); i++) {
      slaveGlobalIDs[0][i] = _mesh->vertices()[i].getGlobalIndex();
      if (_mesh->vertices()[i].isTagged()) {
        masterAtInterface = true;
        slaveTags[0][i]   = 1;
      } else {
        slaveTags[0][i] = 0;
      }
    }
    PRECICE_DEBUG("My tags: " << slaveTags[0]);

    // receive slave data
    int ranksAtInterface = 0;
    if (masterAtInterface)
      ranksAtInterface++;

    for (int rank = 1; rank < utils::MasterSlave::getSize(); rank++) {
      int localNumberOfVertices = -1;
      utils::MasterSlave::_communication->receive(localNumberOfVertices, rank);
      PRECICE_DEBUG("Rank " << rank << " has " << localNumberOfVertices << " vertices.");
      slaveOwnerVecs[rank].resize(localNumberOfVertices, 0);

      if (localNumberOfVertices != 0) {
        PRECICE_DEBUG("Receive tags from slave rank " << rank);
        utils::MasterSlave::_communication->receive(slaveTags[rank], rank);
        utils::MasterSlave::_communication->receive(slaveGlobalIDs[rank], rank);
        PRECICE_DEBUG("Rank " << rank << " has tags " << slaveTags[rank]);
        PRECICE_DEBUG("Rank " << rank << " has global IDs " << slaveGlobalIDs[rank]);
        bool atInterface = false;
        utils::MasterSlave::_communication->receive(atInterface, rank);
        if (atInterface)
          ranksAtInterface++;
      }
    }

    // Decide upon owners,
    PRECICE_DEBUG("Decide owners, first round by rough load balancing");
    int localGuess = _mesh->getGlobalNumberOfVertices() / ranksAtInterface; // Guess for a decent load balancing
    // First round: every slave gets localGuess vertices
    for (int rank = 0; rank < utils::MasterSlave::getSize(); rank++) {
      int counter = 0;
      for (size_t i = 0; i < slaveOwnerVecs[rank].size(); i++) {
         // Vertex has no owner yet and rank could be owner
        if (globalOwnerVec[slaveGlobalIDs[rank][i]] == 0 && slaveTags[rank][i] == 1) {
          slaveOwnerVecs[rank][i]                 = 1; // Now rank is owner
          globalOwnerVec[slaveGlobalIDs[rank][i]] = 1; // Vertex now has owner
          counter++;
          if (counter == localGuess)
            break;
        }
      }
    }

    // Second round: distribute all other vertices in a greedy way
    PRECICE_DEBUG("Decide owners, second round in greedy way");
    for (int rank = 0; rank < utils::MasterSlave::getSize(); rank++) {
      for (size_t i = 0; i < slaveOwnerVecs[rank].size(); i++) {
        if (globalOwnerVec[slaveGlobalIDs[rank][i]] == 0 && slaveTags[rank][i] == 1) {
          slaveOwnerVecs[rank][i]                 = 1;
          globalOwnerVec[slaveGlobalIDs[rank][i]] = rank + 1;
        }
      }
    }

    // Send information back to slaves
    for (int rank = 1; rank < utils::MasterSlave::getSize(); rank++) {
      if (not slaveTags[rank].empty())
        PRECICE_DEBUG("Send owner information to slave rank " << rank);
        utils::MasterSlave::_communication->send(slaveOwnerVecs[rank], rank);
    }
    // Master data
    PRECICE_DEBUG("My owner information: " << slaveOwnerVecs[0]);
    setOwnerInformation(slaveOwnerVecs[0]);

#ifndef NDEBUG
    for (size_t i = 0; i < globalOwnerVec.size(); i++) {
      if (globalOwnerVec[i] == 0) {
        PRECICE_DEBUG("The Vertex with global index " << i << " of mesh: " << _mesh->getName()
              << " was completely filtered out, since it has no influence on any mapping.");
      }
    }
#endif
    auto filteredVertices = std::count(globalOwnerVec.begin(), globalOwnerVec.end(), 0);
    if (filteredVertices)
      PRECICE_WARN(filteredVertices << " of " << _mesh->getGlobalNumberOfVertices()
           << " vertices of mesh " << _mesh->getName() << " have been filtered out "
           << "since they have no influence on the mapping.");

  }
}

bool ReceivedPartition::areProvidedMeshesEmpty() const {
    return (_fromMapping && not _fromMapping->getOutputMesh()->vertices().empty()) ||
           (_toMapping   && not _toMapping->getInputMesh()->vertices().empty());
}

void ReceivedPartition::setOwnerInformation(const std::vector<int> &ownerVec)
{
  size_t i = 0;
  for (mesh::Vertex &vertex : _mesh->vertices()) {
    PRECICE_ASSERT(i < ownerVec.size());
    PRECICE_ASSERT(ownerVec[i] != -1);
    vertex.setOwner(ownerVec[i] == 1);
    i++;
  }
}

} // namespace partition
} // namespace precice

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
  TRACE();
  INFO("Receive global mesh " << _mesh->getName());
  Event e("partition.receiveGlobalMesh." + _mesh->getName(), precice::syncMode);
  if (not utils::MasterSlave::_slaveMode) {
    assertion(_mesh->vertices().empty());
    com::CommunicateMesh(_m2n->getMasterCommunication()).receiveMesh(*_mesh, 0);
  }
}

void ReceivedPartition::compute()
{
  TRACE(_geometricFilter);

  // handle coupling mode first (i.e. serial participant)
  if (not utils::MasterSlave::_slaveMode && not utils::MasterSlave::_masterMode) { //coupling mode
    _mesh->setGlobalNumberOfVertices(_mesh->vertices().size());
    computeVertexOffsets();
    for (mesh::Vertex &v : _mesh->vertices()) {
      v.setOwner(true);
    }
    return;
  }

  // check to prevent false configuration
  if (not utils::MasterSlave::_slaveMode) {
    CHECK(_fromMapping || _toMapping,
          "The received mesh " << _mesh->getName()
          << " needs a mapping, either from it, to it, or both. Maybe you don't want to receive this mesh at all?")
  }


  // To understand the following steps, it is recommended to look at BU's thesis, especially Figure 69 on page 89 
  // for RBF-based filtering. https://mediatum.ub.tum.de/doc/1320661/document.pdf


  // (0) set global number of vertices before filtering
  if (utils::MasterSlave::_masterMode) {
    _mesh->setGlobalNumberOfVertices(_mesh->vertices().size());
  }

  // (1) Bounding-Box-Filter

  if (_geometricFilter == FILTER_FIRST) { //pre-filter-post-filter

    INFO("Pre-filter mesh " << _mesh->getName() << " by bounding-box");
    Event e("partition.preFilterMesh." + _mesh->getName(), precice::syncMode);

    if (utils::MasterSlave::_slaveMode) {
      prepareBoundingBox();
      com::CommunicateMesh(utils::MasterSlave::_communication).sendBoundingBox(_bb, 0);
      com::CommunicateMesh(utils::MasterSlave::_communication).receiveMesh(*_mesh, 0);

      if(areProvidedMeshesEmpty()) {
        std::string msg = "The re-partitioning completely filtered out the mesh " + _mesh->getName() +
          " received on this rank at the coupling interface. "
          "Most probably, the coupling interfaces of your coupled participants do not match geometry-wise. "
          "Please check your geometry setup again. Small overlaps or gaps are no problem. "
          "If your geometry setup is correct and if you have very different mesh resolutions on both sides, increasing the safety-factor "
          "of the decomposition strategy might be necessary.";
        CHECK(not _mesh->vertices().empty(), msg);
      }

    } else { // Master
      assertion(utils::MasterSlave::_rank == 0);
      assertion(utils::MasterSlave::_size > 1);

      for (int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++) {
        com::CommunicateMesh(utils::MasterSlave::_communication).receiveBoundingBox(_bb, rankSlave);

        DEBUG("From slave " << rankSlave << ", bounding mesh: " << _bb[0].first
              << ", " << _bb[0].second << " and " << _bb[1].first << ", " << _bb[1].second);
        mesh::Mesh slaveMesh("SlaveMesh", _dimensions, _mesh->isFlipNormals());
        filterMesh(slaveMesh, true);
        com::CommunicateMesh(utils::MasterSlave::_communication).sendMesh(slaveMesh, rankSlave);
      }

      // Now also filter the remaining master mesh
      prepareBoundingBox();
      mesh::Mesh filteredMesh("FilteredMesh", _dimensions, _mesh->isFlipNormals());
      filterMesh(filteredMesh, true);
      _mesh->clear();
      _mesh->addMesh(filteredMesh);
      _mesh->computeState();
      DEBUG("Master mesh after filtering, #vertices " << _mesh->vertices().size());

      if(areProvidedMeshesEmpty()) {
        std::string msg = "The re-partitioning completely filtered out the mesh " + _mesh->getName() + " received on this rank at the coupling interface. "
          "Most probably, the coupling interfaces of your coupled participants do not match geometry-wise. "
          "Please check your geometry setup again. Small overlaps or gaps are no problem. "
          "If your geometry setup is correct and if you have very different mesh resolutions on both sides, increasing the safety-factor "
          "of the decomposition strategy might be necessary.";
        CHECK(not _mesh->vertices().empty(), msg);
      }
    }
  } else {
    INFO("Broadcast mesh " << _mesh->getName());
    Event e1("partition.broadcastMesh." + _mesh->getName(), precice::syncMode);

    if (utils::MasterSlave::_slaveMode) {
      com::CommunicateMesh(utils::MasterSlave::_communication).broadcastReceiveMesh(*_mesh);
    } else { // Master
      assertion(utils::MasterSlave::_rank == 0);
      assertion(utils::MasterSlave::_size > 1);
      com::CommunicateMesh(utils::MasterSlave::_communication).broadcastSendMesh(*_mesh);
    }

    e1.stop();

    if (_geometricFilter == BROADCAST_FILTER) {

      INFO("Filter mesh " << _mesh->getName() << " by bounding-box");
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
        CHECK(not filteredMesh.vertices().empty(), msg);
      }

      DEBUG("Bounding box filter, filtered from " << _mesh->vertices().size() << " vertices to " << filteredMesh.vertices().size() << " vertices.");
      _mesh->clear();
      _mesh->addMesh(filteredMesh);
      _mesh->computeState();
      e2.stop();
    } else {
      assertion(_geometricFilter == NO_FILTER);
    }
  }

  // (2) Tag vertices 1st round (i.e. who could be owned by this rank)
  DEBUG("Tag vertices for filtering: 1st round.");
  // go to both meshes, vertex is tagged if already one mesh tags him
  if (_fromMapping)
    _fromMapping->tagMeshFirstRound();
  if (_toMapping)
    _toMapping->tagMeshFirstRound();

  // (3) Define which vertices are owned by this rank
  DEBUG("Create owner information.");
  createOwnerInformation();

  // (4) Tag vertices 2nd round (what should be filtered out)
  DEBUG("Tag vertices for filtering: 2nd round.");
  if (_fromMapping)
    _fromMapping->tagMeshSecondRound();
  if (_toMapping)
    _toMapping->tagMeshSecondRound();

  // (5) Filter mesh according to tag
  INFO("Filter mesh " << _mesh->getName() << " by mappings");
  Event e5("partition.filterMeshMappings" + _mesh->getName(), precice::syncMode);
  mesh::Mesh filteredMesh("FilteredMesh", _dimensions, _mesh->isFlipNormals());
  filterMesh(filteredMesh, false);
  DEBUG("Mapping filter, filtered from " << _mesh->vertices().size() << " vertices to " << filteredMesh.vertices().size() << " vertices.");
  _mesh->clear();
  _mesh->addMesh(filteredMesh);
  _mesh->computeState();
  e5.stop();

  // (6) Compute distribution
  INFO("Feedback distribution for mesh " << _mesh->getName());
  Event e6("partition.feedbackMesh." + _mesh->getName(), precice::syncMode);
  if (utils::MasterSlave::_slaveMode) {
    int numberOfVertices = _mesh->vertices().size();
    utils::MasterSlave::_communication->send(numberOfVertices, 0);
    if (numberOfVertices != 0) {
      std::vector<int> vertexIDs(numberOfVertices, -1);
      for (int i = 0; i < numberOfVertices; i++) {
        vertexIDs[i] = _mesh->vertices()[i].getGlobalIndex();
      }
      utils::MasterSlave::_communication->send(vertexIDs, 0);
    }
    int globalNumberOfVertices = -1;
    utils::MasterSlave::_communication->broadcast(globalNumberOfVertices, 0);
    assertion(globalNumberOfVertices != -1);
    _mesh->setGlobalNumberOfVertices(globalNumberOfVertices);
  } else { // Master
    int              numberOfVertices = _mesh->vertices().size();
    std::vector<int> vertexIDs(numberOfVertices, -1);
    for (int i = 0; i < numberOfVertices; i++) {
      vertexIDs[i] = _mesh->vertices()[i].getGlobalIndex();
    }
    _mesh->getVertexDistribution()[0] = vertexIDs;

    for (int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++) {
      int numberOfSlaveVertices = -1;
      utils::MasterSlave::_communication->receive(numberOfSlaveVertices, rankSlave);
      std::vector<int> slaveVertexIDs;
      if (numberOfSlaveVertices != 0) {
        utils::MasterSlave::_communication->receive(slaveVertexIDs, rankSlave);
      }
      _mesh->getVertexDistribution()[rankSlave] = slaveVertexIDs;
    }
    utils::MasterSlave::_communication->broadcast(_mesh->getGlobalNumberOfVertices());
  }
  e6.stop();

  computeVertexOffsets();
}

void ReceivedPartition::filterMesh(mesh::Mesh &filteredMesh, const bool filterByBB)
{
  TRACE(filterByBB);

  DEBUG("Bounding mesh. #vertices: " << _mesh->vertices().size()
        << ", #edges: " << _mesh->edges().size()
        << ", #triangles: " << _mesh->triangles().size()
        << ", rank: " << utils::MasterSlave::_rank);

  std::map<int, mesh::Vertex *> vertexMap;
  std::map<int, mesh::Edge *>   edgeMap;
  int                           vertexCounter = 0;

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

  DEBUG("Filtered mesh. #vertices: " << filteredMesh.vertices().size()
        << ", #edges: " << filteredMesh.edges().size()
        << ", #triangles: " << filteredMesh.triangles().size()
        << ", rank: " << utils::MasterSlave::_rank);
}

void ReceivedPartition::prepareBoundingBox()
{
  TRACE(_safetyFactor);

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
  assertion(_safetyFactor >= 0.0);

  double maxSideLength = 1e-6; // we need some minimum > 0 here

  for (int d = 0; d < _dimensions; d++) {
    maxSideLength = std::max(maxSideLength, _bb[d].second - _bb[d].first);
  }
  for (int d = 0; d < _dimensions; d++) {
    _bb[d].second += _safetyFactor * maxSideLength;
    _bb[d].first -= _safetyFactor * maxSideLength;
    DEBUG("Merged BoundingBox, dim: " << d << ", first: " << _bb[d].first << ", second: " << _bb[d].second);
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
  TRACE();

  if (utils::MasterSlave::_slaveMode) {
    int numberOfVertices = _mesh->vertices().size();
    utils::MasterSlave::_communication->send(numberOfVertices, 0);

    if (numberOfVertices != 0) {
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
      DEBUG("My tags: " << tags);
      DEBUG("My global IDs: " << globalIDs);
      utils::MasterSlave::_communication->send(tags, 0);
      utils::MasterSlave::_communication->send(globalIDs, 0);
      utils::MasterSlave::_communication->send(atInterface, 0);

      std::vector<int> ownerVec(numberOfVertices, -1);
      utils::MasterSlave::_communication->receive(ownerVec, 0);
      DEBUG("My owner information: " << ownerVec);
      setOwnerInformation(ownerVec);
    }
  }

  else if (utils::MasterSlave::_masterMode) {
    // To temporary store which vertices already have an owner
    std::vector<int> globalOwnerVec(_mesh->getGlobalNumberOfVertices(), 0);
    // The same per rank
    std::vector<std::vector<int>> slaveOwnerVecs(utils::MasterSlave::_size);
    // Global IDs per rank
    std::vector<std::vector<int>> slaveGlobalIDs(utils::MasterSlave::_size);
    // Tag information per rank
    std::vector<std::vector<int>> slaveTags(utils::MasterSlave::_size);

    // Fill master data
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
    DEBUG("My tags: " << slaveTags[0]);

    // receive slave data
    int ranksAtInterface = 0;
    if (masterAtInterface)
      ranksAtInterface++;

    for (int rank = 1; rank < utils::MasterSlave::_size; rank++) {
      int localNumberOfVertices = -1;
      utils::MasterSlave::_communication->receive(localNumberOfVertices, rank);
      DEBUG("Rank " << rank << " has " << localNumberOfVertices << " vertices.");
      slaveOwnerVecs[rank].resize(localNumberOfVertices, 0);

      if (localNumberOfVertices != 0) {
        utils::MasterSlave::_communication->receive(slaveTags[rank], rank);
        utils::MasterSlave::_communication->receive(slaveGlobalIDs[rank], rank);
        DEBUG("Rank " << rank << " has this tags " << slaveTags[rank]);
        DEBUG("Rank " << rank << " has this global IDs " << slaveGlobalIDs[rank]);
        bool atInterface = false;
        utils::MasterSlave::_communication->receive(atInterface, rank);
        if (atInterface)
          ranksAtInterface++;
      }
    }

    // Decide upon owners,
    int localGuess = _mesh->getGlobalNumberOfVertices() / ranksAtInterface; // Guess for a decent load balancing
    // First round: every slave gets localGuess vertices
    for (int rank = 0; rank < utils::MasterSlave::_size; rank++) {
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
    for (int rank = 0; rank < utils::MasterSlave::_size; rank++) {
      for (size_t i = 0; i < slaveOwnerVecs[rank].size(); i++) {
        if (globalOwnerVec[slaveGlobalIDs[rank][i]] == 0 && slaveTags[rank][i] == 1) {
          slaveOwnerVecs[rank][i]                 = 1;
          globalOwnerVec[slaveGlobalIDs[rank][i]] = rank + 1;
        }
      }
    }

    // Send information back to slaves
    for (int rank = 1; rank < utils::MasterSlave::_size; rank++) {
      if (not slaveTags[rank].empty())
        utils::MasterSlave::_communication->send(slaveOwnerVecs[rank], rank);
    }
    // Master data
    DEBUG("My owner information: " << slaveOwnerVecs[0]);
    setOwnerInformation(slaveOwnerVecs[0]);

#ifndef NDEBUG
    for (size_t i = 0; i < globalOwnerVec.size(); i++) {
      if (globalOwnerVec[i] == 0) {
        WARN("The Vertex with global index " << i << " of mesh: " << _mesh->getName()
             << " was completely filtered out, since it has no influence on any mapping.");
      }
    }
#endif
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
    assertion(i < ownerVec.size());
    assertion(ownerVec[i] != -1);
    vertex.setOwner(ownerVec[i] == 1);
    i++;
  }
}

} // namespace partition
} // namespace precice

#include "partition/ReceivedPartition.hpp"
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
#include "mapping/Mapping.hpp"
#include "mapping/SharedPointer.hpp"
#include "mesh/BoundingBox.hpp"
#include "mesh/Filter.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "partition/Partition.hpp"
#include "utils/Event.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/assertion.hpp"
#include "utils/fmt.hpp"

using precice::utils::Event;

namespace precice {
extern bool syncMode;

namespace partition {

ReceivedPartition::ReceivedPartition(
    const mesh::PtrMesh &mesh, GeometricFilter geometricFilter, double safetyFactor, bool allowDirectAccess)
    : Partition(mesh),
      _geometricFilter(geometricFilter),
      _bb(mesh->getDimensions()),
      _dimensions(mesh->getDimensions()),
      _safetyFactor(safetyFactor),
      _allowDirectAccess(allowDirectAccess)
{
}

void ReceivedPartition::communicate()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_mesh->vertices().empty());

  // for two-level initialization, receive mesh partitions
  if (m2n().usesTwoLevelInitialization()) {
    PRECICE_INFO("Receive mesh partitions for mesh {}", _mesh->getName());
    Event e("partition.receiveMeshPartitions." + _mesh->getName(), precice::syncMode);

    if (utils::MasterSlave::isMaster()) {
      // Master receives remote mesh's global number of vertices
      int globalNumberOfVertices = -1;
      m2n().getMasterCommunication()->receive(globalNumberOfVertices, 0);
      _mesh->setGlobalNumberOfVertices(globalNumberOfVertices);
    }

    // each rank receives max/min global vertex indices from connected remote ranks
    m2n().broadcastReceiveAll(_remoteMinGlobalVertexIDs, *_mesh);
    m2n().broadcastReceiveAll(_remoteMaxGlobalVertexIDs, *_mesh);
    // each rank receives mesh partition from connected remote ranks
    m2n().broadcastReceiveAllMesh(*_mesh);

  } else {
    // for one-level initialization receive complete mesh on master
    PRECICE_INFO("Receive global mesh {}", _mesh->getName());
    Event e("partition.receiveGlobalMesh." + _mesh->getName(), precice::syncMode);

    if (not utils::MasterSlave::isSlave()) {
      // a ReceivedPartition can only have one communication, @todo nicer design
      com::CommunicateMesh(m2n().getMasterCommunication()).receiveMesh(*_mesh, 0);
      _mesh->setGlobalNumberOfVertices(_mesh->vertices().size());
    }
  }

  // for both initialization concepts broadcast and set the global number of vertices
  if (utils::MasterSlave::isMaster()) {
    utils::MasterSlave::getCommunication()->broadcast(_mesh->getGlobalNumberOfVertices());
  }
  if (utils::MasterSlave::isSlave()) {
    int globalNumberOfVertices = -1;
    utils::MasterSlave::getCommunication()->broadcast(globalNumberOfVertices, 0);
    PRECICE_ASSERT(globalNumberOfVertices >= 0);
    _mesh->setGlobalNumberOfVertices(globalNumberOfVertices);
  }
}

void ReceivedPartition::compute()
{
  PRECICE_TRACE();

  // handle coupling mode first (i.e. serial participant)
  if (!utils::MasterSlave::isParallel()) { //coupling mode
    PRECICE_DEBUG("Handle partition data structures for serial participant");

    if (_allowDirectAccess) {
      // Prepare the bounding boxes
      prepareBoundingBox();
      // Filter out vertices not laying in the bounding box
      mesh::Mesh filteredMesh("FilteredMesh", _dimensions, mesh::Mesh::MESH_ID_UNDEFINED);
      // To discuss: maybe check this somewhere in the SolverInterfaceImpl, as we have now a similar check for the parallel case
      PRECICE_CHECK(!_bb.empty(), "You are running this participant in serial mode and the bounding box on mesh \"{}\", is empty. Did you call setMeshAccessRegion with valid data?", _mesh->getName());
      unsigned int nFilteredVertices = 0;
      mesh::filterMesh(filteredMesh, *_mesh, [&](const mesh::Vertex &v) { if(!_bb.contains(v))
              ++nFilteredVertices;
          return _bb.contains(v); });

      if (nFilteredVertices > 0) {
        PRECICE_WARN("{} vertices on mesh \"{}\" have been filtered out due to the defined bounding box in \"setMeshAccessRegion\" "
                     "in serial mode. Associated data values of the filtered vertices will be filled with zero values in order to provide valid data for other participants when reading data.",
                     nFilteredVertices, _mesh->getName());
      }

      _mesh->clear();
      _mesh->addMesh(filteredMesh);
    }

    int vertexCounter = 0;
    for (mesh::Vertex &v : _mesh->vertices()) {
      v.setOwner(true);
      _mesh->getVertexDistribution()[0].push_back(vertexCounter);
      vertexCounter++;
    }
    _mesh->getVertexOffsets().push_back(vertexCounter);
    return;
  }

  // check to prevent false configuration
  if (not utils::MasterSlave::isSlave()) {
    PRECICE_CHECK(hasAnyMapping() || _allowDirectAccess,
                  "The received mesh {} needs a mapping, either from it, to it, or both. Maybe you don't want to receive this mesh at all?",
                  _mesh->getName());
  }

  // To better understand steps (2) to (5), it is recommended to look at BU's thesis, especially Figure 69 on page 89
  // for RBF-based filtering. https://mediatum.ub.tum.de/doc/1320661/document.pdf

  // (1) Bounding-Box-Filter
  filterByBoundingBox();

  // (2) Tag vertices 1st round (i.e. who could be owned by this rank)
  PRECICE_DEBUG("Tag vertices for filtering: 1st round.");
  // go to both meshes, vertex is tagged if already one mesh tags him
  tagMeshFirstRound();

  // (3) Define which vertices are owned by this rank
  PRECICE_DEBUG("Create owner information.");
  createOwnerInformation();

  // (4) Tag vertices 2nd round (what should be filtered out)
  PRECICE_DEBUG("Tag vertices for filtering: 2nd round.");
  tagMeshSecondRound();

  // (5) Filter mesh according to tag
  PRECICE_INFO("Filter mesh {} by mappings", _mesh->getName());
  Event      e5("partition.filterMeshMappings" + _mesh->getName(), precice::syncMode);
  mesh::Mesh filteredMesh("FilteredMesh", _dimensions, mesh::Mesh::MESH_ID_UNDEFINED);
  mesh::filterMesh(filteredMesh, *_mesh, [&](const mesh::Vertex &v) { return v.isTagged(); });
  PRECICE_DEBUG("Mapping filter, filtered from {} to {} vertices, {} to {} edges, and {} to {} triangles.",
                _mesh->vertices().size(), filteredMesh.vertices().size(),
                _mesh->edges().size(), filteredMesh.edges().size(),
                _mesh->triangles().size(), filteredMesh.triangles().size());

  _mesh->clear();
  _mesh->addMesh(filteredMesh);
  e5.stop();

  // (6) Compute vertex distribution or local communication map
  if (m2n().usesTwoLevelInitialization()) {

    PRECICE_INFO("Compute communication map for mesh {}", _mesh->getName());
    Event e6("partition.computeCommunicationMap." + _mesh->getName(), precice::syncMode);

    // Fill two data structures: remoteCommunicationMap and this rank's communication map (_mesh->getCommunicationMap()).
    // remoteCommunicationMap: connectedRank -> {remote local vertex index}
    // _mesh->getCommunicationMap(): connectedRank -> {this rank's local vertex index}
    // A vertex belongs to a specific connected rank if its global vertex ID lies within the ranks min and max.
    std::map<int, std::vector<int>> remoteCommunicationMap;

    for (size_t vertexIndex = 0; vertexIndex < _mesh->vertices().size(); ++vertexIndex) {
      for (size_t rankIndex = 0; rankIndex < _mesh->getConnectedRanks().size(); ++rankIndex) {
        int globalVertexIndex = _mesh->vertices()[vertexIndex].getGlobalIndex();
        if (globalVertexIndex <= _remoteMaxGlobalVertexIDs[rankIndex] && globalVertexIndex >= _remoteMinGlobalVertexIDs[rankIndex]) {
          int remoteRank = _mesh->getConnectedRanks()[rankIndex];
          remoteCommunicationMap[remoteRank].push_back(globalVertexIndex - _remoteMinGlobalVertexIDs[rankIndex]); //remote local vertex index
          _mesh->getCommunicationMap()[remoteRank].push_back(vertexIndex);                                        //this rank's local vertex index
        }
      }
    }

    // communicate remote communication map to all remote connected ranks
    m2n().scatterAllCommunicationMap(remoteCommunicationMap, *_mesh);

  } else {

    PRECICE_INFO("Feedback distribution for mesh {}", _mesh->getName());
    Event e6("partition.feedbackMesh." + _mesh->getName(), precice::syncMode);
    if (utils::MasterSlave::isSlave()) {
      int numberOfVertices = _mesh->vertices().size();
      utils::MasterSlave::getCommunication()->send(numberOfVertices, 0);
      if (numberOfVertices != 0) {
        std::vector<int> vertexIDs(numberOfVertices, -1);
        for (int i = 0; i < numberOfVertices; i++) {
          vertexIDs[i] = _mesh->vertices()[i].getGlobalIndex();
        }
        PRECICE_DEBUG("Send partition feedback to master");
        utils::MasterSlave::getCommunication()->send(vertexIDs, 0);
      }
    } else { // Master
      int              numberOfVertices = _mesh->vertices().size();
      std::vector<int> vertexIDs(numberOfVertices, -1);
      for (int i = 0; i < numberOfVertices; i++) {
        vertexIDs[i] = _mesh->vertices()[i].getGlobalIndex();
      }
      _mesh->getVertexDistribution()[0] = std::move(vertexIDs);

      for (int rankSlave : utils::MasterSlave::allSlaves()) {
        int numberOfSlaveVertices = -1;
        utils::MasterSlave::getCommunication()->receive(numberOfSlaveVertices, rankSlave);
        PRECICE_ASSERT(numberOfSlaveVertices >= 0);
        std::vector<int> slaveVertexIDs(numberOfSlaveVertices, -1);
        if (numberOfSlaveVertices != 0) {
          PRECICE_DEBUG("Receive partition feedback from slave rank {}", rankSlave);
          utils::MasterSlave::getCommunication()->receive(slaveVertexIDs, rankSlave);
        }
        _mesh->getVertexDistribution()[rankSlave] = std::move(slaveVertexIDs);
      }
    }
  }

  // (7) Compute vertex offsets
  PRECICE_DEBUG("Compute vertex offsets");
  if (utils::MasterSlave::isSlave()) {

    // send number of vertices
    PRECICE_DEBUG("Send number of vertices: {}", _mesh->vertices().size());
    int numberOfVertices = _mesh->vertices().size();
    utils::MasterSlave::getCommunication()->send(numberOfVertices, 0);

    // set vertex offsets
    utils::MasterSlave::getCommunication()->broadcast(_mesh->getVertexOffsets(), 0);
    PRECICE_DEBUG("My vertex offsets: {}", _mesh->getVertexOffsets());

  } else if (utils::MasterSlave::isMaster()) {

    _mesh->getVertexOffsets().resize(utils::MasterSlave::getSize());
    _mesh->getVertexOffsets()[0] = _mesh->vertices().size();

    // receive number of slave vertices and fill vertex offsets
    for (int rankSlave : utils::MasterSlave::allSlaves()) {
      int numberOfSlaveVertices = -1;
      utils::MasterSlave::getCommunication()->receive(numberOfSlaveVertices, rankSlave);
      _mesh->getVertexOffsets()[rankSlave] = numberOfSlaveVertices + _mesh->getVertexOffsets()[rankSlave - 1];
    }

    // broadcast vertex offsets
    PRECICE_DEBUG("My vertex offsets: {}", _mesh->getVertexOffsets());
    utils::MasterSlave::getCommunication()->broadcast(_mesh->getVertexOffsets());
  }
}

namespace {
auto errorMeshFilteredOut(const std::string &meshName, const int rank)
{
  return fmt::format("The re-partitioning completely filtered out the mesh \"{0}\" received on rank {1} "
                     "at the coupling interface, although the provided mesh partition on this rank is "
                     "non-empty. Most probably, the coupling interfaces of your coupled participants do "
                     "not match geometry-wise. Please check your geometry setup again. Small overlaps or "
                     "gaps are no problem. If your geometry setup is correct and if you have very different "
                     "mesh resolutions on both sides, you may want to increase the safety-factor: "
                     "\"<use-mesh mesh=\"{0} \" ... safety-factor=\"N\"/> (default value is 0.5) of the "
                     "decomposition strategy or disable the filtering completely: "
                     "\"<use-mesh mesh=\"{0}\" ... geometric-filter=\"no-filter\" />",
                     meshName, rank);
}
} // namespace

void ReceivedPartition::filterByBoundingBox()
{
  PRECICE_TRACE(_geometricFilter);

  if (m2n().usesTwoLevelInitialization()) {
    std::string msg = "The received mesh " + _mesh->getName() +
                      " cannot solely be filtered on the master rank "
                      "(option \"filter-on-master\") if it is communicated by an m2n communication that uses "
                      "two-level initialization. Use \"filter-on-slaves\" or \"no-filter\" instead.";
    PRECICE_CHECK(_geometricFilter != ON_MASTER, msg);
  }

  prepareBoundingBox();

  if (_geometricFilter == ON_MASTER) { //filter on master and communicate reduced mesh then

    PRECICE_ASSERT(not m2n().usesTwoLevelInitialization());
    PRECICE_INFO("Pre-filter mesh {} by bounding box on master", _mesh->getName());
    Event e("partition.preFilterMesh." + _mesh->getName(), precice::syncMode);

    if (utils::MasterSlave::isSlave()) {
      PRECICE_DEBUG("Send bounding box to master");
      com::CommunicateBoundingBox(utils::MasterSlave::getCommunication()).sendBoundingBox(_bb, 0);
      PRECICE_DEBUG("Receive filtered mesh");
      com::CommunicateMesh(utils::MasterSlave::getCommunication()).receiveMesh(*_mesh, 0);

      if (isAnyProvidedMeshNonEmpty()) {
        PRECICE_CHECK(not _mesh->vertices().empty(), errorMeshFilteredOut(_mesh->getName(), utils::MasterSlave::getRank()));
      }

    } else { // Master
      PRECICE_ASSERT(utils::MasterSlave::getRank() == 0);
      PRECICE_ASSERT(utils::MasterSlave::getSize() > 1);

      for (int rankSlave : utils::MasterSlave::allSlaves()) {
        mesh::BoundingBox slaveBB(_bb.getDimension());
        com::CommunicateBoundingBox(utils::MasterSlave::getCommunication()).receiveBoundingBox(slaveBB, rankSlave);

        PRECICE_DEBUG("From slave {}, bounding mesh: {}", rankSlave, slaveBB);
        mesh::Mesh slaveMesh("SlaveMesh", _dimensions, mesh::Mesh::MESH_ID_UNDEFINED);
        mesh::filterMesh(slaveMesh, *_mesh, [&slaveBB](const mesh::Vertex &v) { return slaveBB.contains(v); });
        PRECICE_DEBUG("Send filtered mesh to slave: {}", rankSlave);
        com::CommunicateMesh(utils::MasterSlave::getCommunication()).sendMesh(slaveMesh, rankSlave);
      }

      // Now also filter the remaining master mesh
      mesh::Mesh filteredMesh("FilteredMesh", _dimensions, mesh::Mesh::MESH_ID_UNDEFINED);
      mesh::filterMesh(filteredMesh, *_mesh, [&](const mesh::Vertex &v) { return _bb.contains(v); });
      PRECICE_DEBUG("Master mesh, filtered from {} to {} vertices, {} to {} edges, and {} to {} triangles.",
                    _mesh->vertices().size(), filteredMesh.vertices().size(),
                    _mesh->edges().size(), filteredMesh.edges().size(),
                    _mesh->triangles().size(), filteredMesh.triangles().size());
      _mesh->clear();
      _mesh->addMesh(filteredMesh);

      if (isAnyProvidedMeshNonEmpty()) {
        PRECICE_CHECK(not _mesh->vertices().empty(), errorMeshFilteredOut(_mesh->getName(), utils::MasterSlave::getRank()));
      }
    }
  } else {
    if (not m2n().usesTwoLevelInitialization()) {
      PRECICE_INFO("Broadcast mesh {}", _mesh->getName());
      Event e("partition.broadcastMesh." + _mesh->getName(), precice::syncMode);

      if (utils::MasterSlave::isSlave()) {
        com::CommunicateMesh(utils::MasterSlave::getCommunication()).broadcastReceiveMesh(*_mesh);
      } else { // Master
        PRECICE_ASSERT(utils::MasterSlave::isMaster());
        com::CommunicateMesh(utils::MasterSlave::getCommunication()).broadcastSendMesh(*_mesh);
      }
    }
    if (_geometricFilter == ON_SLAVES) {

      PRECICE_INFO("Filter mesh {} by bounding box on slaves", _mesh->getName());
      Event e("partition.filterMeshBB." + _mesh->getName(), precice::syncMode);

      mesh::Mesh filteredMesh("FilteredMesh", _dimensions, mesh::Mesh::MESH_ID_UNDEFINED);
      mesh::filterMesh(filteredMesh, *_mesh, [&](const mesh::Vertex &v) { return _bb.contains(v); });

      PRECICE_DEBUG("Bounding box filter, filtered from {} to {} vertices, {} to {} edges, and {} to {} triangles.",
                    _mesh->vertices().size(), filteredMesh.vertices().size(),
                    _mesh->edges().size(), filteredMesh.edges().size(),
                    _mesh->triangles().size(), filteredMesh.triangles().size());

      _mesh->clear();
      _mesh->addMesh(filteredMesh);
      if (isAnyProvidedMeshNonEmpty()) {
        PRECICE_CHECK(not _mesh->vertices().empty(), errorMeshFilteredOut(_mesh->getName(), utils::MasterSlave::getRank()));
      }
    } else {
      PRECICE_ASSERT(_geometricFilter == NO_FILTER);
    }
  }
}

void ReceivedPartition::compareBoundingBoxes()
{
  PRECICE_TRACE();

  _mesh->clear();
  _mesh->clearPartitioning();
  _boundingBoxPrepared = false;
  _remoteMinGlobalVertexIDs.clear();
  _remoteMaxGlobalVertexIDs.clear();

  // @todo handle coupling mode (i.e. serial participant)
  // @todo treatment of multiple m2ns
  PRECICE_ASSERT(_m2ns.size() == 1);
  if (not m2n().usesTwoLevelInitialization())
    return;

  // receive and broadcast number of remote ranks
  int numberOfRemoteRanks = -1;
  if (utils::MasterSlave::isMaster()) {
    m2n().getMasterCommunication()->receive(numberOfRemoteRanks, 0);
    utils::MasterSlave::getCommunication()->broadcast(numberOfRemoteRanks);
  } else {
    PRECICE_ASSERT(utils::MasterSlave::isSlave());
    utils::MasterSlave::getCommunication()->broadcast(numberOfRemoteRanks, 0);
  }

  // define and initialize remote bounding box map
  mesh::Mesh::BoundingBoxMap remoteBBMap;
  mesh::BoundingBox          initialBB(_mesh->getDimensions());

  for (int remoteRank = 0; remoteRank < numberOfRemoteRanks; remoteRank++) {
    remoteBBMap.emplace(remoteRank, initialBB);
  }

  // receive and broadcast remote bounding box map
  if (utils::MasterSlave::isMaster()) {
    com::CommunicateBoundingBox(m2n().getMasterCommunication()).receiveBoundingBoxMap(remoteBBMap, 0);
    com::CommunicateBoundingBox(utils::MasterSlave::getCommunication()).broadcastSendBoundingBoxMap(remoteBBMap);
  } else {
    PRECICE_ASSERT(utils::MasterSlave::isSlave());
    com::CommunicateBoundingBox(utils::MasterSlave::getCommunication()).broadcastReceiveBoundingBoxMap(remoteBBMap);
  }

  // prepare local bounding box
  prepareBoundingBox();

  if (utils::MasterSlave::isMaster()) {                 // Master
    std::map<int, std::vector<int>> connectionMap;      //local ranks -> {remote ranks}
    std::vector<int>                connectedRanksList; // local ranks with any connection

    // connected ranks for master
    _mesh->getConnectedRanks().clear();
    for (auto &remoteBB : remoteBBMap) {
      if (_bb.overlapping(remoteBB.second)) {
        _mesh->getConnectedRanks().push_back(remoteBB.first); //connected remote ranks for this rank
      }
    }
    if (not _mesh->getConnectedRanks().empty()) {
      connectionMap[0] = _mesh->getConnectedRanks();
      connectedRanksList.push_back(0);
    }

    // receive connected ranks from slaves and add them to the connection map
    for (int rank : utils::MasterSlave::allSlaves()) {
      std::vector<int> slaveConnectedRanks;
      int              connectedRanksSize = -1;
      utils::MasterSlave::getCommunication()->receive(connectedRanksSize, rank);
      if (connectedRanksSize != 0) {
        connectedRanksList.push_back(rank);
        utils::MasterSlave::getCommunication()->receive(slaveConnectedRanks, rank);
        connectionMap[rank] = slaveConnectedRanks;
      }
    }

    // send connectionMap to other master
    m2n().getMasterCommunication()->send(connectedRanksList, 0);
    PRECICE_CHECK(not connectionMap.empty(),
                  "The mesh \"{}\" of this participant seems to have no partitions at the coupling interface. "
                  "Check that both mapped meshes are describing the same geometry. "
                  "If you deal with very different mesh resolutions, consider increasing the safety-factor in the <use-mesh /> tag.",
                  _mesh->getName());
    com::CommunicateBoundingBox(m2n().getMasterCommunication()).sendConnectionMap(connectionMap, 0);
  } else {
    PRECICE_ASSERT(utils::MasterSlave::isSlave());

    _mesh->getConnectedRanks().clear();
    for (const auto &remoteBB : remoteBBMap) {
      if (_bb.overlapping(remoteBB.second)) {
        _mesh->getConnectedRanks().push_back(remoteBB.first);
      }
    }

    // send connected ranks to master
    utils::MasterSlave::getCommunication()->send(static_cast<int>(_mesh->getConnectedRanks().size()), 0);
    if (not _mesh->getConnectedRanks().empty()) {
      utils::MasterSlave::getCommunication()->send(_mesh->getConnectedRanks(), 0);
    }
  }
}

void ReceivedPartition::prepareBoundingBox()
{
  PRECICE_TRACE(_safetyFactor);

  if (_boundingBoxPrepared)
    return;

  PRECICE_DEBUG("Merge bounding boxes and increase by safety factor");

  // Reset the BoundingBox
  _bb = mesh::BoundingBox{_dimensions};

  // Create BB around all "other" meshes
  for (mapping::PtrMapping &fromMapping : _fromMappings) {
    auto other_bb = fromMapping->getOutputMesh()->getBoundingBox();
    _bb.expandBy(other_bb);
    _bb.scaleBy(_safetyFactor);
    _boundingBoxPrepared = true;
  }
  for (mapping::PtrMapping &toMapping : _toMappings) {
    auto other_bb = toMapping->getInputMesh()->getBoundingBox();
    _bb.expandBy(other_bb);
    _bb.scaleBy(_safetyFactor);
    _boundingBoxPrepared = true;
  }

  // Expand by user-defined bounding box in case a direct access is desired
  if (_allowDirectAccess) {
    auto &other_bb = _mesh->getBoundingBox();
    _bb.expandBy(other_bb);

    // The safety factor is for mapping based partitionings applied, as usual.
    // For the direct access, however, we don't apply any safety factor scaling.
    // If the user defines a safety factor and the partitioning is entirely based
    // on the defined access region (setMeshAccessRegion), we raise a warning
    // to inform the user
    const float defaultSafetyFactor = 0.5;
    if (utils::MasterSlave::isMaster() && !hasAnyMapping() && (_safetyFactor != defaultSafetyFactor)) {
      PRECICE_WARN("The received mesh \"{}\" was entirely partitioned based on the defined access region "
                   "(setMeshAccessRegion) and a safety-factor was defined. However, the safety factor "
                   "will be ignored in this case. You may want to modify the access region by modifying "
                   "the specified region in the function itself.",
                   _mesh->getName());
    }
    _boundingBoxPrepared = true;
  }
}

void ReceivedPartition::createOwnerInformation()
{
  PRECICE_TRACE();
  Event e("partition.createOwnerInformation." + _mesh->getName(), precice::syncMode);

  /*
    We follow different approaches for two-level and one-level methods. For 1LI, a centeralized
    approach is followed, while the 2LI employs a local/parallel scheme for assigning vertices
    to corresponding ranks.
  */

  if (m2n().usesTwoLevelInitialization()) {
    /*
    This function ensures that each vertex is owned by only a single rank and
    is not shared among ranks. Initially, the vertices are checked against the 
    bounding box of each rank. If a vertex fits into only a single bounding box, 
    the vertex is assigned to that rank. If it fits to various bbs, the rank with the 
    lowest number of vertices gets ownership to keep the load as balanced as 
    possible. 
    
    Following steps are taken:

    1- receive local bb map from master
    2- filter bb map to keep the connected ranks
    3- own the vertices that only fit into this rank's bb
    4- send number of owned vertices and the list of shared vertices to neighbors
    5- for the remaining vertices: check if we have less vertices -> own it!
    */

    // #1: receive local bb map from master
    // Define and initialize localBBMap to save local bbs

    mesh::Mesh::BoundingBoxMap localBBMap;
    for (Rank rank = 0; rank < utils::MasterSlave::getSize(); rank++) {
      localBBMap.emplace(rank, mesh::BoundingBox(_dimensions));
    }

    // Define a bb map to save the local connected ranks and respective boundingboxes
    mesh::Mesh::BoundingBoxMap localConnectedBBMap;

    // store global IDs and list of possible shared vertices

    // Global IDs to be saved in:
    std::vector<int> sharedVerticesGlobalIDs;
    // Local IDs to be saved in:
    std::vector<int> sharedVerticesLocalIDs;

    // store possible shared vertices in a map to communicate with neighbors, map: rank -> vertex_global_id
    std::map<int, std::vector<int>> sharedVerticesSendMap;

    // receive list of possible shared vertices from neighboring ranks
    std::map<int, std::vector<int>> sharedVerticesReceiveMap;

    if (utils::MasterSlave::isMaster()) {

      // Insert masters bounding box
      localBBMap.at(0) = _bb;

      // master receives local bb from each slave rank
      for (int rankSlave = 1; rankSlave < utils::MasterSlave::getSize(); rankSlave++) {
        com::CommunicateBoundingBox(utils::MasterSlave::getCommunication()).receiveBoundingBox(localBBMap.at(rankSlave), rankSlave);
      }

      // master broadcast localBBMap to all slaves
      com::CommunicateBoundingBox(utils::MasterSlave::getCommunication()).broadcastSendBoundingBoxMap(localBBMap);
    } else if (utils::MasterSlave::isSlave()) {
      // slaves send local bb to master
      com::CommunicateBoundingBox(utils::MasterSlave::getCommunication()).sendBoundingBox(_bb, 0);
      // slaves receive localBBMap from master
      com::CommunicateBoundingBox(utils::MasterSlave::getCommunication()).broadcastReceiveBoundingBoxMap(localBBMap);
    }

    // #2: filter bb map to keep the connected ranks
    // remove the own bb from the map since we compare the own bb only with other ranks bb.
    localBBMap.erase(utils::MasterSlave::getRank());
    // find and store local connected ranks
    for (const auto &localBB : localBBMap) {
      if (_bb.overlapping(localBB.second)) {
        localConnectedBBMap.emplace(localBB.first, localBB.second);
      }
    }

    // #3: check vertices and keep only those that fit into the current rank's bb
    const int numberOfVertices = _mesh->vertices().size();
    PRECICE_DEBUG("Tag vertices, number of vertices {}", numberOfVertices);
    std::vector<int> tags(numberOfVertices, -1);
    std::vector<int> globalIDs(numberOfVertices, -1);
    bool             atInterface        = false;
    int              ownedVerticesCount = 0; // number of vertices owned by this rank
    for (int i = 0; i < numberOfVertices; i++) {
      globalIDs[i] = _mesh->vertices()[i].getGlobalIndex();
      if (_mesh->vertices()[i].isTagged()) {
        bool vertexIsShared = false;
        for (const auto &neighborRank : localConnectedBBMap) {
          if (neighborRank.second.contains(_mesh->vertices()[i])) {
            vertexIsShared = true;
            sharedVerticesSendMap[neighborRank.first].push_back(globalIDs[i]);
            sharedVerticesGlobalIDs.push_back(globalIDs[i]);
            sharedVerticesLocalIDs.push_back(i);
          }
        }

        if (not vertexIsShared) {
          tags[i]     = 1;
          atInterface = true;
          ownedVerticesCount++;
        }
      }

      else {
        tags[i] = 0;
      }
    }

    // #4: Exchange number of already owned vertices with the neighbors

    // to store receive requests.
    std::vector<com::PtrRequest> vertexNumberRequests;

    // Define and initialize for load balancing
    std::map<int, int> neighborRanksVertexCount;
    for (auto &neighborRank : localConnectedBBMap) {
      neighborRanksVertexCount.emplace(neighborRank.first, 0);
    }

    // Asynchronous recieve number of owned vertices from neighbor ranks
    for (auto &neighborRank : localConnectedBBMap) {
      auto request = utils::MasterSlave::getCommunication()->aReceive(neighborRanksVertexCount.at(neighborRank.first), neighborRank.first);
      vertexNumberRequests.push_back(request);
    }

    // Synchronous send number of owned vertices to neighbor ranks
    for (auto &neighborRank : localConnectedBBMap) {
      utils::MasterSlave::getCommunication()->send(ownedVerticesCount, neighborRank.first);
    }

    // wait until all aReceives are complete.
    for (auto &rqst : vertexNumberRequests) {
      rqst->wait();
    }

    // Exchange list of shared vertices with the neighbor
    // to store send requests.
    std::vector<com::PtrRequest> vertexListRequests;

    for (auto &receivingRank : sharedVerticesSendMap) {
      int  sendSize = receivingRank.second.size();
      auto request  = utils::MasterSlave::getCommunication()->aSend(sendSize, receivingRank.first);
      vertexListRequests.push_back(request);
      if (sendSize != 0) {
        auto request = utils::MasterSlave::getCommunication()->aSend(receivingRank.second, receivingRank.first);
        vertexListRequests.push_back(request);
      }
    }

    for (auto &neighborRank : sharedVerticesSendMap) {
      int receiveSize = 0;
      utils::MasterSlave::getCommunication()->receive(receiveSize, neighborRank.first);
      if (receiveSize != 0) {
        std::vector<int> receivedSharedVertices;
        utils::MasterSlave::getCommunication()->receive(receivedSharedVertices, neighborRank.first);
        sharedVerticesReceiveMap.insert(std::make_pair(neighborRank.first, receivedSharedVertices));
      }
    }

    // wait until all aSends are complete.
    for (auto &rqst : vertexListRequests) {
      rqst->wait();
    }

    // #5: Second round assignment according to the number of owned vertices

    /* In case that a vertex can be shared between two ranks, the rank with lower
       vertex count will own the vertex. 
       If both ranks have same vertex count, the lower rank will own the vertex.
    */

    for (int i = 0; i < sharedVerticesGlobalIDs.size(); i++) {
      bool owned = true;

      for (auto &sharingRank : sharedVerticesReceiveMap) {
        std::vector<int> vec = sharingRank.second;
        if (std::find(vec.begin(), vec.end(), sharedVerticesGlobalIDs[i]) != vec.end()) {
          if ((ownedVerticesCount > neighborRanksVertexCount[sharingRank.first]) ||
              (ownedVerticesCount == neighborRanksVertexCount[sharingRank.first] && utils::MasterSlave::getRank() > sharingRank.first)) {
            owned = false;

            // // Decide upon owners,
            // PRECICE_DEBUG("Decide owners, first round by rough load balancing");
            // // Provide a more descriptive error message if direct access was enabled
            // PRECICE_CHECK(!(ranksAtInterface == 0 && _allowDirectAccess),
            //               "After repartitioning of mesh \"{}\" all ranks are empty. "
            //               "Please check the dimensions of the provided bounding box "
            //               "(in \"setMeshAccessRegion\") and verify that it covers vertices "
            //               "in the mesh or check the definition of the provided meshes.",
            //               _mesh->getName());
            // PRECICE_ASSERT(ranksAtInterface != 0);
            // int localGuess = _mesh->getGlobalNumberOfVertices() / ranksAtInterface; // Guess for a decent load balancing
            // // First round: every slave gets localGuess vertices
            // for (Rank rank : utils::MasterSlave::allRanks()) {
            //   int counter = 0;
            //   for (size_t i = 0; i < slaveOwnerVecs[rank].size(); i++) {
            //     // Vertex has no owner yet and rank could be owner
            //     if (globalOwnerVec[slaveGlobalIDs[rank][i]] == 0 && slaveTags[rank][i] == 1) {
            //       slaveOwnerVecs[rank][i]                 = 1; // Now rank is owner
            //       globalOwnerVec[slaveGlobalIDs[rank][i]] = 1; // Vertex now has owner
            //       counter++;
            //       if (counter == localGuess)

            break;
          }
        }
      }
      tags[sharedVerticesLocalIDs[i]] = owned ? 1 : 0;
    }

    setOwnerInformation(tags);
    auto filteredVertices = std::count(tags.begin(), tags.end(), 0);
    if (filteredVertices)
      PRECICE_WARN("{} of {} vertices of mesh {} have been filtered out since they have no influence on the mapping.",
                   filteredVertices, _mesh->getGlobalNumberOfVertices(), _mesh->getName());
    // end of two-level initialization section
  } else {
    if (utils::MasterSlave::isSlave()) {
      int numberOfVertices = _mesh->vertices().size();
      utils::MasterSlave::getCommunication()->send(numberOfVertices, 0);

      if (numberOfVertices != 0) {
        PRECICE_DEBUG("Tag vertices, number of vertices {}", numberOfVertices);
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
        PRECICE_DEBUG("My tags: {}", tags);
        PRECICE_DEBUG("My global IDs: {}", globalIDs);
        PRECICE_DEBUG("Send tags and global IDs");
        utils::MasterSlave::getCommunication()->send(tags, 0);
        utils::MasterSlave::getCommunication()->send(globalIDs, 0);
        utils::MasterSlave::getCommunication()->send(atInterface, 0);

        PRECICE_DEBUG("Receive owner information");
        std::vector<int> ownerVec(numberOfVertices, -1);
        utils::MasterSlave::getCommunication()->receive(ownerVec, 0);
        PRECICE_DEBUG("My owner information: {}", ownerVec);
        PRECICE_ASSERT(ownerVec.size() == static_cast<std::size_t>(numberOfVertices));
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
      PRECICE_DEBUG("My tags: {}", slaveTags[0]);

      // receive slave data
      Rank ranksAtInterface = 0;
      if (masterAtInterface)
        ranksAtInterface++;

      for (Rank rank : utils::MasterSlave::allSlaves()) {
        int localNumberOfVertices = -1;
        utils::MasterSlave::getCommunication()->receive(localNumberOfVertices, rank);
        PRECICE_DEBUG("Rank {} has {} vertices.", rank, localNumberOfVertices);
        slaveOwnerVecs[rank].resize(localNumberOfVertices, 0);

        if (localNumberOfVertices != 0) {
          PRECICE_DEBUG("Receive tags from slave rank {}", rank);
          utils::MasterSlave::getCommunication()->receive(slaveTags[rank], rank);
          utils::MasterSlave::getCommunication()->receive(slaveGlobalIDs[rank], rank);
          PRECICE_DEBUG("Rank {} has tags {}", rank, slaveTags[rank]);
          PRECICE_DEBUG("Rank {} has global IDs {}", rank, slaveGlobalIDs[rank]);
          bool atInterface = false;
          utils::MasterSlave::getCommunication()->receive(atInterface, rank);
          if (atInterface)
            ranksAtInterface++;
        }
      }

      // Decide upon owners,
      PRECICE_DEBUG("Decide owners, first round by rough load balancing");
      // Provide a more descriptive error message if direct access was enabled
      PRECICE_CHECK(!(ranksAtInterface == 0 && _allowDirectAccess),
                    "After repartitioning of mesh \"{}\" all ranks are empty. "
                    "Please check the dimensions of the provided bounding box "
                    "(in \"setMeshAccessRegion\") and verify that it covers vertices "
                    "in the mesh or check the definition of the provided meshes.",
                    _mesh->getName());
      PRECICE_ASSERT(ranksAtInterface != 0);
      int localGuess = _mesh->getGlobalNumberOfVertices() / ranksAtInterface; // Guess for a decent load balancing
      // First round: every slave gets localGuess vertices
      for (Rank rank : utils::MasterSlave::allRanks()) {
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
      for (Rank rank : utils::MasterSlave::allRanks()) {
        for (size_t i = 0; i < slaveOwnerVecs[rank].size(); i++) {
          if (globalOwnerVec[slaveGlobalIDs[rank][i]] == 0 && slaveTags[rank][i] == 1) {
            slaveOwnerVecs[rank][i]                 = 1;
            globalOwnerVec[slaveGlobalIDs[rank][i]] = rank + 1;
          }
        }
      }

      // Send information back to slaves
      for (Rank rank : utils::MasterSlave::allSlaves()) {
        if (not slaveTags[rank].empty()) {
          PRECICE_DEBUG("Send owner information to slave rank {}", rank);
          utils::MasterSlave::getCommunication()->send(slaveOwnerVecs[rank], rank);
        }
      }
      // Master data
      PRECICE_DEBUG("My owner information: {}", slaveOwnerVecs[0]);
      setOwnerInformation(slaveOwnerVecs[0]);

#ifndef NDEBUG
      for (size_t i = 0; i < globalOwnerVec.size(); i++) {
        if (globalOwnerVec[i] == 0) {
          PRECICE_DEBUG("The Vertex with global index {} of mesh: {} was completely filtered out, since it has no influence on any mapping.",
                        i, _mesh->getName());
        }
      }
#endif
      auto filteredVertices = std::count(globalOwnerVec.begin(), globalOwnerVec.end(), 0);
      if (filteredVertices) {
        PRECICE_WARN("{} of {} vertices of mesh {} have been filtered out since they have no influence on the mapping.{}",
                     filteredVertices, _mesh->getGlobalNumberOfVertices(), _mesh->getName(),
                     _allowDirectAccess ? " Associated data values of the filtered vertices will be filled with zero values in order to "
                                          "provide valid data for other participants when reading data."
                                        : "");
      }
    }
  }
}

bool ReceivedPartition::isAnyProvidedMeshNonEmpty() const
{
  for (const auto &fromMapping : _fromMappings) {
    if (not fromMapping->getOutputMesh()->vertices().empty()) {
      return true;
    }
  }
  for (const auto &toMapping : _toMappings) {
    if (not toMapping->getInputMesh()->vertices().empty()) {
      return true;
    }
  }
  return false;
}

bool ReceivedPartition::hasAnyMapping() const
{
  return not(_fromMappings.empty() && _toMappings.empty());
}

void ReceivedPartition::tagMeshFirstRound()
{
  // We want to have every vertex within the box if we access the mesh directly
  if (_allowDirectAccess) {
    _mesh->tagAll();
    return;
  }

  for (const mapping::PtrMapping &fromMapping : _fromMappings) {
    fromMapping->tagMeshFirstRound();
  }
  for (const mapping::PtrMapping &toMapping : _toMappings) {
    toMapping->tagMeshFirstRound();
  }
}

void ReceivedPartition::tagMeshSecondRound()
{
  // We have already tagged every node in this case in the first round
  if (_allowDirectAccess) {
    return;
  }

  for (const mapping::PtrMapping &fromMapping : _fromMappings) {
    fromMapping->tagMeshSecondRound();
  }
  for (const mapping::PtrMapping &toMapping : _toMappings) {
    toMapping->tagMeshSecondRound();
  }
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

m2n::M2N &ReceivedPartition::m2n()
{
  PRECICE_ASSERT(_m2ns.size() == 1);
  return *_m2ns[0];
}

} // namespace partition
} // namespace precice

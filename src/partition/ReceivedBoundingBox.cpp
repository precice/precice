#include "partition/ReceivedBoundingBox.hpp"
#include <map>
#include <vector>
#include "com/CommunicateBoundingBox.hpp"
#include "com/Communication.hpp"
#include "m2n/M2N.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "utils/Helpers.hpp"
#include "utils/MasterSlave.hpp"

namespace precice
{
namespace partition
{

ReceivedBoundingBox::ReceivedBoundingBox(
    mesh::PtrMesh mesh, double safetyFactor)
    : Partition(mesh),
      _bb(mesh->getDimensions(), std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest())),
      _dimensions(mesh->getDimensions()),
      _safetyFactor(safetyFactor)
{
}

void ReceivedBoundingBox::communicateBoundingBox()
{
  TRACE();

  if (not utils::MasterSlave::_slaveMode) {
    _m2n->getMasterCommunication()->receive(_remoteParComSize, 0);

    // construct and initialize _remoteBBM
    mesh::Mesh::BoundingBox initialBB;
    for (int i = 0; i < _dimensions; i++) {
      initialBB.push_back(std::make_pair(-1, -1));
    }
    for (int remoteRank = 0; remoteRank < _remoteParComSize; remoteRank++) {
      _remoteBBM[remoteRank] = initialBB;
    }

    // master receives global_bb from other master
    com::CommunicateBoundingBox(_m2n->getMasterCommunication()).receiveBoundingBoxMap(_remoteBBM, 0);
  }
}

void ReceivedBoundingBox::computeBoundingBox()
{
  TRACE();

  /// @todo handle coupling mode (i.e. serial participant)

  prepareBoundingBox();

  if (utils::MasterSlave::_masterMode) { // Master
    assertion(utils::MasterSlave::_rank == 0);
    assertion(utils::MasterSlave::_size > 1);

    // broadcast _remoteBBM to all slaves
    utils::MasterSlave::_communication->broadcast(_remoteParComSize);
    com::CommunicateBoundingBox(utils::MasterSlave::_communication).broadcastSendBoundingBoxMap(_remoteBBM);

    std::map<int, std::vector<int>> connectionMap;
    std::vector<int> connectedRanksList;

    // connected ranks for master
    for (auto &remoteBB : _remoteBBM) {
      if (overlapping(_bb, remoteBB.second)) {
        _mesh->getConnectedRanks().push_back(remoteBB.first);
      }
    }
    if (_mesh->getConnectedRanks().size() != 0)
    {
     connectionMap[0] = _mesh->getConnectedRanks();
     connectedRanksList.push_back(0);
    }
      
    // receive connected ranks from slaves and add them to the connection map
    std::vector<int> slaveConnectedRanks;
    for (int rank = 1; rank < utils::MasterSlave::_size; rank++) {
      int connectedRanksSize = 0;
      utils::MasterSlave::_communication->receive(connectedRanksSize, rank);
      if (connectedRanksSize != 0) {
        connectedRanksList.push_back(rank);
        connectionMap[rank].push_back(-1);        
        utils::MasterSlave::_communication->receive(slaveConnectedRanks, rank);
        connectionMap[rank] = slaveConnectedRanks;
        slaveConnectedRanks.clear();
      }
    }

    // send connectionMap to other master
    _m2n->getMasterCommunication()->send(connectedRanksList, 0);
    if (connectionMap.size() != 0) { // @todo we need an error message here instea
      com::CommunicateBoundingBox(_m2n->getMasterCommunication()).sendConnectionMap(connectionMap, 0);
    } else
    {
      ERROR("This participant has no rank in the interface! Please check your test case and make sure that the mesh partition given to preCICE is loacted in the interface");
    }

  } else if (utils::MasterSlave::_slaveMode) {
    utils::MasterSlave::_communication->broadcast(_remoteParComSize, 0);

    // construct and initialize _remoteBBM
    mesh::Mesh::BoundingBox initialBB;
    for (int i = 0; i < _dimensions; i++) {
      initialBB.push_back(std::make_pair(-1, -1));
    }
    for (int remoteRank = 0; remoteRank < _remoteParComSize; remoteRank++) {
      _remoteBBM[remoteRank] = initialBB;
    }

    // receive _remoteBBM from master
    com::CommunicateBoundingBox(utils::MasterSlave::_communication).broadcastReceiveBoundingBoxMap(_remoteBBM);

    for (auto &remoteBB : _remoteBBM) {
      if (overlapping(_bb, remoteBB.second)) {
        _mesh->getConnectedRanks().push_back(remoteBB.first);
      }
    }    

    // send feedback size to master
    utils::MasterSlave::_communication->send((int) _mesh->getConnectedRanks().size(), 0);

    // to prevent sending empty vector!
    if (_mesh->getConnectedRanks().size() != 0)
      utils::MasterSlave::_communication->send(_mesh->getConnectedRanks(), 0);
  }
}

void ReceivedBoundingBox::communicate()
{
  // each rank receives max/min global vertex indexes from connected remote ranks
  _m2n->broadcastReceiveLCM(_maxVertexGlobalIndexDomain, *_mesh);
  
  // each rank receives mesh partition from connected ranks
  _m2n->broadcastReceiveLocalMesh(*_mesh);
}

void ReceivedBoundingBox::compute()
{
  if (not utils::MasterSlave::_slaveMode) {
    CHECK(_fromMapping.use_count() > 0 || _toMapping.use_count() > 0,
          "The received mesh " << _mesh->getName()
          << " needs a mapping, either from it, to it, or both. Maybe you don't want to receive this mesh at all?")
  }

  // (1) Bounding Box Filter

  INFO("Filter mesh " << _mesh->getName() << " by bounding-box");    
  mesh::Mesh filteredMesh("FilteredMesh", _dimensions, _mesh->isFlipNormals());
  filterMesh(filteredMesh, true);
    
  if ((_fromMapping.use_count() > 0 && _fromMapping->getOutputMesh()->vertices().size() > 0) ||
      (_toMapping.use_count() > 0 && _toMapping->getInputMesh()->vertices().size() > 0))
  {
      // this rank has vertices at the coupling interface
      // then, also the filtered mesh should still have vertices
    std::string msg = "The re-partitioning completely filtered out the mesh " + _mesh->getName() + " received on this rank at the coupling interface. "
      "Most probably, the coupling interfaces of your coupled participants do not match geometry-wise. "
      "Please check your geometry setup again. Small overlaps or gaps are no problem. "
      "If your geometry setup is correct and if you have very different mesh resolutions on both sides, increasing the safety-factor "
      "of the decomposition strategy might be necessary.";
    CHECK(filteredMesh.vertices().size() > 0, msg);
  }
  DEBUG("Bounding box filter, filtered from " << _mesh->vertices().size() << " vertices to " << filteredMesh.vertices().size() << " vertices.");
  _mesh->clear();
  _mesh->addMesh(filteredMesh);
  _mesh->computeState();
  
  // (2) Tag vertices 1st round (i.e. who could be owned by this rank)

  DEBUG("Tag vertices for filtering: 1st round.");
  // go to both meshes, vertex is tagged if already one mesh tags him
  if (_fromMapping.use_count() > 0)
    _fromMapping->tagMeshFirstRound();
  if (_toMapping.use_count() > 0)
    _toMapping->tagMeshFirstRound();
  
  // (3) Define which vertices are owned by this rank
  DEBUG("Create owner information.");
  createOwnerInformation();
  
  // (4) Tag vertices 2nd round (what should be filtered out)
  DEBUG("Tag vertices for filtering: 2nd round.");
  if (_fromMapping.use_count() > 0)
    _fromMapping->tagMeshSecondRound();
  if (_toMapping.use_count() > 0)
    _toMapping->tagMeshSecondRound();
  
  
  // (5) Filter mesh according to tag
  filteredMesh.clear();
  INFO("Filter mesh " << _mesh->getName() << " by mappings");
  filterMesh(filteredMesh, false);
  DEBUG("Mapping filter, filtered from " << _mesh->vertices().size() << " vertices to " << filteredMesh.vertices().size() << " vertices.");
  _mesh->clear();
  _mesh->addMesh(filteredMesh);
  _mesh->computeState();
  
  // (6) Compute and feedback local communication map
  INFO("Feedback Communicatin Map "); 
  std::map<int, std::vector<int>> localCommunicationMap;

  /*
   * The following nested loops fill in the localCommunicationbMap which shows this local
   * rank needs which vertices from which rank of the other particpant
   * for example:
   * 4->{6, 8, 9}
   * 8->{45, 49, 55}
   * this rank needs vertices 6 ,8 and 9 from rank 4 and vertices 45, 49 and 55
   * from rank 8
   * If global index of a vertex is between min/max global index of a specific 
   * remore rank, this vertex belongs to that rank
   */
  for (auto &remoteVertex : _mesh->vertices()) {
    for (auto &remoteRank : _maxVertexGlobalIndexDomain) {
      if (remoteVertex.getGlobalIndex() <= remoteRank.second[1] && remoteVertex.getGlobalIndex() >= remoteRank.second[0]) {
        localCommunicationMap[remoteRank.first].push_back(remoteVertex.getGlobalIndex());
        break;
      }
    }
  }
  
  _m2n->broadcastSendLCM(localCommunicationMap, *_mesh);   
}

bool ReceivedBoundingBox::overlapping(mesh::Mesh::BoundingBox currentBB, mesh::Mesh::BoundingBox receivedBB)
{
  /*
   * Here two bounding boxes are compared to check whether they overlap or not!
   * Comparison is done for all dimensions and, of course, to have a proper overlap,
   * each dimension must overlap.
   * We need to check if first AND second is smaller than first of the other BB to prevent false negatives
   * due to empty bounding boxes.
   */
  for (int i = 0; i < _dimensions; i++) {
    if ((currentBB[i].first < receivedBB[i].first && currentBB[i].second < receivedBB[i].first) ||
        (receivedBB[i].first < currentBB[i].first && receivedBB[i].second < currentBB[i].first)) {
      return false;
    }
  }
  return true;
}

void ReceivedBoundingBox::prepareBoundingBox()
{
  TRACE(_safetyFactor);

  _bb.resize(_dimensions, std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()));

  //create BB around both "other" meshes
  if (_fromMapping.use_count() > 0) {
    auto other_bb = _fromMapping->getOutputMesh()->getBoundingBox();
    for (int d = 0; d < _dimensions; d++) {
      if (_bb[d].first > other_bb[d].first)
        _bb[d].first = other_bb[d].first;
      if (_bb[d].second < other_bb[d].second)
        _bb[d].second = other_bb[d].second;
    }
  }
  if (_toMapping.use_count() > 0) {
    auto other_bb = _toMapping->getInputMesh()->getBoundingBox();
    for (int d = 0; d < _dimensions; d++) {
      if (_bb[d].first > other_bb[d].first)
        _bb[d].first = other_bb[d].first;
      if (_bb[d].second < other_bb[d].second)
        _bb[d].second = other_bb[d].second;
    }
  }

  //enlarge BB
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

void ReceivedBoundingBox:: filterMesh(mesh::Mesh& filteredMesh, const bool filterByBB) {
  TRACE(filterByBB);

  DEBUG("Bounding mesh. #vertices: " << _mesh->vertices().size()
               <<", #edges: " << _mesh->edges().size()
               <<", #triangles: " << _mesh->triangles().size() << ", rank: " << utils::MasterSlave::_rank);

  std::map<int, mesh::Vertex*> vertexMap;
  std::map<int, mesh::Edge*> edgeMap;
  int vertexCounter = 0;

  for (const mesh::Vertex& vertex : _mesh->vertices())
  {
    if ((filterByBB && isVertexInBB(vertex)) || (not filterByBB && vertex.isTagged()))
    {
      mesh::Vertex& v = filteredMesh.createVertex(vertex.getCoords());
      v.setGlobalIndex(vertex.getGlobalIndex());
      if(vertex.isTagged()) v.tag();
      v.setOwner(vertex.isOwner());
      vertexMap[vertex.getID()] = &v;
    }
    vertexCounter++;
  }

  // Add all edges formed by the contributing vertices
  for (mesh::Edge& edge : _mesh->edges())
  {
    int vertexIndex1 = edge.vertex(0).getID();
    int vertexIndex2 = edge.vertex(1).getID();
    if (utils::contained(vertexIndex1, vertexMap) &&
        utils::contained(vertexIndex2, vertexMap)) {
      mesh::Edge& e = filteredMesh.createEdge(*vertexMap[vertexIndex1], *vertexMap[vertexIndex2]);
      edgeMap[edge.getID()] = &e;
    }
  }

  // Add all triangles formed by the contributing edges
  if (_dimensions==3)
  {
    for (mesh::Triangle& triangle : _mesh->triangles() )
    {
      int edgeIndex1 = triangle.edge(0).getID();
      int edgeIndex2 = triangle.edge(1).getID();
      int edgeIndex3 = triangle.edge(2).getID();
      if (utils::contained(edgeIndex1, edgeMap) &&
          utils::contained(edgeIndex2, edgeMap) &&
          utils::contained(edgeIndex3, edgeMap))
      {
        filteredMesh.createTriangle(*edgeMap[edgeIndex1],*edgeMap[edgeIndex2],*edgeMap[edgeIndex3]);
      }
    }
  }

  DEBUG("Filtered mesh. #vertices: " << filteredMesh.vertices().size()
               <<", #edges: " << filteredMesh.edges().size()
               <<", #triangles: " << filteredMesh.triangles().size() << ", rank: " << utils::MasterSlave::_rank);
}

bool ReceivedBoundingBox::isVertexInBB(const mesh::Vertex& vertex) {
  for (int d=0; d<_dimensions; d++) {
    if (vertex.getCoords()[d] < _bb[d].first || vertex.getCoords()[d] > _bb[d].second ) {
      return false;
    }
  }
  return true;
}

void ReceivedBoundingBox::createOwnerInformation()
{ 

  int numberOfVertices = _mesh->vertices().size();    
  if (numberOfVertices != 0)
  {
    std::vector<int> tags(numberOfVertices, -1);
    std::vector<int> globalIDs(numberOfVertices, -1);
    for (int i = 0; i < numberOfVertices; i++)
    {
      globalIDs[i] = _mesh->vertices()[i].getGlobalIndex();
      if (_mesh->vertices()[i].isTagged())
      {
        tags[i]     = 1;
      }
      else
      {
        tags[i] = 0;
      }
    }
    DEBUG("My tags: " << tags);
    DEBUG("My global IDs: " << globalIDs); 
  }    
}

} // namespace partition
} // namespace precice

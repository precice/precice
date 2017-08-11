#include "partition/ReceivedPartition.hpp"
#include "com/CommunicateMesh.hpp"
#include "utils/MasterSlave.hpp"
#include "m2n/M2N.hpp"
#include "utils/EventTimings.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "utils/Helpers.hpp"

using precice::utils::Event;

namespace precice {
namespace partition {

logging::Logger ReceivedPartition:: _log ( "precice::partition::ReceivedPartition" );

ReceivedPartition::ReceivedPartition
(
    bool filterFirst, int dimensions, double safetyFactor )
:
  Partition (),
  _filterFirst(filterFirst),
  _bb(dimensions, std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest())),
  _dimensions(dimensions),
  _safetyFactor(safetyFactor)
{}

void ReceivedPartition::communicate()
{
  TRACE();
  //TODO sind globalIDs beim gathern vergeben worden? generell überlegen, wo das passieren soll
  INFO("Receive global mesh " << _mesh->getName());
  if (not utils::MasterSlave::_slaveMode) {
    assertion ( _mesh->vertices().size() == 0 );
    com::CommunicateMesh(_m2n->getMasterCommunication()).receiveMesh ( *_mesh, 0 );
    _mesh->setGlobalNumberOfVertices(_mesh->vertices().size()); //TODO muss hier?
  }
}

void ReceivedPartition::compute()
{
  TRACE();

  //TODO coupling mode abfangen


  // (1) Bounding-Box-Filter

  if(_filterFirst){ //pre-filter-post-filter

    INFO("Pre-filter mesh " << _mesh->getName() << " by bounding-box");
    Event e("pre-filter mesh by bounding box");

    if (utils::MasterSlave::_slaveMode) {
      prepareBoundingBox();
      com::CommunicateMesh(utils::MasterSlave::_communication).sendBoundingBox (_bb, 0);
      com::CommunicateMesh(utils::MasterSlave::_communication).receiveMesh (*_mesh, 0);

      if((_fromMapping.use_count()>0 && _fromMapping->getOutputMesh()->vertices().size()>0) ||
         (_toMapping.use_count()>0 && _toMapping->getInputMesh()->vertices().size()>0)){
           // this rank has vertices at the coupling interface
           // then, also the filtered mesh should still have vertices
        std::string msg = "The re-partitioning completely filtered out the mesh received on this rank at the coupling interface. "
            "Most probably, the coupling interfaces of your coupled participants do not match geometry-wise. "
            "Please check your geometry setup again. Small overlaps or gaps are no problem. "
            "If your geometry setup is correct and if you have very different mesh resolutions on both sides, increasing the safety-factor "
            "of the decomposition strategy might be necessary.";
        CHECK(_mesh->vertices().size()>0, msg);
      }

    }
    else{ // Master
      assertion(utils::MasterSlave::_rank==0);
      assertion(utils::MasterSlave::_size>1);

      for (int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++) {
        com::CommunicateMesh(utils::MasterSlave::_communication).receiveBoundingBox ( _bb, rankSlave);

        DEBUG("From slave " << rankSlave << ", bounding mesh: " << _bb[0].first
                     << ", " << _bb[0].second << " and " << _bb[1].first << ", " << _bb[1].second);
        mesh::Mesh slaveMesh("SlaveMesh", _dimensions, _mesh->isFlipNormals());
        std::vector<int> boundingVertexDistribution = filterMesh(slaveMesh, true); //TODO return here still needed?
        com::CommunicateMesh(utils::MasterSlave::_communication).sendMesh ( slaveMesh, rankSlave );
      }

      // Now also filter the remaining master mesh
      prepareBoundingBox();
      mesh::Mesh filteredMesh("FilteredMesh", _dimensions, _mesh->isFlipNormals());
      std::vector<int> tmpVertexPostitions = filterMesh(filteredMesh, true); //TODO still needed?
      _mesh->clear();
      _mesh->addMesh(filteredMesh);
      _mesh->computeState();
      DEBUG("Master mesh after filtering, #vertices " << _mesh->vertices().size());

      if((_fromMapping.use_count()>0 && _fromMapping->getOutputMesh()->vertices().size()>0) ||
         (_toMapping.use_count()>0 && _toMapping->getInputMesh()->vertices().size()>0)){
           // this rank has vertices at the coupling interface
           // then, also the filtered mesh should still have vertices
        std::string msg = "The re-partitioning completely filtered out the mesh received on this rank at the coupling interface. "
            "Most probably, the coupling interfaces of your coupled participants do not match geometry-wise. "
            "Please check your geometry setup again. Small overlaps or gaps are no problem. "
            "If your geometry setup is correct and if you have very different mesh resolutions on both sides, increasing the safety-factor "
            "of the decomposition strategy might be necessary.";
        CHECK(_mesh->vertices().size()>0, msg);
      }

    }
  }
  else{ //broadcast-filter
    INFO("Broadcast mesh " << _mesh->getName() );
    Event e1("broadcast mesh");

    if (utils::MasterSlave::_slaveMode) {
      com::CommunicateMesh(utils::MasterSlave::_communication).broadcastReceiveMesh (*_mesh);
    }
    else{ // Master
      assertion(utils::MasterSlave::_rank==0);
      assertion(utils::MasterSlave::_size>1);
      com::CommunicateMesh(utils::MasterSlave::_communication).broadcastSendMesh (*_mesh);
    }

    e1.stop();

    INFO("Filter mesh " << _mesh->getName() << " by bounding-box");
    Event e2("filter mesh by bounding box");

    prepareBoundingBox();
    mesh::Mesh filteredMesh("FilteredMesh", _dimensions, _mesh->isFlipNormals());
    std::vector<int> tmpVertexPostitions = filterMesh(filteredMesh, true); //TODO vll brauch man vertexPos hier nicht mehr

    if((_fromMapping.use_count()>0 && _fromMapping->getOutputMesh()->vertices().size()>0) ||
       (_toMapping.use_count()>0 && _toMapping->getInputMesh()->vertices().size()>0)){
         // this rank has vertices at the coupling interface
         // then, also the filtered mesh should still have vertices
      std::string msg = "The re-partitioning completely filtered out the mesh received on this rank at the coupling interface. "
          "Most probably, the coupling interfaces of your coupled participants do not match geometry-wise. "
          "Please check your geometry setup again. Small overlaps or gaps are no problem. "
          "If your geometry setup is correct and if you have very different mesh resolutions on both sides, increasing the safety-factor "
          "of the decomposition strategy might be necessary.";
      CHECK(filteredMesh.vertices().size()>0, msg);
    }

    DEBUG("Bounding box filter, filtered from " << _mesh->vertices().size() << " vertices to " << filteredMesh.vertices().size() << " vertices.");
    _mesh->clear();
    _mesh->addMesh(filteredMesh);
    _mesh->computeState();
    e2.stop();
  }

  // (2) Tag vertices 1st round (i.e. who could be owned by this rank)
  DEBUG("Tag vertices for filtering: 1st round.");
  // go to both meshes, vertex is tagged if already one mesh tags him
  if (_fromMapping.use_count() > 0) _fromMapping->tagMeshFirstRound();
  if (_toMapping.use_count() > 0) _toMapping->tagMeshFirstRound();


  // (3) Define which vertices are owned by this rank
  DEBUG("Create owner information.");
  createOwnerInformation();

  // (4) Tag vertices 2nd round (what should be filtered out)
  DEBUG("Tag vertices for filtering: 2nd round.");
  if (_fromMapping.use_count() > 0) _fromMapping->tagMeshSecondRound();
  if (_toMapping.use_count() > 0) _toMapping->tagMeshSecondRound();

  // (5) Filter mesh according to tag

  INFO("Filter mesh " << _mesh->getName() << " by mappings");
  Event e5("filter mesh by mappings");
  mesh::Mesh filteredMesh("FilteredMesh", _dimensions, _mesh->isFlipNormals());
  std::vector<int> tmpVertexPostitions = filterMesh(filteredMesh, false);
  DEBUG("Mapping filter, filtered from " << _mesh->vertices().size() << " vertices to " << filteredMesh.vertices().size() << " vertices.");
  _mesh->clear();
  _mesh->addMesh(filteredMesh);
  _mesh->computeState();
  e5.stop();

  // (6) Compute distribution

  //TODO old feedback step

}

std::vector<int> ReceivedPartition:: filterMesh(mesh::Mesh& filteredMesh, const bool filterByBB){
  TRACE();

  DEBUG("Bounding mesh. #vertices: " << _mesh->vertices().size()
               <<", #edges: " << _mesh->edges().size()
               <<", #triangles: " << _mesh->triangles().size() << ", rank: " << utils::MasterSlave::_rank);

  std::vector<int> vertexPositions;
  std::map<int, mesh::Vertex*> vertexMap;
  std::map<int, mesh::Edge*> edgeMap;
  int vertexCounter = 0;

  for (const mesh::Vertex& vertex : _mesh->vertices()) {

    if ((filterByBB && isVertexInBB(vertex)) || (not filterByBB && vertex.isTagged())){
      mesh::Vertex& v = filteredMesh.createVertex(vertex.getCoords());
      v.setGlobalIndex(vertex.getGlobalIndex());
      vertexPositions.push_back(vertexCounter);
      vertexMap[vertex.getID()] = &v;
    }
    vertexCounter++;
  }

  // Add all edges formed by the contributing vertices
  for (mesh::Edge& edge : _mesh->edges()) {
    int vertexIndex1 = edge.vertex(0).getID();
    int vertexIndex2 = edge.vertex(1).getID();
    if (utils::contained(vertexIndex1, vertexMap) &&
        utils::contained(vertexIndex2, vertexMap)) {
      mesh::Edge& e = filteredMesh.createEdge(*vertexMap[vertexIndex1], *vertexMap[vertexIndex2]);
      edgeMap[edge.getID()] = &e;
    }
  }

  // Add all triangles formed by the contributing edges
  if (_dimensions==3) {
    for (mesh::Triangle& triangle : _mesh->triangles() ) {
      int edgeIndex1 = triangle.edge(0).getID();
      int edgeIndex2 = triangle.edge(1).getID();
      int edgeIndex3 = triangle.edge(2).getID();
      if (utils::contained(edgeIndex1, edgeMap) &&
          utils::contained(edgeIndex2, edgeMap) &&
          utils::contained(edgeIndex3, edgeMap)) {
        filteredMesh.createTriangle(*edgeMap[edgeIndex1],*edgeMap[edgeIndex2],*edgeMap[edgeIndex3]);
      }
    }
  }

  DEBUG("Filtered mesh. #vertices: " << filteredMesh.vertices().size()
               <<", #edges: " << filteredMesh.edges().size()
               <<", #triangles: " << filteredMesh.triangles().size() << ", rank: " << utils::MasterSlave::_rank);

  return vertexPositions;
}

void ReceivedPartition::prepareBoundingBox(){

  _bb.resize(_dimensions, std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()));

  //create BB around both "other" meshes
  if (_fromMapping.use_count()>0) {
    auto other_bb = _fromMapping->getOutputMesh()->getBoundingBox();
    for (int d=0; d < _dimensions; d++) {
      if (_bb[d].first > other_bb[d].first) _bb[d].first = other_bb[d].first;
      if (_bb[d].second < other_bb[d].second) _bb[d].second = other_bb[d].second;
    }
  }
  if (_toMapping.use_count()>0) {
    auto other_bb = _toMapping->getInputMesh()->getBoundingBox();
    for (int d=0; d<_dimensions; d++) {
      if (_bb[d].first > other_bb[d].first) _bb[d].first = other_bb[d].first;
      if (_bb[d].second < other_bb[d].second) _bb[d].second = other_bb[d].second;
    }
  }

  //enlarge BB
  assertion(_safetyFactor>=0.0);
  for (int d=0; d<_dimensions; d++) {
    if (_bb[d].second > _bb[d].first){
      double sideLength = _bb[d].second - _bb[d].first;
      _bb[d].second += _safetyFactor * sideLength;
      _bb[d].first -= _safetyFactor * sideLength;
    }
  }
}

bool ReceivedPartition::isVertexInBB(const mesh::Vertex& vertex){
  for (int d=0; d<_dimensions; d++) {
    if (vertex.getCoords()[d] < _bb[d].first || vertex.getCoords()[d] > _bb[d].second ) {
      return false;
    }
  }
  return true;
}

}} // namespace precice, partition

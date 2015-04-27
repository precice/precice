#include "CommunicatedGeometry.hpp"
#include "com/CommunicateMesh.hpp"
#include "com/Communication.hpp"
#include "m2n/M2N.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/Globals.hpp"
#include "utils/Helpers.hpp"
#include "utils/EventTimings.hpp"

using precice::utils::Event;

namespace precice {
namespace geometry {

tarch::logging::Log CommunicatedGeometry:: _log ( "precice::geometry::CommunicatedGeometry" );

CommunicatedGeometry:: CommunicatedGeometry
(
  const utils::DynVector& offset,
  const std::string&      accessor,
  const std::string&      provider,
  int                     dimensions)
  :
  Geometry ( offset ),
  _accessorName ( accessor ),
  _providerName ( provider ),
  _receivers (),
  _dimensions(dimensions),
  _boundingFromMapping(),
  _boundingToMapping(),
  _bb(),
  _safetyGap(0),
  _safetyFactor(-1.0)
{
  preciceTrace2 ( "CommunicatedGeometry()", accessor, provider );
}

void CommunicatedGeometry:: addReceiver
(
  const std::string&     receiver,
  m2n::M2N::SharedPointer m2n)
{
  preciceTrace1 ( "addReceiver()", receiver );
  assertion ( m2n.get() != NULL );
  preciceCheck ( ! utils::contained(receiver, _receivers),
                 "addReceiver()", "Receiver \"" << receiver
                 << "\" has been added already to communicated geometry!" );
  preciceCheck ( receiver != _providerName, "addReceiver()",
                 "Receiver \"" << receiver << "\" cannot be the same as "
                 << "provider in communicated geometry!" );
  _receivers[receiver] = m2n;
}

void CommunicatedGeometry:: specializedCreate
(
  mesh::Mesh& seed )
{
  preciceTrace1 ( "specializedCreate()", seed.getName() );
  preciceCheck ( _receivers.size() > 0, "specializedCreate()",
                 "No receivers specified for communicated geometry to create "
                 << "mesh \"" << seed.getName() << "\"!" );
  if ( _accessorName == _providerName ) {
    sendMesh(seed);
  }
  else if ( utils::contained(_accessorName, _receivers) ) {
    receiveMesh(seed);
  }
  else {
    preciceError( "specializedCreate()", "Participant \"" << _accessorName
                  << "\" uses a communicated geometry to create mesh \""
                  << seed.getName()
                  << "\" but is neither provider nor receiver!" );
  }
}


void CommunicatedGeometry:: sendMesh(
  mesh::Mesh& seed)
{
  preciceTrace1 ( "sendMesh()", utils::MasterSlave::_rank );
  // Temporary globalMesh such that the master also keeps his local mesh (seed)
  mesh::Mesh globalMesh(seed.getName(), _dimensions, seed.isFlipNormals());

  if( not utils::MasterSlave::_slaveMode ){
    globalMesh.addMesh(seed); //add local master mesh to global mesh
  }

  // Gather Mesh
  preciceInfo("sendMesh()", "Gather mesh " << seed.getName() );
  if (utils::MasterSlave::_slaveMode || utils::MasterSlave::_masterMode ) {
    Event e("gather mesh");
    if (utils::MasterSlave::_slaveMode) {
      com::CommunicateMesh(utils::MasterSlave::_communication).sendMesh( seed, 0 );
    }
    else{ // Master
      assertion(utils::MasterSlave::_rank==0);
      assertion(utils::MasterSlave::_size>1);

      int numberOfVertices = 0;
      // Vertices of master mesh part do already exist
      for (int i = 0; i < seed.vertices().size(); i++) {
        seed.getVertexDistribution()[0].push_back(numberOfVertices);
        numberOfVertices++;
      }

      for (int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++) {
        int vertexCount1 = globalMesh.vertices().size();
        com::CommunicateMesh(utils::MasterSlave::_communication).receiveMesh ( globalMesh, rankSlave);
        int vertexCount2 = globalMesh.vertices().size();
        int vertexCountDiff = vertexCount2 - vertexCount1;
        preciceDebug("Received sub-mesh, from slave: " << rankSlave <<", vertexCount: " << vertexCountDiff);
        for (int i = 0; i < vertexCountDiff; i++) {
          seed.getVertexDistribution()[rankSlave].push_back(numberOfVertices);
          numberOfVertices++;
        }
      }
      seed.setGlobalNumberOfVertices(numberOfVertices);
    }
  }

  // Send (global) Mesh
  preciceInfo("sendMesh()", "Send global mesh " << seed.getName());
  Event e("send global mesh");
  if (not utils::MasterSlave::_slaveMode) {
    preciceCheck ( globalMesh.vertices().size() > 0,
                   "specializedCreate()", "Participant \"" << _accessorName
                   << "\" provides an invalid (possibly empty) mesh \""
                   << globalMesh.getName() << "\"!" );
    for (auto &pair : _receivers) {
      com::CommunicateMesh(pair.second->getMasterCommunication()).sendMesh ( globalMesh, 0 );
    }
  }
}

void CommunicatedGeometry:: receiveMesh(
  mesh::Mesh& seed)
{
  preciceInfo("receiveMesh()", "Receive global mesh " << seed.getName() );
  if (not utils::MasterSlave::_slaveMode) {
    Event e("receive global mesh");
    assertion ( seed.vertices().size() == 0 );
    assertion ( utils::contained(_accessorName, _receivers) );
    m2n::M2N::SharedPointer m2n ( _receivers[_accessorName] );
    com::CommunicateMesh(m2n->getMasterCommunication()).receiveMesh ( seed, 0 );
  }
  if (utils::MasterSlave::_slaveMode || utils::MasterSlave::_masterMode){
    scatterMesh(seed);
  }
}

void CommunicatedGeometry:: scatterMesh(
  mesh::Mesh& seed)
{
  preciceTrace1 ( "scatterMesh()", utils::MasterSlave::_rank );
  using tarch::la::raw;

  // Send bounding boxes from all slaves to master, then send according part of the global mesh back
  preciceInfo("scatterMesh()", "Scatter bounding-box-filtered meshes for " << seed.getName() );
  Event e1("scatter bounding-box-filtered meshes");
  std::map<int,std::vector<int> > boundingVertexDistribution;
  if (utils::MasterSlave::_slaveMode) {
    mesh::Mesh::BoundingBox bb = mesh::Mesh::BoundingBox (_dimensions,
                                                          std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()));
    mergeBoundingBoxes(bb);
    com::CommunicateMesh(utils::MasterSlave::_communication).sendBoundingBox (bb, 0);
    com::CommunicateMesh(utils::MasterSlave::_communication).receiveMesh (seed, 0);
  }
  else{ // Master
    assertion(utils::MasterSlave::_rank==0);
    assertion(utils::MasterSlave::_size>1);
    seed.setGlobalNumberOfVertices(seed.vertices().size());
    for (int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++) {
      _bb = mesh::Mesh::BoundingBox (_dimensions, std::make_pair(0.0,0.0));
      com::CommunicateMesh(utils::MasterSlave::_communication).receiveBoundingBox ( _bb, rankSlave);
      for (int d=0; d<_dimensions; d++) {
        if (_bb[d].second > _bb[d].first && _safetyGap < _bb[d].second - _bb[d].first)
          _safetyGap = _bb[d].second - _bb[d].first;
      }
      assertion(_safetyFactor>=0.0);
      _safetyGap *= _safetyFactor;
      preciceDebug("From slave " << rankSlave << ", bounding mesh: " << _bb[0].first
                   << ", " << _bb[0].second << " and " << _bb[1].first << ", " << _bb[1].second);
      mesh::Mesh slaveMesh("SlaveMesh", _dimensions, seed.isFlipNormals());
      boundingVertexDistribution[rankSlave] = filterMesh(seed, slaveMesh, false);
      com::CommunicateMesh(utils::MasterSlave::_communication).sendMesh ( slaveMesh, rankSlave );
    }
    // Now also filter the remaining master mesh
    _bb = mesh::Mesh::BoundingBox (_dimensions,
                                   std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()));
    mergeBoundingBoxes(_bb);
    for (int d=0; d < _dimensions; d++) {
      if (_bb[d].second > _bb[d].first && _safetyGap < _bb[d].second - _bb[d].first)
        _safetyGap = _bb[d].second - _bb[d].first;
    }
    _safetyGap *= _safetyFactor;
    mesh::Mesh filteredMesh("FilteredMesh", _dimensions, seed.isFlipNormals());
    boundingVertexDistribution[0] = filterMesh(seed, filteredMesh, false);
    seed.clear();
    seed.addMesh(filteredMesh);
    preciceDebug("Master mesh after filtering, #vertices " << seed.vertices().size());
  }
  e1.stop();

  preciceInfo("scatterMesh()", "Compute bounding mappings for mesh " << seed.getName() );
  Event e2("compute bounding mappings and filter");
  seed.computeState();
  computeBoundingMappings();

  preciceInfo("scatterMesh()", "Filter mesh " << seed.getName() );
  mesh::Mesh filteredMesh("FilteredMesh", _dimensions, seed.isFlipNormals());
  auto filteredVertexPositions = filterMesh(seed, filteredMesh, true);
  seed.clear();
  seed.addMesh(filteredMesh);
  clearBoundingMappings();
  e2.stop();

  int numberOfVertices = filteredVertexPositions.size();

  preciceInfo("scatterMesh()", "Gather vertex distribution for mesh " << seed.getName() );
  Event e3("gather vertex distribution");
  if (utils::MasterSlave::_slaveMode) {
    utils::MasterSlave::_communication->send(numberOfVertices,0);
    if (numberOfVertices!=0) {
      utils::MasterSlave::_communication->send(raw(filteredVertexPositions),numberOfVertices,0);
    }
  }
  else { // Master
    assertion(utils::MasterSlave::_rank==0);
    assertion(utils::MasterSlave::_size>1);

    //we need to merge the 2 filtering steps, each slave only holds local IDs
    std::vector<int> globalVertexIDs(numberOfVertices, -1);
    for (int i=0; i < numberOfVertices; i++) {
      globalVertexIDs[i] = boundingVertexDistribution[0][filteredVertexPositions[i]];
    }
    seed.getVertexDistribution()[0] = globalVertexIDs;

    for (int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++){
      int numberOfVertices = -1;
      utils::MasterSlave::_communication->receive(numberOfVertices,rankSlave);
      std::vector<int> slaveVertexIDs(numberOfVertices,-1);
      if (numberOfVertices!=0) {
        utils::MasterSlave::_communication->receive(raw(slaveVertexIDs),numberOfVertices,rankSlave);
      }

      //we need to merge the 2 filtering steps, each slave only holds local IDs
      std::vector<int> globalVertexIDs(numberOfVertices,-1);
      for(int i=0;i<numberOfVertices;i++){
        globalVertexIDs[i] = boundingVertexDistribution[rankSlave][slaveVertexIDs[i]];
      }
      seed.getVertexDistribution()[rankSlave] = globalVertexIDs;
    }
  }
}

void CommunicatedGeometry:: setBoundingFromMapping(
  mapping::PtrMapping mapping)
{
  _boundingFromMapping = mapping;
}

void CommunicatedGeometry:: setBoundingToMapping(
  mapping::PtrMapping mapping)
{
  _boundingToMapping = mapping;
}

void CommunicatedGeometry:: setSafetyFactor(
  double safetyFactor)
{
  _safetyFactor = safetyFactor;
}

void CommunicatedGeometry:: computeBoundingMappings()
{
  if (_boundingFromMapping.use_count() > 0) {
    _boundingFromMapping->computeMapping();
  }
  if (_boundingToMapping.use_count() > 0) {
    _boundingToMapping->computeMapping();
  }
}

void CommunicatedGeometry:: clearBoundingMappings()
{
  if (_boundingFromMapping.use_count() > 0) {
    _boundingFromMapping->clear();
  }
  if (_boundingToMapping.use_count() > 0) {
    _boundingToMapping->clear();
  }
}

void CommunicatedGeometry:: mergeBoundingBoxes(mesh::Mesh::BoundingBox& bb){
  if (_boundingFromMapping.use_count()>0) {
    auto bb1 = _boundingFromMapping->getOutputMesh()->getBoundingBox();
    for (int d=0; d < _dimensions; d++) {
      if (bb[d].first > bb1[d].first) bb[d].first = bb1[d].first;
      if (bb[d].second < bb1[d].second) bb[d].second = bb1[d].second;
    }
  }
  if (_boundingToMapping.use_count()>0) {
    auto bb2 = _boundingToMapping->getInputMesh()->getBoundingBox();
    for (int d=0; d<_dimensions; d++) {
      if (bb[d].first > bb2[d].first) bb[d].first = bb2[d].first;
      if (bb[d].second < bb2[d].second) bb[d].second = bb2[d].second;
    }
  }
}

bool CommunicatedGeometry:: doesVertexContribute(
  const mesh::Vertex& vertex, bool filterByMapping)
{
  if (filterByMapping) {
    //works as easy as this since only read-consistent and write-conservative are allowed
    assertion(_boundingFromMapping.use_count()>0 || _boundingToMapping.use_count()>0);
    bool exit = false;
    if (_boundingFromMapping.use_count() > 0) {
      exit = exit || _boundingFromMapping->doesVertexContribute(vertex.getID());
    }
    if (_boundingToMapping.use_count() > 0) {
      exit = exit || _boundingToMapping->doesVertexContribute(vertex.getID());
    }
    return exit;
  }
  else { //filter by bounding box
    for (int d=0; d<_dimensions; d++) {
      if (vertex.getCoords()[d] < _bb[d].first - _safetyGap || vertex.getCoords()[d] > _bb[d].second + _safetyGap) {
        return false;
      }
    }
    return true;
  }
}

std::vector<int> CommunicatedGeometry:: filterMesh(mesh::Mesh& seed, mesh::Mesh& filteredMesh, bool filterByMapping){
  preciceTrace1 ( "filterMesh()", utils::MasterSlave::_rank );

  preciceDebug("Bounding mesh. #vertices: " << seed.vertices().size()
               <<", #edges: " << seed.edges().size()
               <<", #triangles: " << seed.triangles().size() << ", rank: " << utils::MasterSlave::_rank);

  std::vector<int> vertexPositions;
  std::map<int, mesh::Vertex*> vertexMap;
  std::map<int, mesh::Edge*> edgeMap;
  int vertexCounter = 0;

  for (const mesh::Vertex& vertex : seed.vertices()) {
    if (doesVertexContribute(vertex, filterByMapping)){
      mesh::Vertex& v = filteredMesh.createVertex(vertex.getCoords());
      vertexPositions.push_back(vertexCounter);
      vertexMap[vertex.getID()] = &v;
    }
    vertexCounter++;
  }

  // Add all edges formed by the contributing vertices
  for (mesh::Edge& edge : seed.edges()) {
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
    for (mesh::Triangle& triangle : seed.triangles() ) {
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

  preciceDebug("Filtered mesh. #vertices: " << filteredMesh.vertices().size()
               <<", #edges: " << filteredMesh.edges().size()
               <<", #triangles: " << filteredMesh.triangles().size() << ", rank: " << utils::MasterSlave::_rank);

  return vertexPositions;
}

}} // namespace precice, geometry

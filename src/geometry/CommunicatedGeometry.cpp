// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "CommunicatedGeometry.hpp"
#include "com/CommunicateMesh.hpp"
#include "com/Communication.hpp"
#include "m2n/GlobalCommunication.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/Globals.hpp"
#include "utils/Helpers.hpp"

namespace precice {
namespace geometry {

tarch::logging::Log CommunicatedGeometry:: _log ( "precice::geometry::CommunicatedGeometry" );

CommunicatedGeometry:: CommunicatedGeometry
(
  const utils::DynVector&  offset,
  const std::string&       accessor,
  const std::string&       provider,
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
  _safetyGap(0)
{
  preciceTrace2 ( "CommunicatedGeometry()", accessor, provider );
}

void CommunicatedGeometry:: addReceiver
(
  const std::string&     receiver,
  m2n::PtrGlobalCommunication com )
{
  preciceTrace1 ( "addReceiver()", receiver );
  assertion ( com.get() != NULL );
  preciceCheck ( ! utils::contained(receiver, _receivers),
                 "addReceiver()", "Receiver \"" << receiver
                 << "\" has been added already to communicated geometry!" );
  preciceCheck ( receiver != _providerName, "addReceiver()",
                 "Receiver \"" << receiver << "\" cannot be the same as "
                 << "provider in communicated geometry!" );
  _receivers[receiver] = com;
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
  // temporary globalMesh such that the master also keeps his local mesh (seed)
  mesh::Mesh globalMesh("GlobalMesh", _dimensions, seed.isFlipNormals());

  if( not utils::MasterSlave::_slaveMode ){
    globalMesh.addMesh(seed); //add local master mesh to global mesh
  }

  //gather Mesh
  preciceInfo("sendMesh()", "Gather mesh " << seed.getName() );
  if(utils::MasterSlave::_slaveMode || utils::MasterSlave::_masterMode ){
    if(utils::MasterSlave::_rank>0){ //slave
      com::CommunicateMesh(utils::MasterSlave::_communication).sendMesh ( seed, 0 );
    }
    else{ //master
      assertion(utils::MasterSlave::_rank==0);
      assertion(utils::MasterSlave::_size>1);

      int numberOfVertices = 0;
      //vertices of master mesh part do already exist
      foreach ( const mesh::Vertex& vertex, seed.vertices() ){
        seed.getVertexDistribution()[0].push_back(numberOfVertices);
        numberOfVertices++;
      }

      mesh::Mesh slaveMesh("SlaveMesh", _dimensions, seed.isFlipNormals());
      mesh::Mesh& rSlaveMesh = slaveMesh;

      for(int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++){
        rSlaveMesh.clear();
        com::CommunicateMesh(utils::MasterSlave::_communication).receiveMesh ( rSlaveMesh, rankSlave);
        utils::DynVector coord(_dimensions);
        globalMesh.addMesh(rSlaveMesh); //add slave mesh to global mesh

        for(int i = 0; i < rSlaveMesh.vertices().size(); i++){
          seed.getVertexDistribution()[rankSlave].push_back(numberOfVertices);
          numberOfVertices++;
        }
      }
      seed.setGlobalNumberOfVertices(numberOfVertices);
    }
  }

  //send (global) Mesh
  preciceInfo("sendMesh()", "Send global mesh " << seed.getName() );
  if(not utils::MasterSlave::_slaveMode){
    preciceCheck ( globalMesh.vertices().size() > 0,
                   "specializedCreate()", "Participant \"" << _accessorName
                   << "\" provides an invalid (possibly empty) mesh \""
                   << globalMesh.getName() << "\"!" );
    typedef std::map<std::string,m2n::PtrGlobalCommunication>::value_type Pair;
    foreach ( Pair & pair, _receivers ) {
      if ( ! pair.second->isConnected() ) {
        pair.second->acceptConnection ( _providerName, pair.first, 0, 1 );
      }
      com::CommunicateMesh(pair.second->getMasterCommunication()).sendMesh ( globalMesh, 0 );
    }
  }
}

void CommunicatedGeometry:: receiveMesh(
  mesh::Mesh& seed)
{
  preciceInfo("receiveMesh()", "Receive global mesh " << seed.getName() );
  if(not utils::MasterSlave::_slaveMode){
    assertion ( seed.vertices().size() == 0 );
    assertion ( utils::contained(_accessorName, _receivers) );
    m2n::PtrGlobalCommunication com ( _receivers[_accessorName] );
    if ( ! com->isConnected() ) {
      com->requestConnection ( _providerName, _accessorName, 0, 1 );
    }
    com::CommunicateMesh(com->getMasterCommunication()).receiveMesh ( seed, 0 );
  }
  if(utils::MasterSlave::_slaveMode || utils::MasterSlave::_masterMode){
    scatterMesh(seed);
  }
}

void CommunicatedGeometry:: scatterMesh(
    mesh::Mesh& seed)
{
  preciceTrace1 ( "scatterMesh()", utils::MasterSlave::_rank );
  using tarch::la::raw;

  // send bounding boxes from all slaves to master, then send according part of the global mesh back
  preciceInfo("scatterMesh()", "Broadcast global mesh " << seed.getName() );
  std::map<int,std::vector<int> > boundingVertexDistribution;
  if(utils::MasterSlave::_rank>0){ //slave
    mesh::Mesh::BoundingBox bb = mesh::Mesh::BoundingBox (_dimensions,
        std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()));
    mergeBoundingBoxes(bb);
    com::CommunicateMesh(utils::MasterSlave::_communication).sendBoundingBox ( bb, 0);
    com::CommunicateMesh(utils::MasterSlave::_communication).receiveMesh ( seed, 0);
  }
  else{ //master
    assertion(utils::MasterSlave::_rank==0);
    assertion(utils::MasterSlave::_size>1);
    seed.setGlobalNumberOfVertices(seed.vertices().size());
    for(int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++){
      _bb = mesh::Mesh::BoundingBox (_dimensions, std::make_pair(0,0));
      com::CommunicateMesh(utils::MasterSlave::_communication).receiveBoundingBox ( _bb, rankSlave);
      for (int d=0; d<_dimensions; d++){
        if(_safetyGap < _bb[d].second - _bb[d].first) _safetyGap = _bb[d].second - _bb[d].first;
      }
      _safetyGap *= 0.1;
      preciceDebug("From slave " << rankSlave << ", bounding mesh: " << _bb[0].first
          << ", " << _bb[0].second << " and " << _bb[1].first << ", " << _bb[1].second);
      mesh::Mesh slaveMesh("SlaveMesh", _dimensions, seed.isFlipNormals());
      boundingVertexDistribution[rankSlave] = filterMesh(seed, slaveMesh, false);
      com::CommunicateMesh(utils::MasterSlave::_communication).sendMesh ( slaveMesh, rankSlave );
    }
    //now also filter the remaining master mesh
    _bb = mesh::Mesh::BoundingBox (_dimensions,
    std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()));
    mergeBoundingBoxes(_bb);
    for (int d=0; d<_dimensions; d++){
      if(_safetyGap < _bb[d].second - _bb[d].first) _safetyGap = _bb[d].second - _bb[d].first;
    }
    _safetyGap *= 0.1;
    mesh::Mesh filteredMesh("FilteredMesh", _dimensions, seed.isFlipNormals());
    boundingVertexDistribution[0] = filterMesh(seed, filteredMesh, false);
    seed.clear();
    seed.addMesh(filteredMesh);
    preciceDebug("Master mesh after filtering, #vertices " << seed.vertices().size());
  }

  preciceInfo("scatterMesh()", "Compute bounding mappings for mesh " << seed.getName() );
  seed.computeState();
  computeBoundingMappings();

  preciceInfo("scatterMesh()", "Filter mesh " << seed.getName() );
  mesh::Mesh filteredMesh("FilteredMesh", _dimensions, seed.isFlipNormals());
  auto filteredVertexIDs = filterMesh(seed, filteredMesh, true);
  seed.clear();
  seed.addMesh(filteredMesh);
  clearBoundingMappings();

  int numberOfVertices = filteredVertexIDs.size();

  preciceInfo("scatterMesh()", "Gather vertex distribution for mesh " << seed.getName() );
  if(utils::MasterSlave::_rank>0){ //slave
    utils::MasterSlave::_communication->send(numberOfVertices,0);
    if(numberOfVertices!=0){
      utils::MasterSlave::_communication->send(raw(filteredVertexIDs),numberOfVertices,0);
    }
  }
  else{ //master
    assertion(utils::MasterSlave::_rank==0);
    assertion(utils::MasterSlave::_size>1);
    seed.getVertexDistribution()[0] = filteredVertexIDs;
    for(int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++){
      int numberOfVertices = -1;
      utils::MasterSlave::_communication->receive(numberOfVertices,rankSlave);
      std::vector<int> slaveVertexIDs(numberOfVertices,-1);
      if(numberOfVertices!=0){
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

void CommunicatedGeometry:: computeBoundingMappings()
{
  if(_boundingFromMapping.use_count()>0){
    _boundingFromMapping->computeMapping();
  }
  if(_boundingToMapping.use_count()>0){
    _boundingToMapping->computeMapping();
  }
}

void CommunicatedGeometry:: clearBoundingMappings()
{
  if(_boundingFromMapping.use_count()>0){
    _boundingFromMapping->clear();
  }
  if(_boundingToMapping.use_count()>0){
    _boundingToMapping->clear();
  }
}

void CommunicatedGeometry:: mergeBoundingBoxes(mesh::Mesh::BoundingBox& bb){
  if(_boundingFromMapping.use_count()>0){
    auto bb1 = _boundingFromMapping->getOutputMesh()->getBoundingBox();
    for(int d=0; d<_dimensions; d++){
      if(bb[d].first > bb1[d].first) bb[d].first = bb1[d].first;
      if(bb[d].second < bb1[d].second) bb[d].second = bb1[d].second;
    }
  }
  if(_boundingToMapping.use_count()>0){
    auto bb2 = _boundingToMapping->getInputMesh()->getBoundingBox();
    for(int d=0; d<_dimensions; d++){
      if(bb[d].first > bb2[d].first) bb[d].first = bb2[d].first;
      if(bb[d].second < bb2[d].second) bb[d].second = bb2[d].second;
    }
  }
}

bool CommunicatedGeometry:: doesVertexContribute(
    const mesh::Vertex& vertex, bool filterByMapping)
{
  if(filterByMapping){
    //works as easy as this since only read-consistent and write-conservative are allowed
    assertion(_boundingFromMapping.use_count()>0||_boundingToMapping.use_count()>0);
    bool exit = false;
    if(_boundingFromMapping.use_count()>0){
      exit = exit || _boundingFromMapping->doesVertexContribute(vertex.getID());
    }
    if(_boundingToMapping.use_count()>0){
      exit = exit || _boundingToMapping->doesVertexContribute(vertex.getID());
    }
    return exit;
  }
  else{ //filter by bounding box
    for (int d=0; d<_dimensions; d++){
      if(vertex.getCoords()[d] < _bb[d].first - _safetyGap || vertex.getCoords()[d] > _bb[d].second + _safetyGap){
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

  std::vector<int> vertexIDs;
  std::map<int, mesh::Vertex*> vertexMap;
  std::map<int, mesh::Edge*> edgeMap;

  foreach ( const mesh::Vertex& vertex, seed.vertices() ){
    if(doesVertexContribute(vertex, filterByMapping)){
      mesh::Vertex& v = filteredMesh.createVertex(vertex.getCoords());
      vertexIDs.push_back(vertex.getID());
      vertexMap[vertex.getID()] = &v;
    }
  }

  //add all edges formed by contributing vertices
  foreach ( mesh::Edge& edge, seed.edges() ){
    int vertexIndex1 = edge.vertex(0).getID();
    int vertexIndex2 = edge.vertex(1).getID();
    if(vertexMap.find(vertexIndex1) != vertexMap.end() &&
       vertexMap.find(vertexIndex2) != vertexMap.end()){
      mesh::Edge& e = filteredMesh.createEdge(*vertexMap[vertexIndex1], *vertexMap[vertexIndex2]);
      edgeMap[edge.getID()] = &e;
    }
  }

  //add all triangles formed by contributing edges
  if(_dimensions==3){
    foreach (mesh::Triangle& triangle, seed.triangles() ){
      int edgeIndex1 = triangle.edge(0).getID();
      int edgeIndex2 = triangle.edge(1).getID();
      int edgeIndex3 = triangle.edge(2).getID();
      if(edgeMap.find(edgeIndex1) != edgeMap.end() &&
         edgeMap.find(edgeIndex2) != edgeMap.end() &&
         edgeMap.find(edgeIndex3) != edgeMap.end() ){
        filteredMesh.createTriangle(*edgeMap[edgeIndex1],*edgeMap[edgeIndex2],*edgeMap[edgeIndex3]);
      }
    }
  }

  preciceDebug("Filtered mesh. #vertices: " << filteredMesh.vertices().size()
         <<", #edges: " << filteredMesh.edges().size()
         <<", #triangles: " << filteredMesh.triangles().size() << ", rank: " << utils::MasterSlave::_rank);

  return vertexIDs;
}

}} // namespace precice, geometry

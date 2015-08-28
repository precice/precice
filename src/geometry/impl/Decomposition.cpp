// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "Decomposition.hpp"
#include "utils/EventTimings.hpp"
#include "utils/Helpers.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mapping/Mapping.hpp"

using precice::utils::Event;

namespace precice {
namespace geometry {
namespace impl {

tarch::logging::Log Decomposition:: _log ( "precice::geometry::Decomposition" );

Decomposition:: Decomposition
(
  int dimensions, double safetyFactor )
:
  _dimensions (dimensions),
  _boundingFromMapping(),
  _boundingToMapping(),
  _bb(),
  _safetyGap(0.0),
  _safetyFactor(safetyFactor),
  _filterByMapping(false)
{}

void Decomposition:: setBoundingFromMapping(
  mapping::PtrMapping mapping)
{
  _boundingFromMapping = mapping;
}

void Decomposition:: setBoundingToMapping(
  mapping::PtrMapping mapping)
{
  _boundingToMapping = mapping;
}

void Decomposition:: computeBoundingMappings()
{
  if (_boundingFromMapping.use_count() > 0 && _boundingFromMapping->isProjectionMapping()) {
    _boundingFromMapping->computeMapping();
  }
  if (_boundingToMapping.use_count() > 0 && _boundingToMapping->isProjectionMapping()) {
    _boundingToMapping->computeMapping();
  }
}

void Decomposition:: clearBoundingMappings()
{
  if (_boundingFromMapping.use_count() > 0) {
    _boundingFromMapping->clear();
  }
  if (_boundingToMapping.use_count() > 0) {
    _boundingToMapping->clear();
  }
}

std::vector<int> Decomposition:: filterMesh(mesh::Mesh& seed, mesh::Mesh& filteredMesh){
  preciceTrace1 ( "filterMesh()", utils::MasterSlave::_rank );

  preciceDebug("Bounding mesh. #vertices: " << seed.vertices().size()
               <<", #edges: " << seed.edges().size()
               <<", #triangles: " << seed.triangles().size() << ", rank: " << utils::MasterSlave::_rank);

  std::vector<int> vertexPositions;
  std::map<int, mesh::Vertex*> vertexMap;
  std::map<int, mesh::Edge*> edgeMap;
  int vertexCounter = 0;

  for (const mesh::Vertex& vertex : seed.vertices()) {
    if (doesVertexContribute(vertex)){
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

void Decomposition:: mergeBoundingBoxes(mesh::Mesh::BoundingBox& bb){
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


bool Decomposition:: doesVertexContribute(
  const mesh::Vertex& vertex)
{
  if (_filterByMapping) {
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
  else if(_boundingToMapping->isProjectionMapping() && _boundingFromMapping->isProjectionMapping()){ //filter by bounding box
    for (int d=0; d<_dimensions; d++) {
      if (vertex.getCoords()[d] < _bb[d].first - _safetyGap || vertex.getCoords()[d] > _bb[d].second + _safetyGap) {
        return false;
      }
    }
    return true;
  }
  else{ //no filtering here for non-projection mappings (i.e. RBF mappings)
    return true;
  }
}

}}} // namespace precice, geometry, impl

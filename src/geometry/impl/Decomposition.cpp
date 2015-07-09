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
  int dimensions )
:
  _dimensions (dimensions),
  _boundingFromMapping(),
  _boundingToMapping()
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
  if (_boundingFromMapping.use_count() > 0) {
    _boundingFromMapping->computeMapping();
  }
  if (_boundingToMapping.use_count() > 0) {
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

void Decomposition:: feedback(
  mesh::Mesh& seed,
  std::map<int,std::vector<int> >& boundingVertexDistribution,
  std::vector<int>& filteredVertexPositions)
{
  preciceTrace1 ( "feedback()", utils::MasterSlave::_rank );
  preciceInfo("feedback()", "Feedback mesh " << seed.getName() );
  Event e("feedback");

  int numberOfVertices = filteredVertexPositions.size();

  preciceInfo("scatterMesh()", "Gather vertex distribution for mesh " << seed.getName() );
  if (utils::MasterSlave::_slaveMode) {
    utils::MasterSlave::_communication->send(numberOfVertices,0);
    if (numberOfVertices!=0) {
      utils::MasterSlave::_communication->send(filteredVertexPositions.data(),numberOfVertices,0);
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
        utils::MasterSlave::_communication->receive(slaveVertexIDs.data(),numberOfVertices,rankSlave);
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

}}} // namespace precice, geometry, impl

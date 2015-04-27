// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "NearestNeighborMapping.hpp"
#include "query/FindClosestVertex.hpp"

namespace precice {
namespace mapping {

tarch::logging::Log NearestNeighborMapping::
  _log ( "precice::mapping::NearestNeighborMapping" );

NearestNeighborMapping:: NearestNeighborMapping
(
  Constraint constraint,
  int        dimensions)
:
  Mapping(constraint, dimensions),
  _hasComputedMapping(false),
  _vertexIndices()
{
  setInputRequirement(VERTEX);
  setOutputRequirement(VERTEX);
}

void NearestNeighborMapping:: computeMapping()
{
  preciceTrace1("computeMapping()", input()->vertices().size());
  assertion(input().get() != NULL);
  assertion(output().get() != NULL);
  if (getConstraint() == CONSISTENT){
    preciceDebug("Compute consistent mapping");
    size_t verticesSize = output()->vertices().size();
    _vertexIndices.resize(verticesSize);
    const mesh::Mesh::VertexContainer& outputVertices = output()->vertices();
    for ( size_t i=0; i < verticesSize; i++ ){
      const utils::DynVector& coords = outputVertices[i].getCoords();
      query::FindClosestVertex find(coords);
      find(*input());
      assertion(find.hasFound());
      _vertexIndices[i] = find.getClosestVertex().getID();
    }
  }
  else {
    assertion1(getConstraint() == CONSERVATIVE, getConstraint());
    preciceDebug("Compute conservative mapping");
    size_t verticesSize = input()->vertices().size();
    _vertexIndices.resize(verticesSize);
    const mesh::Mesh::VertexContainer& inputVertices = input()->vertices();
    for ( size_t i=0; i < verticesSize; i++ ){
      const utils::DynVector& coords = inputVertices[i].getCoords();
      query::FindClosestVertex find(coords);
      find(*output());
      assertion(find.hasFound());
      _vertexIndices[i] = find.getClosestVertex().getID();
    }
  }
  _hasComputedMapping = true;
}

bool NearestNeighborMapping:: hasComputedMapping()
{
  preciceTrace1("hasComputedMapping()", _hasComputedMapping);
  return _hasComputedMapping;
}

void NearestNeighborMapping:: clear()
{
  preciceTrace("clear()");
  _vertexIndices.clear();
  _hasComputedMapping = false;
}

void NearestNeighborMapping:: map
(
  int inputDataID,
  int outputDataID )
{
  preciceTrace2 ( "map()", inputDataID, outputDataID );
  const utils::DynVector& inputValues = input()->data(inputDataID)->values();
  utils::DynVector& outputValues = output()->data(outputDataID)->values();
  //assign(outputValues) = 0.0;
  int valueDimensions = input()->data(inputDataID)->getDimensions();
  assertion2 ( valueDimensions == output()->data(outputDataID)->getDimensions(),
              valueDimensions, output()->data(outputDataID)->getDimensions() );
  assertion3 ( inputValues.size() / valueDimensions == (int)input()->vertices().size(),
               inputValues.size(), valueDimensions, input()->vertices().size() );
  assertion3 ( outputValues.size() / valueDimensions == (int)output()->vertices().size(),
               outputValues.size(), valueDimensions, output()->vertices().size() );
  if (getConstraint() == CONSISTENT){
    preciceDebug("Map consistent");
    size_t outSize = output()->vertices().size();
    for ( size_t i=0; i < outSize; i++ ){
      int inputIndex = _vertexIndices[i] * valueDimensions;
      for ( int dim=0; dim < valueDimensions; dim++ ){
        outputValues[(i*valueDimensions)+dim] = inputValues[inputIndex+dim];
      }
    }
  }
  else {
    assertion1(getConstraint() == CONSERVATIVE, getConstraint());
    preciceDebug("Map conservative");
    size_t inSize = input()->vertices().size();
    for ( size_t i=0; i < inSize; i++ ){
      int outputIndex = _vertexIndices[i] * valueDimensions;
      for ( int dim=0; dim < valueDimensions; dim++ ){
        outputValues[outputIndex+dim] += inputValues[(i*valueDimensions)+dim];
      }
    }
  }
}

bool NearestNeighborMapping::doesVertexContribute(
  int vertexID)
{
  return utils::contained(vertexID,_vertexIndices);
}

}} // namespace precice, mapping

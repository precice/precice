#include "NearestNeighborMapping.hpp"
#include "query/FindClosestVertex.hpp"
#include "logging/LogMakros.hpp"

namespace precice {
namespace mapping {

logging::Logger NearestNeighborMapping::
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
  ppreciceTrace1("computeMapping()", input()->vertices().size());
  assertion(input().get() != nullptr);
  assertion(output().get() != nullptr);
  if (getConstraint() == CONSISTENT){
    ppreciceDebug("Compute consistent mapping");
    size_t verticesSize = output()->vertices().size();
    _vertexIndices.resize(verticesSize);
    const mesh::Mesh::VertexContainer& outputVertices = output()->vertices();
    for ( size_t i=0; i < verticesSize; i++ ){
      const utils::DynVector& coords = outputVertices[i].getCoords();
      query::FindClosestVertex find(coords); // Search for the output vertex ...
      find(*input()); // ... inside the input mesh
      assertion(find.hasFound());
      _vertexIndices[i] = find.getClosestVertex().getID();
    }
  }
  else {
    assertion1(getConstraint() == CONSERVATIVE, getConstraint());
    ppreciceDebug("Compute conservative mapping");
    size_t verticesSize = input()->vertices().size();
    _vertexIndices.resize(verticesSize);
    const mesh::Mesh::VertexContainer& inputVertices = input()->vertices();
    for ( size_t i=0; i < verticesSize; i++ ){
      const utils::DynVector& coords = inputVertices[i].getCoords();
      query::FindClosestVertex find(coords); // Search for the input vertex ...
      find(*output()); // ... inside the output mesh
      assertion(find.hasFound());
      _vertexIndices[i] = find.getClosestVertex().getID();
    }
  }
  _hasComputedMapping = true;
}

bool NearestNeighborMapping:: hasComputedMapping() const
{
  ppreciceTrace1("hasComputedMapping()", _hasComputedMapping);
  return _hasComputedMapping;
}

void NearestNeighborMapping:: clear()
{
  ppreciceTrace("clear()");
  _vertexIndices.clear();
  _hasComputedMapping = false;
}

void NearestNeighborMapping:: map
(
  int inputDataID,
  int outputDataID )
{
  ppreciceTrace2 ( "map()", inputDataID, outputDataID );
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
    ppreciceDebug("Map consistent");
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
    ppreciceDebug("Map conservative");
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
  int vertexID) const
{
  return utils::contained(vertexID,_vertexIndices);
}

bool NearestNeighborMapping:: isProjectionMapping() const
{
  return true;
}

}} // namespace precice, mapping

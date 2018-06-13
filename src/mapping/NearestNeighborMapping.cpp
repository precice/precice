#include "NearestNeighborMapping.hpp"
#include "query/FindClosestVertex.hpp"
#include "utils/Helpers.hpp"
#include "mesh/RTree.hpp"
#include <Eigen/Core>
#include <boost/function_output_iterator.hpp>

namespace precice {
namespace mapping {

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
  TRACE(input()->vertices().size());
  assertion(input().get() != nullptr);
  assertion(output().get() != nullptr);
  
  if (getConstraint() == CONSISTENT){
    DEBUG("Compute consistent mapping");
    mesh::rtree::PtrRTree rtree = mesh::rtree::getVertexRTree(input());
    size_t verticesSize = output()->vertices().size();
    _vertexIndices.resize(verticesSize);
    const mesh::Mesh::VertexContainer& outputVertices = output()->vertices();
    for ( size_t i=0; i < verticesSize; i++ ) {
        const Eigen::VectorXd& coords = outputVertices[i].getCoords();
        // Search for the output vertex inside the input mesh and add index to _vertexIndices
        rtree->query(boost::geometry::index::nearest(coords, 1),
                     boost::make_function_output_iterator([&](size_t const& val) {
                         _vertexIndices[i] =  input()->vertices()[val].getID();
                       }));
    }
  }
  else {
    assertion(getConstraint() == CONSERVATIVE, getConstraint());
    DEBUG("Compute conservative mapping");
    mesh::rtree::PtrRTree rtree = mesh::rtree::getVertexRTree(output());
    size_t verticesSize = input()->vertices().size();
    _vertexIndices.resize(verticesSize);
    const mesh::Mesh::VertexContainer& inputVertices = input()->vertices();
    for ( size_t i=0; i < verticesSize; i++ ){
      const Eigen::VectorXd& coords = inputVertices[i].getCoords();
      // Search for the input vertex inside the output mesh and add index to _vertexIndices
      rtree->query(boost::geometry::index::nearest(coords, 1),
                   boost::make_function_output_iterator([&](size_t const& val) {
                       _vertexIndices[i] =  output()->vertices()[val].getID();
                     }));
    }
  }
  _hasComputedMapping = true;
}

bool NearestNeighborMapping:: hasComputedMapping() const
{
  TRACE(_hasComputedMapping);
  return _hasComputedMapping;
}

void NearestNeighborMapping:: clear()
{
  TRACE();
  _vertexIndices.clear();
  _hasComputedMapping = false;
}

void NearestNeighborMapping:: map
(
  int inputDataID,
  int outputDataID )
{
  TRACE(inputDataID, outputDataID);
  const Eigen::VectorXd& inputValues = input()->data(inputDataID)->values();
  Eigen::VectorXd& outputValues = output()->data(outputDataID)->values();
  //assign(outputValues) = 0.0;
  int valueDimensions = input()->data(inputDataID)->getDimensions();
  assertion ( valueDimensions == output()->data(outputDataID)->getDimensions(),
              valueDimensions, output()->data(outputDataID)->getDimensions() );
  assertion ( inputValues.size() / valueDimensions == (int)input()->vertices().size(),
               inputValues.size(), valueDimensions, input()->vertices().size() );
  assertion ( outputValues.size() / valueDimensions == (int)output()->vertices().size(),
               outputValues.size(), valueDimensions, output()->vertices().size() );
  if (getConstraint() == CONSISTENT){
    DEBUG("Map consistent");
    size_t outSize = output()->vertices().size();
    for ( size_t i=0; i < outSize; i++ ){
      int inputIndex = _vertexIndices[i] * valueDimensions;
      for ( int dim=0; dim < valueDimensions; dim++ ){
        outputValues((i*valueDimensions)+dim) = inputValues(inputIndex+dim);
      }
    }
  }
  else {
    assertion(getConstraint() == CONSERVATIVE, getConstraint());
    DEBUG("Map conservative");
    size_t inSize = input()->vertices().size();
    for ( size_t i=0; i < inSize; i++ ){
      int outputIndex = _vertexIndices[i] * valueDimensions;
      for ( int dim=0; dim < valueDimensions; dim++ ){
        outputValues(outputIndex+dim) += inputValues((i*valueDimensions)+dim);
      }
    }
  }
}

void NearestNeighborMapping::tagMeshFirstRound()
{
  TRACE();

  computeMapping();

  if (getConstraint() == CONSISTENT){
    for(mesh::Vertex& v : input()->vertices()){
      if(utils::contained(v.getID(),_vertexIndices)) v.tag();
    }
  }
  else {
    assertion(getConstraint() == CONSERVATIVE, getConstraint());
    for(mesh::Vertex& v : output()->vertices()){
      if(utils::contained(v.getID(),_vertexIndices)) v.tag();
    }
  }

  clear();
}

void NearestNeighborMapping::tagMeshSecondRound()
{
  TRACE();
  // for NN mapping no operation needed here
}

}} // namespace precice, mapping

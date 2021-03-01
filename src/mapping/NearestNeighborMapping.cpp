#include "NearestNeighborMapping.hpp"

#include <Eigen/Core>
#include <boost/container/flat_set.hpp>
#include <functional>
#include <memory>
#include "logging/LogMacros.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Vertex.hpp"
#include "query/Index.hpp"
#include "utils/Event.hpp"
#include "utils/Statistics.hpp"
#include "utils/assertion.hpp"

namespace precice {
extern bool syncMode;

namespace mapping {

NearestNeighborMapping::NearestNeighborMapping(
    Constraint constraint,
    int        dimensions)
    : Mapping(constraint, dimensions)
{
  if (hasConstraint(SCALEDCONSISTENT)) {
    setInputRequirement(Mapping::MeshRequirement::FULL);
    setOutputRequirement(Mapping::MeshRequirement::FULL);
  } else {
    setInputRequirement(Mapping::MeshRequirement::VERTEX);
    setOutputRequirement(Mapping::MeshRequirement::VERTEX);
  }
}

void NearestNeighborMapping::computeMapping()
{
  PRECICE_TRACE(input()->vertices().size());

  PRECICE_ASSERT(input().get() != nullptr);
  PRECICE_ASSERT(output().get() != nullptr);

  const std::string     baseEvent = "map.nn.computeMapping.From" + input()->getName() + "To" + output()->getName();
  precice::utils::Event e(baseEvent, precice::syncMode);

  if (hasConstraint(CONSERVATIVE)) {
    PRECICE_DEBUG("Compute conservative mapping");
    precice::utils::Event e2(baseEvent + ".getIndexOnVertices", precice::syncMode);
    query::Index          indexTree(output());
    e2.stop();
    size_t verticesSize = input()->vertices().size();
    _vertexIndices.resize(verticesSize);
    utils::statistics::DistanceAccumulator distanceStatistics;
    const mesh::Mesh::VertexContainer &    inputVertices = input()->vertices();
    // Search for the output vertex inside the input mesh and add index to _vertexIndices
    for (size_t i = 0; i < verticesSize; i++) {
      auto matchedVertex = indexTree.getClosestVertex(inputVertices[i]);
      _vertexIndices[i]  = matchedVertex.index;
      distanceStatistics(matchedVertex.distance);
    }
    if (distanceStatistics.empty()) {
      PRECICE_INFO("Mapping distance not available due to empty partition.");
    } else {
      PRECICE_INFO("Mapping distance " << distanceStatistics);
    }
  } else {
    PRECICE_DEBUG("Compute consistent mapping");
    precice::utils::Event e2(baseEvent + ".getIndexOnVertices", precice::syncMode);
    query::Index          indexTree(input());
    e2.stop();
    size_t verticesSize = output()->vertices().size();
    _vertexIndices.resize(verticesSize);
    utils::statistics::DistanceAccumulator distanceStatistics;
    const mesh::Mesh::VertexContainer &    outputVertices = output()->vertices();
    for (size_t i = 0; i < verticesSize; i++) {
      auto matchedVertex = indexTree.getClosestVertex(outputVertices[i]);
      _vertexIndices[i]  = matchedVertex.index;
      distanceStatistics(matchedVertex.distance);
    }
    if (distanceStatistics.empty()) {
      PRECICE_INFO("Mapping distance not available due to empty partition.");
    } else {
      PRECICE_INFO("Mapping distance " << distanceStatistics);
    }
  }
  _hasComputedMapping = true;
}

bool NearestNeighborMapping::hasComputedMapping() const
{
  PRECICE_TRACE(_hasComputedMapping);
  return _hasComputedMapping;
}

void NearestNeighborMapping::clear()
{
  PRECICE_TRACE();
  _vertexIndices.clear();
  _hasComputedMapping = false;
  if (getConstraint() == CONSISTENT) {
    query::clearCache(input()->getID());
  } else {
    query::clearCache(output()->getID());
  }
}

void NearestNeighborMapping::map(
    int inputDataID,
    int outputDataID)
{
  PRECICE_TRACE(inputDataID, outputDataID);

  precice::utils::Event e("map.nn.mapData.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);

  const Eigen::VectorXd &inputValues  = input()->data(inputDataID)->values();
  Eigen::VectorXd &      outputValues = output()->data(outputDataID)->values();
  //assign(outputValues) = 0.0;
  int valueDimensions = input()->data(inputDataID)->getDimensions();
  PRECICE_ASSERT(valueDimensions == output()->data(outputDataID)->getDimensions(),
                 valueDimensions, output()->data(outputDataID)->getDimensions());
  PRECICE_ASSERT(inputValues.size() / valueDimensions == (int) input()->vertices().size(),
                 inputValues.size(), valueDimensions, input()->vertices().size());
  PRECICE_ASSERT(outputValues.size() / valueDimensions == (int) output()->vertices().size(),
                 outputValues.size(), valueDimensions, output()->vertices().size());
  if (hasConstraint(CONSERVATIVE)) {
    PRECICE_DEBUG("Map conservative");
    size_t const inSize = input()->vertices().size();
    for (size_t i = 0; i < inSize; i++) {
      int const outputIndex = _vertexIndices[i] * valueDimensions;
      for (int dim = 0; dim < valueDimensions; dim++) {
        outputValues(outputIndex + dim) += inputValues((i * valueDimensions) + dim);
      }
    }
  } else {
    PRECICE_DEBUG("Map consistent");
    size_t const outSize = output()->vertices().size();
    for (size_t i = 0; i < outSize; i++) {
      int inputIndex = _vertexIndices[i] * valueDimensions;
      for (int dim = 0; dim < valueDimensions; dim++) {
        outputValues((i * valueDimensions) + dim) = inputValues(inputIndex + dim);
      }
    }
    if (hasConstraint(SCALEDCONSISTENT)) {
      scaleConsistentMapping(inputDataID, outputDataID);
    }
  }
}

void NearestNeighborMapping::tagMeshFirstRound()
{
  PRECICE_TRACE();
  precice::utils::Event e("map.nn.tagMeshFirstRound.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);

  computeMapping();

  // Lookup table of all indices used in the mapping
  const boost::container::flat_set<int> indexSet(_vertexIndices.begin(), _vertexIndices.end());

  if (hasConstraint(CONSERVATIVE)) {
    PRECICE_ASSERT(getConstraint() == CONSERVATIVE, getConstraint());
    for (mesh::Vertex &v : output()->vertices()) {
      if (indexSet.count(v.getID()) != 0)
        v.tag();
    }
  } else {
    for (mesh::Vertex &v : input()->vertices()) {
      if (indexSet.count(v.getID()) != 0)
        v.tag();
    }
  }

  clear();
}

void NearestNeighborMapping::tagMeshSecondRound()
{
  PRECICE_TRACE();
  // for NN mapping no operation needed here
}

} // namespace mapping
} // namespace precice

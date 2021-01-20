#include "NearestNeighborMapping.hpp"

#include <Eigen/Core>
#include <boost/container/flat_set.hpp>
#include <functional>
#include <memory>

#include <boost/version.hpp>
#if BOOST_VERSION < 106600
#include <boost/function_output_iterator.hpp>
#else
#include <boost/iterator/function_output_iterator.hpp>
#endif

#include "logging/LogMacros.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/RTree.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Vertex.hpp"
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
  setInputRequirement(Mapping::MeshRequirement::VERTEX);
  setOutputRequirement(Mapping::MeshRequirement::VERTEX);
}

void NearestNeighborMapping::computeMapping()
{
  PRECICE_TRACE(input()->vertices().size());

  PRECICE_ASSERT(input().get() != nullptr);
  PRECICE_ASSERT(output().get() != nullptr);

  const std::string     baseEvent = "map.nn.computeMapping.From" + input()->getName() + "To" + output()->getName();
  precice::utils::Event e(baseEvent, precice::syncMode);

  if (getConstraint() == CONSISTENT) {
    PRECICE_DEBUG("Compute consistent mapping");
    precice::utils::Event e2(baseEvent + ".getIndexOnVertices", precice::syncMode);
    auto                  rtree = mesh::rtree::getVertexRTree(input());
    e2.stop();
    size_t verticesSize = output()->vertices().size();
    _vertexIndices.resize(verticesSize);
    utils::statistics::DistanceAccumulator distanceStatistics;
    const mesh::Mesh::VertexContainer &    outputVertices = output()->vertices();
    for (size_t i = 0; i < verticesSize; i++) {
      const Eigen::VectorXd &coords = outputVertices[i].getCoords();
      // Search for the output vertex inside the input mesh and add index to _vertexIndices
      rtree->query(boost::geometry::index::nearest(coords, 1),
                   boost::make_function_output_iterator([&](size_t const &val) {
                     const auto &match = input()->vertices()[val];
                     _vertexIndices[i] = match.getID();
                     distanceStatistics(bg::distance(match, coords));
                   }));
    }
    if (distanceStatistics.empty()) {
      PRECICE_INFO("Mapping distance not available due to empty partition.");
    } else {
      PRECICE_INFO("Mapping distance " << distanceStatistics);
    }
  } else {
    PRECICE_ASSERT(getConstraint() == CONSERVATIVE, getConstraint());
    PRECICE_DEBUG("Compute conservative mapping");
    precice::utils::Event e2(baseEvent + ".getIndexOnVertices", precice::syncMode);
    auto                  rtree = mesh::rtree::getVertexRTree(output());
    e2.stop();
    size_t verticesSize = input()->vertices().size();
    _vertexIndices.resize(verticesSize);
    utils::statistics::DistanceAccumulator distanceStatistics;
    const mesh::Mesh::VertexContainer &    inputVertices = input()->vertices();
    for (size_t i = 0; i < verticesSize; i++) {
      const Eigen::VectorXd &coords = inputVertices[i].getCoords();
      // Search for the input vertex inside the output mesh and add index to _vertexIndices
      rtree->query(boost::geometry::index::nearest(coords, 1),
                   boost::make_function_output_iterator([&](size_t const &val) {
                     const auto &match = output()->vertices()[val];
                     _vertexIndices[i] = match.getID();
                     distanceStatistics(bg::distance(match, coords));
                   }));
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
    mesh::rtree::clear(*input());
  } else {
    mesh::rtree::clear(*output());
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
  if (getConstraint() == CONSISTENT) {
    PRECICE_DEBUG("Map consistent");
    size_t const outSize = output()->vertices().size();
    for (size_t i = 0; i < outSize; i++) {
      int inputIndex = _vertexIndices[i] * valueDimensions;
      for (int dim = 0; dim < valueDimensions; dim++) {
        outputValues((i * valueDimensions) + dim) = inputValues(inputIndex + dim);
      }
    }
  } else {
    PRECICE_ASSERT(getConstraint() == CONSERVATIVE, getConstraint());
    PRECICE_DEBUG("Map conservative");
    size_t const inSize = input()->vertices().size();
    for (size_t i = 0; i < inSize; i++) {
      int const outputIndex = _vertexIndices[i] * valueDimensions;
      for (int dim = 0; dim < valueDimensions; dim++) {
        outputValues(outputIndex + dim) += inputValues((i * valueDimensions) + dim);
      }
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

  if (getConstraint() == CONSISTENT) {
    for (mesh::Vertex &v : input()->vertices()) {
      if (indexSet.count(v.getID()) != 0)
        v.tag();
    }
  } else {
    PRECICE_ASSERT(getConstraint() == CONSERVATIVE, getConstraint());
    for (mesh::Vertex &v : output()->vertices()) {
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

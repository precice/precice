
#include "NearestNeighborBaseMapping.hpp"

#include <boost/container/flat_set.hpp>
#include <functional>
#include <iostream>
#include "logging/LogMacros.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Vertex.hpp"
#include "query/Index.hpp"
#include "utils/Event.hpp"
#include "utils/Parallel.hpp"
#include "utils/Statistics.hpp"
#include "utils/assertion.hpp"

namespace precice {
extern bool syncMode;

namespace mapping {

NearestNeighborBaseMapping::NearestNeighborBaseMapping(
    Constraint  constraint,
    int         dimensions,
    bool        hasGradient,
    std::string mappingName,
    std::string mappingNameShort)
    : Mapping(constraint, dimensions),
      _hasGradient(hasGradient),
      MAPPING_NAME(mappingName),
      MAPPING_NAME_SHORT(mappingNameShort)
{
}

bool NearestNeighborBaseMapping::hasGradient()
{
  return _hasGradient;
}

void NearestNeighborBaseMapping::computeMapping()
{
  PRECICE_TRACE(input()->vertices().size());

  PRECICE_ASSERT(input().get() != nullptr);
  PRECICE_ASSERT(output().get() != nullptr);

  const std::string     baseEvent = "map." + MAPPING_NAME_SHORT + ".computeMapping.From" + input()->getName() + "To" + output()->getName();
  precice::utils::Event e(baseEvent, precice::syncMode);

  // Setup Direction of Mapping
  mesh::PtrMesh origins, searchSpace;
  if (hasConstraint(CONSERVATIVE)) {
    PRECICE_DEBUG("Compute conservative mapping");
    origins     = input();
    searchSpace = output();
  } else {
    PRECICE_DEBUG((hasConstraint(CONSISTENT) ? "Compute consistent mapping" : "Compute scaled-consistent mapping"));
    origins     = output();
    searchSpace = input();
  }

  // For log ouputs
  precice::utils::Event e2(baseEvent + ".getIndexOnVertices", precice::syncMode);
  query::Index          indexTree(searchSpace);
  e2.stop();

  // Set up of output arrays
  const size_t verticesSize   = origins->vertices().size();
  const auto & sourceVertices = origins->vertices();
  _vertexIndices.resize(verticesSize);
  if (hasGradient()) {
    _distancesMatched.resize(verticesSize);
  }

  // Needed for error studies
  utils::statistics::DistanceAccumulator distanceStatistics;

  for (size_t i = 0; i < verticesSize; ++i) {

    const auto &matchedVertex = indexTree.getClosestVertex(sourceVertices[i].getCoords());

    // Match the difference vector between the source vector and the matched one (relevant for gradient)
    if (hasGradient()) {

      const auto &matchedVertexCoords = searchSpace.get()->vertices()[matchedVertex.index].getCoords();

      // We calculate the distances uniformly for both mapping constraints as the difference (input - output)
      if (hasConstraint(CONSERVATIVE)) {
        _distancesMatched[i] = sourceVertices[i].getCoords() - matchedVertexCoords;
      } else {
        _distancesMatched[i] = matchedVertexCoords - sourceVertices[i].getCoords();
      }
    }

    _vertexIndices[i] = matchedVertex.index;
    distanceStatistics(matchedVertex.distance);
  }

  // This is the distance object between the coordinates of the vertices and its match in the mesh.
  // It contains the min, max, average and count of accumulated values.
  if (distanceStatistics.empty()) {
    PRECICE_INFO("Mapping distance not available due to empty partition.");
  } else {
    PRECICE_INFO("Mapping distance {}", distanceStatistics);
  }

  _hasComputedMapping = true;
}

bool NearestNeighborBaseMapping::hasComputedMapping() const
{
  PRECICE_TRACE(_hasComputedMapping);
  return _hasComputedMapping;
}

void NearestNeighborBaseMapping::clear()
{
  PRECICE_TRACE();
  _vertexIndices.clear();
  _hasComputedMapping = false;

  if (hasGradient())
    _distancesMatched.clear();

  if (getConstraint() == CONSISTENT) {
    query::clearCache(input()->getID());
  } else {
    query::clearCache(output()->getID());
  }
}

void NearestNeighborBaseMapping::tagMeshFirstRound()
{
  PRECICE_TRACE();
  precice::utils::Event e("map." + MAPPING_NAME_SHORT + ".tagMeshFirstRound.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);

  computeMapping();

  // Lookup table of all indices used in the mapping
  const boost::container::flat_set<int> indexSet(_vertexIndices.begin(), _vertexIndices.end());

  // Get the source mesh depending on the constraint
  const mesh::PtrMesh &source = hasConstraint(CONSERVATIVE) ? output() : input();

  // Tag all vertices used in the mapping
  for (mesh::Vertex &v : source->vertices()) {
    if (indexSet.count(v.getID()) != 0) {
      v.tag();
    }
  }

  clear();
}

void NearestNeighborBaseMapping::tagMeshSecondRound()
{
  PRECICE_TRACE();
  // for NN mapping no operation needed here
}

} // namespace mapping
} // namespace precice

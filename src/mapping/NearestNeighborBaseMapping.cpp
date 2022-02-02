
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
    bool        requireGradient,
    std::string mappingName,
    std::string mappingNameShort)
    : Mapping(constraint, dimensions, requireGradient),
      mappingName(mappingName),
      mappingNameShort(mappingNameShort),
      _requireGradient(requireGradient)
{
}

bool NearestNeighborBaseMapping::requireGradient()
{
  return _requireGradient;
}

void NearestNeighborBaseMapping::computeMapping()
{
  PRECICE_TRACE(input()->vertices().size());

  PRECICE_ASSERT(input().get() != nullptr);
  PRECICE_ASSERT(output().get() != nullptr);

  const std::string     baseEvent = "map." + mappingNameShort + ".computeMapping.From" + input()->getName() + "To" + output()->getName();
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

  // Needed for error calculations
  utils::statistics::DistanceAccumulator distanceStatistics;

  for (size_t i = 0; i < verticesSize; ++i) {
    const auto &matchedVertex = indexTree.getClosestVertex(sourceVertices[i].getCoords());
    _vertexIndices[i]         = matchedVertex.index;
    distanceStatistics(matchedVertex.distance);
  }

  // For gradient mapping, the calculation of offsets between source and matched vertex necessary
  onMappingComputed(origins, searchSpace);

  // This is the distance object between the coordinates of the vertices and its match in the mesh.
  // This prints min, max, average and count of the distances.
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

  if (requireGradient())
    _offsetsMatched.clear();

  if (getConstraint() == CONSISTENT) {
    query::clearCache(input()->getID());
  } else {
    query::clearCache(output()->getID());
  }
}

void NearestNeighborBaseMapping::onMappingComputed(mesh::PtrMesh origins, mesh::PtrMesh searchSpace)
{
  // Does nothing by default
}

void NearestNeighborBaseMapping::tagMeshFirstRound()
{
  PRECICE_TRACE();
  precice::utils::Event e("map." + mappingNameShort + ".tagMeshFirstRound.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);

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

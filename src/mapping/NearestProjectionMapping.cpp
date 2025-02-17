#include "NearestProjectionMapping.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include <ostream>
#include <utility>

#include "logging/LogMacros.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/Polation.hpp"
#include "math/differences.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Vertex.hpp"
#include "profiling/Event.hpp"
#include "query/Index.hpp"
#include "utils/IntraComm.hpp"
#include "utils/Statistics.hpp"
#include "utils/assertion.hpp"

namespace precice::mapping {

NearestProjectionMapping::NearestProjectionMapping(
    Constraint constraint,
    int        dimensions)
    : BarycentricBaseMapping(constraint, dimensions)
{
  if (constraint == CONSISTENT) {
    setInputRequirement(Mapping::MeshRequirement::FULL);
    setOutputRequirement(Mapping::MeshRequirement::VERTEX);
  } else if (constraint == CONSERVATIVE) {
    setInputRequirement(Mapping::MeshRequirement::VERTEX);
    setOutputRequirement(Mapping::MeshRequirement::FULL);
  } else {
    PRECICE_ASSERT(isScaledConsistent());
    setInputRequirement(Mapping::MeshRequirement::FULL);
    setOutputRequirement(Mapping::MeshRequirement::FULL);
  }

  PRECICE_CHECK(constraint != SCALED_CONSISTENT_VOLUME, "Nearest-projection can't be used with volume version of the scaled-consistent mapping. Use scaled-consistent instead.");
}

void NearestProjectionMapping::computeMapping()
{
  PRECICE_TRACE(input()->nVertices(), output()->nVertices());
  const std::string         baseEvent = "map.np.computeMapping.From" + input()->getName() + "To" + output()->getName();
  precice::profiling::Event e(baseEvent, profiling::Synchronize);

  // Setup Direction of Mapping
  mesh::PtrMesh origins, searchSpace;
  if (hasConstraint(CONSERVATIVE)) {
    PRECICE_DEBUG("Compute conservative mapping");
    origins     = input();
    searchSpace = output();
  } else {
    PRECICE_DEBUG("Compute consistent mapping");
    origins     = output();
    searchSpace = input();
  }

  const auto &fVertices = origins->vertices();

  if (getDimensions() == 2) {
    PRECICE_WARN_IF(!fVertices.empty() && searchSpace->edges().empty(),
                    "2D Mesh \"{}\" does not contain edges. "
                    "Nearest projection mapping falls back to nearest neighbor mapping.",
                    searchSpace->getName());
  } else {
    PRECICE_WARN_IF(!fVertices.empty() && searchSpace->triangles().empty(),
                    "3D Mesh \"{}\" does not contain triangles. "
                    "Nearest projection mapping will map to primitives of lower dimension.",
                    searchSpace->getName());
  }

  // Amount of nearest elements to fetch for detailed comparison.
  // This safety margin results in a candidate set which forms the base for the
  // local nearest projection and counters the loss of detail due to bounding box generation.
  // @TODO Add a configuration option for this factor
  constexpr int nnearest = 4;

  utils::statistics::DistanceAccumulator distanceStatistics;
  std::size_t                            toTriangles{0}, toEdges{0}, toVertices{0};

  _interpolations.clear();
  _interpolations.reserve(fVertices.size());

  auto &index = searchSpace->index();
  for (const auto &fVertex : fVertices) {
    // Nearest projection element is edge for 2d if exists, if not, it is the nearest vertex
    // Nearest projection element is triangle for 3d if exists, if not the edge and at the worst case it is the nearest vertex
    auto match = index.findNearestProjection(fVertex.getCoords(), nnearest);
    distanceStatistics(match.polation.distance());
    switch (match.polation.nElements()) {
    case 1:
      ++toVertices;
      break;
    case 2:
      ++toEdges;
      break;
    case 3:
      ++toTriangles;
      break;
    default:
      PRECICE_UNREACHABLE("");
    }

    _interpolations.push_back(std::move(match.polation));
  }

  if (distanceStatistics.empty()) {
    PRECICE_INFO("Mapping distance not available due to empty partition.");
  } else {
    PRECICE_INFO("Mapping distance {}", distanceStatistics);
    PRECICE_INFO("Nearest-projections are {} triangles, {} edges, and {} vertices", toTriangles, toEdges, toVertices);
  }

  _hasComputedMapping = true;
}

std::string NearestProjectionMapping::getName() const
{
  return "nearest-projection";
}

} // namespace precice::mapping

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
#include "query/Index.hpp"
#include "utils/Event.hpp"
#include "utils/Statistics.hpp"
#include "utils/assertion.hpp"

namespace precice {
extern bool syncMode;

namespace mapping {

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
    PRECICE_ASSERT(isScaledConsistent(), constraint);
    setInputRequirement(Mapping::MeshRequirement::FULL);
    setOutputRequirement(Mapping::MeshRequirement::FULL);
  }

  PRECICE_CHECK(constraint != SCALED_CONSISTENT_VOLUME,
                ::precice::MappingError, "Nearest-projection can't be used with volume version of the scaled-consistent mapping. Use scaled-consistent instead.");
}

void NearestProjectionMapping::computeMapping()
{
  PRECICE_TRACE(input()->vertices().size(), output()->vertices().size());
  const std::string     baseEvent = "map.np.computeMapping.From" + input()->getName() + "To" + output()->getName();
  precice::utils::Event e(baseEvent, precice::syncMode);

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
    if (!fVertices.empty() && searchSpace->edges().empty()) {
      PRECICE_WARN("2D Mesh \"{}\" does not contain edges. "
                   "Nearest projection mapping falls back to nearest neighbor mapping.",
                   searchSpace->getName());
    }
  } else {
    if (!fVertices.empty() && searchSpace->triangles().empty()) {
      PRECICE_WARN("3D Mesh \"{}\" does not contain triangles. "
                   "Nearest projection mapping will map to primitives of lower dimension.",
                   searchSpace->getName());
    }
  }

  // Amount of nearest elements to fetch for detailed comparison.
  // This safety margin results in a candidate set which forms the base for the
  // local nearest projection and counters the loss of detail due to bounding box generation.
  // @TODO Add a configuration option for this factor
  constexpr int nnearest = 4;

  utils::statistics::DistanceAccumulator distanceStatistics;

  _interpolations.clear();
  _interpolations.reserve(fVertices.size());

  auto &index = searchSpace->index();
  for (const auto &fVertex : fVertices) {
    // Nearest projection element is edge for 2d if exists, if not, it is the nearest vertex
    // Nearest projection element is triangle for 3d if exists, if not the edge and at the worst case it is the nearest vertex
    auto match = index.findNearestProjection(fVertex.getCoords(), nnearest);
    distanceStatistics(match.polation.distance());
    _interpolations.push_back(std::move(match.polation));
  }

  if (distanceStatistics.empty()) {
    PRECICE_INFO("Mapping distance not available due to empty partition.");
  } else {
    PRECICE_INFO("Mapping distance {}", distanceStatistics);
  }

  _hasComputedMapping = true;
}

} // namespace mapping
} // namespace precice

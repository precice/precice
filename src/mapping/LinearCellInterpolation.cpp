#include "LinearCellInterpolation.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <deque>
#include <memory>
#include <ostream>
#include <unordered_set>
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

LinearCellInterpolation::LinearCellInterpolation(
    Constraint constraint,
    int        dimensions)
    : BarycentricBaseMapping(constraint, dimensions)
{
  PRECICE_CHECK(getDimensions() == 2, "Volume mapping not available in 3D.");
  if (constraint == CONSISTENT) {
    setInputRequirement(Mapping::MeshRequirement::FULL);
    setOutputRequirement(Mapping::MeshRequirement::VERTEX);
  } else if (constraint == CONSERVATIVE) {
    setInputRequirement(Mapping::MeshRequirement::VERTEX);
    setOutputRequirement(Mapping::MeshRequirement::FULL);
  } else {
    PRECICE_ASSERT(constraint == SCALEDCONSISTENT, constraint);
    setInputRequirement(Mapping::MeshRequirement::FULL);
    setOutputRequirement(Mapping::MeshRequirement::FULL);
  }

  PRECICE_CHECK(constraint != SCALEDCONSISTENT, "Volume mapping doesn't support scaled-consistent mappings.");
}

void LinearCellInterpolation::computeMapping()
{
  PRECICE_TRACE(input()->vertices().size(), output()->vertices().size());
  const std::string     baseEvent = "map.vci.computeMapping.From" + input()->getName() + "To" + output()->getName();
  precice::utils::Event e(baseEvent, precice::syncMode);
  PRECICE_ASSERT(getDimensions() == 2, "Volume mapping not available in 3D.");

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
    if (!fVertices.empty() && searchSpace->triangles().empty()) {
      PRECICE_WARN("2D Mesh \"{}\" does not contain triangles. "
                   "Linear cell interpolation falls back to nearest projection mapping.",
                   searchSpace->getName());
    }
  } else {
    if (!fVertices.empty() && searchSpace->triangles().empty()) {
      PRECICE_WARN("3D Mesh \"{}\" does not contain tetrahedra. "
                   "Linear cell interpolation falls back to nearest projection mapping.",
                   searchSpace->getName());
    }
  }

  // Amount of nearest elements to fetch for detailed comparison.
  // This safety margin results in a candidate set which forms the base for the
  // local nearest projection and counters the loss of detail due to bounding box generation.
  // @TODO Add a configuration option for this factor
  constexpr int nnearest = 4;

  auto                                  &index = searchSpace->index();
  utils::statistics::DistanceAccumulator fallbackStatistics;

  _interpolations.clear();
  _interpolations.reserve(fVertices.size());

  for (const auto &fVertex : fVertices) {
    // Find tetrahedra (3D) or triangle (2D) or fall-back on NP
    auto match = index.findCellInterpolation(fVertex.getCoords(), nnearest);
    _interpolations.push_back(std::move(match.polation));
    if (match.distance != 0.0)
    {
      // Only push when fall-back occurs, so the number of entries is the number of vertices outside the domain
      fallbackStatistics(match.distance);
    }
  }

  if (!fallbackStatistics.empty()) {
    PRECICE_WARN("Some points are outisde of the domain defined by connectivity. Fall-back on Nearest-Projection occured."
                " Statistics: {}", fallbackStatistics);
    
  }

  _hasComputedMapping = true;
}


} // namespace mapping
} // namespace precice

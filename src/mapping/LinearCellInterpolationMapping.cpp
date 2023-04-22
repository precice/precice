#include "LinearCellInterpolationMapping.hpp"
#include "logging/LogMacros.hpp"
#include "query/Index.hpp"
#include "utils/Event.hpp"
#include "utils/Statistics.hpp"
#include "utils/assertion.hpp"

namespace precice {
extern bool syncMode;

namespace mapping {

LinearCellInterpolationMapping::LinearCellInterpolationMapping(
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

  PRECICE_CHECK(constraint != SCALED_CONSISTENT_SURFACE, "Volume mapping doesn't support scaled-consistent-surface mappings. Use \"scaled-consistent-volume\" instead.");
}

void LinearCellInterpolationMapping::computeMapping()
{
  PRECICE_TRACE(input()->vertices().size(), output()->vertices().size());
  const std::string     baseEvent = "map.vci.computeMapping.From" + input()->getName() + "To" + output()->getName();
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

  const auto &fVertices       = origins->vertices();
  bool        hasConnectivity = true;

  if (getDimensions() == 2) {
    if (!fVertices.empty() && searchSpace->triangles().empty()) {
      const bool hasEdges = searchSpace->hasEdges();
      PRECICE_WARN("2D Mesh \"{}\" does not contain triangles{} "
                   "Linear cell interpolation falls back to nearest-{} mapping.",
                   searchSpace->getName(), hasEdges ? "." : " or edges.", hasEdges ? "projection" : "neighbor");
      hasConnectivity = false;
    }
  } else {
    if (!fVertices.empty() && searchSpace->tetrahedra().empty()) {
      const bool hasTriangles = searchSpace->hasTriangles();
      const bool hasEdges     = searchSpace->hasEdges();
      PRECICE_WARN("3D Mesh \"{}\" does not contain tetrahedra{}{} "
                   "Linear cell interpolation falls back to nearest-{} mapping.",
                   searchSpace->getName(), hasEdges ? "" : ", edges", hasTriangles ? "." : " or triangles.", (hasTriangles || hasEdges) ? "projection" : "neighbor");
      hasConnectivity = false;
    }
  }

  // Amount of nearest elements to fetch for detailed comparison.
  // This safety margin results in a candidate set which forms the base for the
  // local nearest projection and counters the loss of detail due to bounding box generation.
  // @TODO Add a configuration option for this factor
  constexpr int nnearest = 4;

  auto &                                 index = searchSpace->index();
  utils::statistics::DistanceAccumulator fallbackStatistics;

  _interpolations.clear();
  _interpolations.reserve(fVertices.size());

  for (const auto &fVertex : fVertices) {
    // Find tetrahedra (3D) or triangle (2D) or fall-back on NP
    auto match    = index.findCellOrProjection(fVertex.getCoords(), nnearest);
    auto distance = match.polation.distance();
    _interpolations.push_back(std::move(match.polation));
    if (!math::equals(distance, 0.0)) {
      // Only push when fall-back occurs, so the number of entries is the number of vertices outside the domain
      fallbackStatistics(distance);
    }
  }

  if (!fallbackStatistics.empty()) {
    if (hasConnectivity) {
      // We have the connectivity, but some fallbacks occurred
      PRECICE_INFO(
          "Linear Cell Interpolation is used, but some points from {} don't lie in the domain defined by the {}. "
          "These points have been projected on the domain boundary. This could come from non-matching discrete geometries or erroneous connectivity information."
          "If distances seem too large, please check your mesh. "
          "The fallback statistics are: {} ",
          searchSpace->getName(), getDimensions() == 2 ? "triangles" : "tetrahedra",
          fallbackStatistics);
    } else {
      // Fallback and no connectivity provided
      PRECICE_INFO("Fallback mapping statistics: {}", fallbackStatistics);
    }
  } else {
    if (!hasConnectivity) {
      // Not all connectivity provided, but no fallback applied
      PRECICE_INFO("All vertices are inside cells, no fallback required");
    } else {
      // No fallback and we have connectivity
      PRECICE_ASSERT(hasConnectivity)
      PRECICE_INFO("Successfully computed linear-cell-interpolation mapping.");
    }
  }

  _hasComputedMapping = true;
}

std::string LinearCellInterpolationMapping::getName() const
{
  return "linear-cell interpolation";
}
} // namespace mapping
} // namespace precice

#include "VolumeCellInterpolation.hpp"

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

VolumeCellInterpolation::VolumeCellInterpolation(
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
    PRECICE_ASSERT(constraint == SCALEDCONSISTENT, constraint);
    setInputRequirement(Mapping::MeshRequirement::FULL);
    setOutputRequirement(Mapping::MeshRequirement::FULL);
  }

  PRECICE_ASSERT(constraint != SCALEDCONSISTENT, "Volume mapping doesn't support scaled-consistent mappings.");
}

void VolumeCellInterpolation::computeMapping()
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
                   "Volume cell interpolation falls back to nearest neighbor mapping.",
                   searchSpace->getName());
    }
  } else {
    if (!fVertices.empty() && searchSpace->triangles().empty()) {
      PRECICE_WARN("3D Mesh \"{}\" does not contain tetrahedrons. "
                   "Volume cell interpolation will map to primitives of lower dimension.",
                   searchSpace->getName());
    }
  }

  // Amount of nearest elements to fetch for detailed comparison.
  // This safety margin results in a candidate set which forms the base for the
  // local nearest projection and counters the loss of detail due to bounding box generation.
  // @TODO Add a configuration option for this factor
  constexpr int nnearest = 4;

  auto &index = searchSpace->index();
  utils::statistics::DistanceAccumulator distanceStatistics;

  _interpolations.clear();
  _interpolations.reserve(fVertices.size());

  for (const auto &fVertex : fVertices) {
    // Find tetrahedra (3D) or triangle (2D) or fall-back on a vertex
    auto match = index.findNearestVolume(fVertex.getCoords(), nnearest);
    _interpolations.push_back(std::move(match.polation));
    distanceStatistics(match.distance);
  }

  if (distanceStatistics.empty()) {
    PRECICE_INFO("Mapping distance not available due to empty partition.");
  } else {
    PRECICE_INFO("Mapping distance {}", distanceStatistics);
  }

  _hasComputedMapping = true;
}

void VolumeCellInterpolation::tagMeshFirstRound()
{
  PRECICE_TRACE();
  precice::utils::Event e("map.vci.tagMeshFirstRound.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);
  PRECICE_DEBUG("Compute Mapping for Tagging");

  computeMapping();
  PRECICE_DEBUG("Tagging First Round");

  // Determine the Mesh to Tag
  mesh::PtrMesh origins;
  if (hasConstraint(CONSERVATIVE)) {
    origins = output();
  } else {
    origins = input();
  }

  // Gather all vertices to be tagged in a first phase.
  // max_count is used to shortcut if all vertices have been tagged.
  std::unordered_set<int> tagged;
  const std::size_t       max_count = origins->vertices().size();

  for (const Polation &interpolation : _interpolations) {
    for (const auto &elem : interpolation.getWeightedElements()) {
      if (!math::equals(elem.weight, 0.0)) {
        tagged.insert(elem.vertexID);
      }
    }
    // Shortcut if all vertices are tagged
    if (tagged.size() == max_count) {
      break;
    }
  }

  // Now tag all vertices to be tagged in the second phase.
  for (auto &v : origins->vertices()) {
    if (tagged.count(v.getID()) == 1) {
      v.tag();
    }
  }
  PRECICE_DEBUG("First Round Tagged {}/{} Vertices", tagged.size(), max_count);

  clear();
}

void VolumeCellInterpolation::tagMeshSecondRound()
{
  PRECICE_TRACE();
  // for VCI mapping no operation needed here
}

} // namespace mapping
} // namespace precice

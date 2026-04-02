#pragma once

#include "MeshContext.hpp"
#include "logging/Logger.hpp"
#include "mesh/BoundingBox.hpp"
#include "mesh/Vertex.hpp"
#include "partition/ReceivedPartition.hpp"
#include "precice/span.hpp"

#include <functional>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace precice::impl {

/// Context for a mesh received from another participant
struct ReceivedMeshContext : MeshContext {
  /// Members that only exist for received meshes (moved from MeshContext)

  /// In case a mapping done by the solver is favored over a preCICE mapping, apply user-defined
  /// bounding-boxes.
  bool allowDirectAccess = false;

  /// setMeshAccessRegion may only be called once per mesh(context).
  /// putting this into the mesh context means that we can only call
  /// this once, regardless of combinations with just-in-time mappings
  /// or multiple such mappings. If multiples are desired, we would
  /// need to shift this into the MappingContext
  std::shared_ptr<mesh::BoundingBox> userDefinedAccessRegion;

  /// Name of participant that provides the mesh
  std::string receiveMeshFrom;

  /// bounding box to speed up decomposition of received mesh is increased by this safety factor
  double safetyFactor = -1;

  std::shared_ptr<precice::partition::ReceivedPartition> partition;

  /// type of geometric filter
  partition::ReceivedPartition::GeometricFilter geoFilter = partition::ReceivedPartition::GeometricFilter::UNDEFINED;

  /// Checks, that all vertices are within the user-defined access region and throws
  /// an error if vertices are not. The function does not return the result to be able
  /// to log the actual outliers
  void checkVerticesInsideAccessRegion(precice::span<const double> coordinates, int meshDim, std::string_view functionName) const;

  /// Filters vertices to those within the local access region
  std::vector<std::reference_wrapper<const mesh::Vertex>> filterVerticesToLocalAccessRegion(bool requiresBB) const;

  mutable logging::Logger _log{"impl::ReceivedMeshContext"};
};

} // namespace precice::impl

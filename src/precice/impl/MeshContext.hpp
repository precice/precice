#pragma once

#include <vector>
#include "MappingContext.hpp"
#include "SharedPointer.hpp"
#include "com/Communication.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/SharedPointer.hpp"
#include "partition/ReceivedPartition.hpp"
#include "partition/SharedPointer.hpp"

namespace precice::impl {

/// Stores a mesh and related objects and data.
struct MeshContext {
  /** Upgrades the mesh requirement to a more specific level.
    * @param[in] requirement The requirement to upgrade to.
    */
  void require(mapping::Mapping::MeshRequirement requirement);

  /// Mesh holding the geometry data structure.
  mesh::PtrMesh mesh;

  /// Determines which mesh type has to be provided by the accessor.
  mapping::Mapping::MeshRequirement meshRequirement = mapping::Mapping::MeshRequirement::UNDEFINED;

  /// Name of participant that creates the mesh.
  std::string receiveMeshFrom;

  /// bounding box to speed up decomposition of received mesh is increased by this safety factor
  double safetyFactor = -1;

  /// In case a mapping done by the solver is favored over a preCICE mapping, apply user-defined
  /// bounding-boxes.
  bool allowDirectAccess = false;

  /// setMeshAccessRegion may only be called once per mesh(context)
  /// putting this into the mesh context means that we can only call
  /// this once, regardless of combinations with just-in-time mappings
  /// or multiple such mappings
  /// if multiples are desired, we would need to shift this into the
  /// MappingContext
  std::shared_ptr<mesh::BoundingBox> userDefinedAccessRegion;

  /// True, if accessor does create the mesh.
  bool provideMesh = false;

  /// type of geometric filter
  partition::ReceivedPartition::GeometricFilter geoFilter = partition::ReceivedPartition::GeometricFilter::UNDEFINED;

  /// Partition creating the parallel decomposition of the mesh
  partition::PtrPartition partition;

  /// Mappings used when mapping data from the mesh. Can be empty.
  std::vector<MappingContext> fromMappingContexts;

  /// Mappings used when mapping data to the mesh. Can be empty.
  std::vector<MappingContext> toMappingContexts;

  void clearMappings()
  {
    for (auto &mc : fromMappingContexts) {
      mc.mapping->clear();
    }
    for (auto &mc : toMappingContexts) {
      mc.mapping->clear();
    }
  }
};

inline void MeshContext::require(mapping::Mapping::MeshRequirement requirement)
{
  meshRequirement = std::max(meshRequirement, requirement);
}

} // namespace precice::impl

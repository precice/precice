#pragma once

#include <vector>
#include "MappingContext.hpp"
#include "SharedPointer.hpp"
#include "com/Communication.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/SharedPointer.hpp"
#include "partition/SharedPointer.hpp"

namespace precice::impl {

/** Base class storing common mesh-related objects and data.
 * ProvidedMeshContext and ReceivedMeshContext derive from this.
 */
struct MeshContext {
  /** Upgrades the mesh requirement to a more specific level.
   * @param[in] requirement The requirement to upgrade to.
   */
  void require(mapping::Mapping::MeshRequirement requirement);

  /// Mesh holding the geometry data structure.
  mesh::PtrMesh mesh;

  /// Determines which mesh type has to be provided by the accessor.
  mapping::Mapping::MeshRequirement meshRequirement = mapping::Mapping::MeshRequirement::UNDEFINED;

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

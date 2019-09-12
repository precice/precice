#pragma once

#include "MappingContext.hpp"
#include "partition/SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"
#include "com/Communication.hpp"
#include "mapping/Mapping.hpp"
#include "SharedPointer.hpp"
#include "partition/ReceivedPartition.hpp"
#include <vector>

namespace precice {
namespace impl {

/// Stores a mesh and related objects and data.
struct MeshContext
{
  MeshContext ( int dimensions )
   :
    localOffset ( Eigen::VectorXd::Zero(dimensions) )
  {}

   /** Upgrades the mesh requirement to a more specific level.
    * @param[in] requirement The requirement to upgrade to.
    */
   void require(mapping::Mapping::MeshRequirement requirement);

   /// Mesh holding the geometry data structure.
   mesh::PtrMesh mesh;

   /// Data IDs of properties the geometry does possess.
   std::vector<int> associatedData;

   /// Determines which mesh type has to be provided by the accessor.
   mapping::Mapping::MeshRequirement meshRequirement = mapping::Mapping::MeshRequirement::UNDEFINED;

   /// Name of participant that creates the mesh.
   std::string receiveMeshFrom;

   /// bounding box to speed up decomposition of received mesh is increased by this safety factor
   double safetyFactor = -1;

   /// True, if accessor does create the mesh.
   bool provideMesh = false;

   /// Name of the patch as called in the adapter
   std::string patchName;

   /// type of geometric filter
   partition::ReceivedPartition::GeometricFilter geoFilter = partition::ReceivedPartition::GeometricFilter::UNDEFINED;

   /// Offset only applied to meshes local to the accessor.
   Eigen::VectorXd localOffset;

   /// Partition creating the parallel decomposition of the mesh
   partition::PtrPartition partition;

   /// Mapping used when mapping data from the mesh. Can be empty.
   MappingContext fromMappingContext;

   /// Mapping used when mapping data to the mesh. Can be empty.
   MappingContext toMappingContext;

};

inline void MeshContext::require(mapping::Mapping::MeshRequirement requirement) {
    meshRequirement = std::max(meshRequirement, requirement);
}

}} // namespace precice, impl

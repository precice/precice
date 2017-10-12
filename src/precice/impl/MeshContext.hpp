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

   // @brief Mesh holding the geometry data structure.
   mesh::PtrMesh mesh;

  // @brief Data IDs of properties the geometry does posses.
   std::vector<int> associatedData;

   // @brief Determines which mesh type has to be provided by the accessor.
   mapping::Mapping::MeshRequirement meshRequirement;

   // @brief Name of participant creating the geometry the mesh.
   std::string receiveMeshFrom;

   // @brief bounding box to speed up decomposition of received mesh is increased by this safety factor
   double safetyFactor;

   // @brief True, if accessor does create the geometry of the mesh.
   bool provideMesh;

   /// type of geometric filter
   partition::ReceivedPartition::GeometricFilter geoFilter;

   /// Offset only applied to meshes local to the accessor.
   Eigen::VectorXd localOffset;

   // @brief Partition creating the parallel decomposition of the mesh
   partition::PtrPartition partition;

   // @brief Mapping used when mapping data from the mesh. Can be empty.
   MappingContext fromMappingContext;

   // @brief Mapping used when mapping data to the mesh. Can be empty.
   MappingContext toMappingContext;

   MeshContext ( int dimensions )
   :
     mesh (),
     associatedData (),
     meshRequirement ( mapping::Mapping::UNDEFINED ),
     receiveMeshFrom ( "" ),
     safetyFactor(-1.0),
     provideMesh ( false ),
     geoFilter (partition::ReceivedPartition::GeometricFilter::UNDEFINED),
     localOffset ( Eigen::VectorXd::Zero(dimensions) ),
     partition (),
     fromMappingContext(),
     toMappingContext()
   {}
};

}} // namespace precice, impl

// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IMPL_MESHCONTEXT_HPP_
#define PRECICE_IMPL_MESHCONTEXT_HPP_

#include "MappingContext.hpp"
#include "geometry/SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"
#include "spacetree/SharedPointer.hpp"
#include "com/SharedPointer.hpp"
#include "mapping/Mapping.hpp"
#include "SharedPointer.hpp"
#include "utils/Dimensions.hpp"
#include <vector>

namespace precice {
namespace impl {

/**
 * @brief Stores a mesh and related objects and data.
 */
struct MeshContext
{

   // @brief Mesh holding the geometry data structure.
   mesh::PtrMesh mesh;

   // @brief Spacetree accelerating access to geometry data structure.
   spacetree::PtrSpacetree spacetree;

   // @brief Data IDs of properties the geometry does posses.
   std::vector<int> associatedData;

   // @brief Determines which mesh type has to be provided by the accessor.
   mapping::Mapping::MeshRequirement meshRequirement;

   // @brief Name of participant creating the geometry the mesh.
   std::string receiveMeshFrom;

   // @brief True, if accessor does create the geometry of the mesh.
   bool provideMesh;

   // @brief True, if accessor has to provide write positions on initialization.
   bool provideWritePositions;

   // @brief True, if accessor has to provide read positions on initialization.
   bool provideReadPositions;

   // @brief Offset only applied to meshes local to the accessor.
   utils::DynVector localOffset;

   // @brief Geometry creating the mesh. Can be empty.
   geometry::PtrGeometry geometry;

   // @brief Mapping used when mapping data from the mesh. Can be empty.
  MappingContext mappingContext;

   /**
    * @brief Constructor.
    */
   MeshContext ( int dimensions )
   :
     mesh (),
     spacetree (),
     associatedData (),
     meshRequirement ( mapping::Mapping::UNDEFINED ),
     receiveMeshFrom ( "" ),
     provideMesh ( false ),
     localOffset ( dimensions, 0.0 ),
     geometry (),
     mappingContext()
   {}
};

}} // namespace precice, impl

#endif /* PRECICE_IMPL_MESHCONTEXT_HPP_ */

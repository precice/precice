// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "VoxelPosition.hpp"
#include "mesh/Group.hpp"
#include "mesh/PropertyContainer.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "utils/Globals.hpp"
#include "spacetree/Spacetree.hpp"

namespace precice {
  namespace impl {

    /**
     * @brief Hides private members from interface of VoxelPosition.
     */
    struct VoxelPositionImplementation
    {
      // @brief Geometry IDs from the geometries contained in the voxel.
      std::vector<int> meshIDs;

      // @brief Position of the query voxel relative to the geometry.
      int position;

      // @brief Content of the voxel.
      mesh::PtrGroup content;

      VoxelPositionImplementation ()
      :
        meshIDs (),
        position (spacetree::Spacetree::positionUndefined()),
        content ( mesh::PtrGroup(new mesh::Group) )
      {}
    };

  }
}

namespace precice {

VoxelPosition:: VoxelPosition ()
:
  _impl ( new impl::VoxelPositionImplementation() )
{}

VoxelPosition:: VoxelPosition
(
  int                   position,
  const mesh::PtrGroup& content)
:
  _impl ( new impl::VoxelPositionImplementation() )
{
  _impl->content = content;
  _impl->position = position;
  if ( position == constants::positionOnGeometry()) {
    assertion ( content->size() > 0 );
    int geoID = mesh::PropertyContainer::INDEX_GEOMETRY_ID;
    std::vector<int> ids;
    for ( mesh::Vertex & vertex : content->vertices() ) {
      vertex.getProperties ( geoID, ids );
    }
    for ( mesh::Edge & edge : content->edges() ) {
      edge.getProperties ( geoID, ids );
    }
    for ( mesh::Triangle & triangle : content->triangles() ) {
      triangle.getProperties ( geoID, ids );
    }
    _impl->meshIDs.clear ();
    for ( int id : ids ) {
      if ( ! utils::contained(id, _impl->meshIDs) ) {
        _impl->meshIDs.push_back ( id );
      }
    }
  }
}

VoxelPosition:: VoxelPosition
(
  int                      position,
  const std::vector<int> & meshIDs )
:
  _impl ( new impl::VoxelPositionImplementation() )
{
  _impl->content = mesh::PtrGroup ( new mesh::Group() );
  _impl->position = position;
  _impl->meshIDs = meshIDs;
}

VoxelPosition:: VoxelPosition
(
  const VoxelPosition & toCopy )
:
  _impl ( new impl::VoxelPositionImplementation(*toCopy._impl) )
{}

VoxelPosition & VoxelPosition:: operator=
(
  const VoxelPosition & toAssign )
{
  _impl = new impl::VoxelPositionImplementation ( *toAssign._impl );
  return *this;
}

VoxelPosition:: ~VoxelPosition ()
{
  assertion ( _impl != nullptr );
  delete _impl;
}

std::vector<int> & VoxelPosition:: meshIDs ()
{
  return _impl->meshIDs;
}

int VoxelPosition:: position ()
{
  return _impl->position;
}

void VoxelPosition:: setPosition
(
  int position )
{
  _impl->position = position;
}

MeshHandle VoxelPosition:: contentHandle ()
{
  return MeshHandle ( *(_impl->content) );
}

} // namespace precice


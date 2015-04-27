// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_VOXELPOSITION_HPP_
#define PRECICE_VOXELPOSITION_HPP_

#include "Constants.hpp"
#include "MeshHandle.hpp"
#include <vector>
#include <memory>

namespace precice {
  namespace impl {
    struct VoxelPositionImplementation;
  }
  namespace mesh {
    class Group;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {

/**
 * @brief Holds information about a query voxel in relation to geometry.
 *
 * Is returned from the function inquireVoxelPosition().
 */
class VoxelPosition
{
public:

  /**
   * @brief Default constructor, creates an empty position.
   */
  VoxelPosition ();

  /**
   * @brief Constructor, does set a content and retrieves mesh IDs from content.
   */
  VoxelPosition (
    int                                    position,
    const std::shared_ptr<mesh::Group> &   content );

  /**
   * @brief Constructor, does not set any content.
   */
  VoxelPosition (
    int                      position,
    const std::vector<int> & meshIDs );

  /**
   * @brief Copy constructor.
   */
  VoxelPosition ( const VoxelPosition & toCopy );

  /**
   * @brief Assignment constructor.
   */
  VoxelPosition& operator= ( const VoxelPosition & toAssign );

  /**
   * @brief Destructor, frees resources.
   */
  ~VoxelPosition ();

  /**
   * @brief Geometry IDs from the geometries contained in the voxel.
   */
  std::vector<int>& meshIDs();

  /**
   * @brief Position of the query voxel relative to the geometry.
   *
   * The meaning of the returned value can be decoded by using the
   * precice::constants::position...() methods.
   */
  int position();

  /**
   * @brief Sets the query voxel position relative to the geometry.
   */
  void setPosition ( int position );

  /**
   * @brief Returns a handle to the content of the voxel.
   *
   * ATTENTION: the handle is invalidated when the VoxelPosition is destroyed.
   */
  MeshHandle contentHandle();

private:

  impl::VoxelPositionImplementation * _impl;
};

} // namespace precice

#endif /* PRECICE_VOXELPOSITION_HPP_ */

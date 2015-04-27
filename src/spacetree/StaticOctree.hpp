// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_NEWSPACETREE_STATICOCTREE_HPP_
#define PRECICE_NEWSPACETREE_STATICOCTREE_HPP_

#include "spacetree/Spacetree.hpp"
#include "spacetree/impl/OctreeCell.hpp"
#include "utils/Dimensions.hpp"
#include "tarch/logging/Log.h"

namespace precice {
namespace spacetree {

/**
 * @brief Quad- (2D), Octree (3D) that is built completely on initialization.
 */
class StaticOctree : public Spacetree
{
public:

  /**
   * @brief Constructor.
   *
   * @param center [IN] Center of the root cell.
   * @param halflength [IN] Half sidelength of the root cell.
   * @param refinementLimit [IN] Fully refined cell is guaranteed to have smaller
   *        or equal sidelengths.
   */
  StaticOctree (
    const utils::DynVector& center,
    double halflength,
    double refinementLimit );

  /**
   * @brief Destructor, empty.
   */
  virtual ~StaticOctree() {}

  /**
   * @brief Adds a mesh to the spacetree.
   *
   * The spacetree is rebuilt on the next search call.
   */
  virtual void addMesh ( const mesh::PtrMesh& content );

  /**
   * @brief Returns the meshes added to the spacetree.
   */
  virtual const std::vector<mesh::PtrMesh>& meshes() const;

  /**
   * @brief (Re)Initializes the spacetree structure.
   *
   * Is called automatically, after meshChanged was called, and a search method
   * is envoked.
   */
  virtual void initialize();

  /**
   * @brief Sets flag, which leads to a rebuild of the tree when issuing search.
   */
  virtual void meshChanged ( mesh::Mesh& mesh );

  /**
   * @brief Searches for the position (inside/outside/on) of a point in space.
   */
  virtual int searchPosition ( const utils::DynVector& point );

  /**
   * @brief Searches for the distance to the next mesh.
   */
  virtual void searchDistance ( query::FindClosest& findClosest );

  /**
   * @brief Searches for the content and position of a hexahedron.
   */
  virtual int searchContent ( query::FindVoxelContent& find );

  /**
   * @brief Visitor entry to visit all cells of the spacetree.
   */
  virtual void accept ( Visitor& visitor );

  /**
   * @brief Clears the spacetree.
   */
  virtual void clear();

private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  // @brief Contained meshes.
  std::vector<mesh::PtrMesh> _meshes;

  // @brief Center of root cell.
  utils::DynVector _center;

  // @brief Half sidelengths of root cell.
  double _halflength;

  // @brief Guaranteed max sidelengths of fully refined cell.
  double _refinementLimit;

  // @brief Root cell data structure.
  impl::OctreeCell _rootCell;

  // @brief Flag to signal that a contained mesh has changed.
  bool _meshChanged;
};

}} // namespace precice, spacetree

#endif /* PRECICE_NEWSPACETREE_STATICOCTREE_HPP_ */

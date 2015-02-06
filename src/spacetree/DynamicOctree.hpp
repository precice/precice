// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_NEWSPACETREE_DYNAMICOCTREE_HPP_
#define PRECICE_NEWSPACETREE_DYNAMICOCTREE_HPP_

#include "spacetree/Spacetree.hpp"
#include "spacetree/impl/OctreeCell.hpp"
#include "utils/Dimensions.hpp"
#include "tarch/logging/Log.h"
#include <vector>

namespace precice {
namespace spacetree {

/**
 * @brief Quad- (2D), DynamicOctree (3D).
 */
class DynamicOctree : public Spacetree
{
public:

  DynamicOctree (
    const utils::DynVector& center,
    double halflength,
    double refinementLimit );

  virtual ~DynamicOctree() {}

  virtual void addMesh (  const mesh::PtrMesh& content );

  virtual const std::vector<mesh::PtrMesh>& meshes() const;

  virtual void initialize();

  virtual void meshChanged ( mesh::Mesh& mesh );

  virtual int searchPosition ( const utils::DynVector& point );

  /**
   * @brief Spacetree accelerated distance search.
   */
  virtual void searchDistance ( query::FindClosest& findClosest );

  /**
   * @brief Spacetree accelerated voxel content search.
   */
  virtual int searchContent ( query::FindVoxelContent& find );

  virtual void accept ( Visitor& visitor );

  virtual void clear();

private:

  static tarch::logging::Log _log;

  // @brief Contained meshes.
  std::vector<mesh::PtrMesh> _meshes;

  utils::DynVector _center;

  double _halflength;

  double _refinementLimit;

  impl::OctreeCell _rootCell;

  bool _meshChanged;
};

}} // namespace precice, spacetree

#endif /* PRECICE_NEWSPACETREE_DYNAMICOCTREE_HPP_ */

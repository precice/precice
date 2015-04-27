// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "DynamicOctree.hpp"
#include "spacetree/impl/DynamicTraversal.hpp"
#include "query/FindVoxelContent.hpp"
#include "query/FindClosest.hpp"

namespace precice {
namespace spacetree {

tarch::logging::Log DynamicOctree:: _log("precice::spacetree::DynamicOctree");

DynamicOctree:: DynamicOctree
(
  const utils::DynVector& center,
  double halflength,
  double refinementLimit )
:
  _meshes(),
  _center(center),
  _halflength(halflength),
  _refinementLimit(refinementLimit),
  _rootCell(),
  _meshChanged(true)
{}

void DynamicOctree:: addMesh
(
  const mesh::PtrMesh& mesh )
{
  preciceTrace1("addMesh()", mesh->getName());
  assertion(_rootCell.content().empty()); // Spacetree is not initialized yet
  _meshes.push_back(mesh);
  mesh->addListener(*this);
}

const std::vector<mesh::PtrMesh>& DynamicOctree:: meshes() const
{
  return _meshes;
}

void DynamicOctree:: initialize()
{
  preciceTrace("initialize()");
  assertion(_rootCell.content().empty()); // Spacetree is not initialized yet
  int dim = _center.size();
  query::FindVoxelContent findVoxel ( _center, utils::DynVector(dim,_halflength),
      query::FindVoxelContent::INCLUDE_BOUNDARY );
  int size = 0;
  foreach (mesh::PtrMesh mesh, _meshes){
    assertion2(mesh->getDimensions() == dim, mesh->getDimensions(), dim);
    size += mesh->content().size();
    findVoxel(*mesh);
  }
  _rootCell.content().add(findVoxel.content());
  _rootCell.setPosition(positionOnGeometry());
  preciceCheck((int)_rootCell.content().size() == size, "initialize()",
               "Not all meshes are contained in the spacetree!");
  _meshChanged = false;
}

void DynamicOctree:: meshChanged
(
  mesh::Mesh& mesh )
{
  preciceTrace1 ( "meshChanged()", mesh.getName() );
  _meshChanged = true;
}

int DynamicOctree:: searchPosition
(
  const utils::DynVector& point )
{
  preciceTrace1("searchPosition()", point);
  if (_meshChanged){
    preciceDebug("A mesh has changed recently, rebuilding spacetree");
    clear();
    initialize();
  }
  impl::DynamicTraversal<impl::OctreeCell> traversal;
  utils::DynVector halflengths(point.size(), _halflength);
  return traversal.searchPosition ( _rootCell, point, _center, halflengths,
                                    _refinementLimit );
}

void DynamicOctree:: searchDistance
(
  query::FindClosest& findClosest )
{
  preciceTrace1 ( "searchDistance()", findClosest.getSearchPoint() );
  if (_meshChanged){
    preciceDebug("A mesh has changed recently, rebuilding spacetree");
    clear();
    initialize();
  }
  impl::DynamicTraversal<impl::OctreeCell> traversal;
  utils::DynVector halflengths(findClosest.getSearchPoint().size(), _halflength);
  traversal.searchDistance ( _rootCell, findClosest, _center, halflengths,
                             _refinementLimit );
}

int DynamicOctree:: searchContent
(
  query::FindVoxelContent& findContent )
{
  preciceTrace2 ( "searchContent()", findContent.getVoxelCenter(),
                  findContent.getVoxelHalflengths() );
  if (_meshChanged){
    preciceDebug("A mesh has changed recently, rebuilding spacetree");
    clear();
    initialize();
  }
  impl::DynamicTraversal<impl::OctreeCell> traversal;
  utils::DynVector halflengths(findContent.getVoxelCenter().size(), _halflength);
  return traversal.searchContent ( _rootCell, findContent, _center, halflengths,
                                   _refinementLimit );
}

void DynamicOctree:: accept ( Visitor& visitor )
{
  preciceTrace("accept()");
  if (_meshChanged){
    preciceDebug("A mesh has changed recently, rebuilding spacetree");
    clear();
    initialize();
  }
  utils::DynVector halflengths(_center.size(), _halflength);
  _rootCell.accept(visitor, _center, halflengths);
}

void DynamicOctree:: clear()
{
  preciceTrace("clear()");
  _rootCell.clear();
  assertion(_rootCell.content().empty());
}

}} // namespace precice, spacetree

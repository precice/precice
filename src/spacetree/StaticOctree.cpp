// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "StaticOctree.hpp"
#include "spacetree/impl/StaticTraversal.hpp"
#include "spacetree/impl/Environment.hpp"
#include "query/FindVoxelContent.hpp"
#include "query/FindClosest.hpp"

namespace precice {
namespace spacetree {

tarch::logging::Log StaticOctree:: _log("precice::spacetree::StaticOctree");

StaticOctree:: StaticOctree
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

void StaticOctree:: addMesh
(
  const mesh::PtrMesh& mesh )
{
  preciceTrace1("addMesh()", mesh->getName());
  assertion(_rootCell.content().empty()); // Spacetree is not initialized yet
  _meshes.push_back(mesh);
  mesh->addListener(*this);
}

const std::vector<mesh::PtrMesh>& StaticOctree:: meshes() const
{
  return _meshes;
}

void StaticOctree:: initialize()
{
  preciceTrace("initialize()");
  assertion(_rootCell.content().empty());
  int dim = _center.size();
  query::FindVoxelContent findVoxel ( _center, utils::DynVector(dim,_halflength),
      query::FindVoxelContent::INCLUDE_BOUNDARY );
  size_t size = 0;
  foreach (mesh::PtrMesh mesh, _meshes){
    assertion2(mesh->getDimensions() == dim, mesh->getDimensions(), dim);
    size += mesh->content().size();
    findVoxel(*mesh);
  }
  _rootCell.content().add(findVoxel.content());
  _rootCell.setPosition(positionOnGeometry());
  preciceCheck(_rootCell.content().size() == size, "initialize()",
               "Not all meshes are contained in the spacetree!");
  impl::StaticTraversal<impl::OctreeCell> traversal;
  utils::DynVector halflengths(dim, _halflength);
  int twoPowerDim = std::pow(2.0, dim);
  int sides = (dim == 2) ? 4 : 6;
  impl::Environment env(twoPowerDim, sides);
  if ( dim == 2 ){
    preciceDebug( "Setting 2D environment cell neighbor indices" );
    tarch::la::DynamicVector<int> indices(2);
    indices[0] = 1; indices[1] = 2; // Neighbors cell 0
    env.setNeighborCellIndices(0, indices);
    indices[0] = 0; indices[1] = 3; // Neighbors cell 1
    env.setNeighborCellIndices(1, indices);
    indices[0] = 3; indices[1] = 0; // Neighbors cell 2
    env.setNeighborCellIndices(2, indices);
    indices[0] = 2; indices[1] = 1; // Neighbors cell 3
    env.setNeighborCellIndices(3, indices);

    assignList(indices) = 1, 3; // Sides cell 0
    env.setNeighborSideIndices(0, indices);
    assignList(indices) = 0, 3; // Sides cell 1
    env.setNeighborSideIndices(1, indices);
    assignList(indices) = 1, 2; // Sides cell 2
    env.setNeighborSideIndices(2, indices);
    assignList(indices) = 0, 2; // Sides cell 3
    env.setNeighborSideIndices(3, indices);
  }
  else {
    preciceDebug( "Setting 3D environment cell neighbor indices" );
    assertion1 ( dim == 3, dim );
    tarch::la::DynamicVector<int> indices(3);
    assignList(indices) = 1, 2, 4; // Cell 0
    env.setNeighborCellIndices(0, indices);
    assignList(indices) = 0, 3, 5; // Cell 1
    env.setNeighborCellIndices(1, indices);
    assignList(indices) = 3, 0, 6; // Cell 2
    env.setNeighborCellIndices(2, indices);
    assignList(indices) = 2, 1, 7; // Cell 3
    env.setNeighborCellIndices(3, indices);
    assignList(indices) = 5, 6, 0; // Cell 4
    env.setNeighborCellIndices(4, indices);
    assignList(indices) = 4, 7, 1; // Cell 5
    env.setNeighborCellIndices(5, indices);
    assignList(indices) = 7, 4, 2; // Cell 6
    env.setNeighborCellIndices(6, indices);
    assignList(indices) = 6, 5, 3; // Cell 7
    env.setNeighborCellIndices(7, indices);

    assignList(indices) = 1, 3, 5; // Side 0
    env.setNeighborSideIndices(0, indices);
    assignList(indices) = 0, 3, 5; // Side 1
    env.setNeighborSideIndices(1, indices);
    assignList(indices) = 1, 2, 5; // Side 2
    env.setNeighborSideIndices(2, indices);
    assignList(indices) = 0, 2, 5; // Side 3
    env.setNeighborSideIndices(3, indices);
    assignList(indices) = 1, 3, 4; // Side 4
    env.setNeighborSideIndices(4, indices);
    assignList(indices) = 0, 3, 4; // Side 5
    env.setNeighborSideIndices(5, indices);
    assignList(indices) = 1, 2, 4; // Side 6
    env.setNeighborSideIndices(6, indices);
    assignList(indices) = 0, 2, 4; // Side 7
    env.setNeighborSideIndices(7, indices);
  }
  traversal.refineAll(_rootCell, _center, halflengths, _refinementLimit, env);
  _meshChanged = false;
}

void StaticOctree:: meshChanged ( mesh::Mesh& mesh )
{
  preciceTrace1("meshChanged()", mesh.getName());
  _meshChanged = true;
}

int StaticOctree:: searchPosition
(
  const utils::DynVector& point )
{
  preciceTrace1 ( "searchPosition()", point );
  if (_meshChanged){
    preciceDebug("A mesh has changed recently, rebuilding spacetree");
    clear();
    initialize();
  }
  impl::StaticTraversal<impl::OctreeCell> traversal;
  utils::DynVector halflengths(point.size(), _halflength);
  return traversal.searchPosition ( _rootCell, point, _center, halflengths );
}

void StaticOctree:: searchDistance
(
  query::FindClosest& findClosest )
{
  preciceTrace1 ( "searchDistance()", findClosest.getSearchPoint() );
  if (_meshChanged){
    preciceDebug("A mesh has changed recently, rebuilding spacetree");
    clear();
    initialize();
  }
  impl::StaticTraversal<impl::OctreeCell> traversal;
  utils::DynVector halflengths(findClosest.getSearchPoint().size(), _halflength);
  traversal.searchDistance ( _rootCell, findClosest, _center, halflengths );
}

int StaticOctree:: searchContent
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
  impl::StaticTraversal<impl::OctreeCell> traversal;
  utils::DynVector halflengths(findContent.getVoxelCenter().size(), _halflength);
  return traversal.searchContent ( _rootCell, findContent, _center, halflengths );
}

void StaticOctree:: accept ( Visitor& visitor )
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

void StaticOctree:: clear()
{
  preciceTrace("clear()");
  _rootCell.clear();
  assertion(_rootCell.content().empty());
}

}} // namespace precice, spacetree

#include "DynamicOctree.hpp"
#include "spacetree/impl/DynamicTraversal.hpp"
#include "query/FindVoxelContent.hpp"
#include "query/FindClosest.hpp"

namespace precice {
namespace spacetree {

logging::Logger DynamicOctree::_log("precice::spacetree::DynamicOctree");

DynamicOctree:: DynamicOctree
(
  const Eigen::VectorXd& center,
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
  TRACE(mesh->getName());
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
  TRACE();
  assertion(_rootCell.content().empty()); // Spacetree is not initialized yet
  int dim = _center.size();
  query::FindVoxelContent findVoxel(_center,
                                    Eigen::VectorXd::Constant(dim,_halflength),
                                    query::FindVoxelContent::INCLUDE_BOUNDARY );
  int size = 0;
  for (mesh::PtrMesh mesh : _meshes){
    assertion(mesh->getDimensions() == dim, mesh->getDimensions(), dim);
    size += mesh->content().size();
    findVoxel(*mesh);
  }
  _rootCell.content().add(findVoxel.content());
  _rootCell.setPosition(positionOnGeometry());
  CHECK((int)_rootCell.content().size() == size, "Not all meshes are contained in the spacetree!");
  _meshChanged = false;
}

void DynamicOctree:: meshChanged
(
  mesh::Mesh& mesh )
{
  TRACE(mesh.getName() );
  _meshChanged = true;
}

int DynamicOctree:: searchPosition
(
  const Eigen::VectorXd& point )
{
  TRACE(point);
  if (_meshChanged){
    DEBUG("A mesh has changed recently, rebuilding spacetree");
    clear();
    initialize();
  }
  impl::DynamicTraversal<impl::OctreeCell> traversal;
  Eigen::VectorXd halflengths = Eigen::VectorXd::Constant(point.size(), _halflength);

  return traversal.searchPosition ( _rootCell, point, _center, halflengths,
                                    _refinementLimit );
}

void DynamicOctree:: searchDistance
(
  query::FindClosest& findClosest )
{
  TRACE(findClosest.getSearchPoint() );
  if (_meshChanged){
    DEBUG("A mesh has changed recently, rebuilding spacetree");
    clear();
    initialize();
  }
  impl::DynamicTraversal<impl::OctreeCell> traversal;
  Eigen::VectorXd halflengths = Eigen::VectorXd::Constant(findClosest.getSearchPoint().size(), _halflength);

  traversal.searchDistance ( _rootCell, findClosest, _center, halflengths,
                             _refinementLimit );
}

int DynamicOctree:: searchContent
(
  query::FindVoxelContent& findContent )
{
  TRACE(findContent.getVoxelCenter(), findContent.getVoxelHalflengths() );
  if (_meshChanged){
    DEBUG("A mesh has changed recently, rebuilding spacetree");
    clear();
    initialize();
  }
  impl::DynamicTraversal<impl::OctreeCell> traversal;
  Eigen::VectorXd halflengths = Eigen::VectorXd::Constant(findContent.getVoxelCenter().size(), _halflength);
  return traversal.searchContent ( _rootCell, findContent, _center, halflengths,
                                   _refinementLimit );
}

void DynamicOctree:: accept ( Visitor& visitor )
{
  TRACE();
  if (_meshChanged){
    DEBUG("A mesh has changed recently, rebuilding spacetree");
    clear();
    initialize();
  }
  Eigen::VectorXd halflengths = Eigen::VectorXd::Constant(_center.size(), _halflength);
  _rootCell.accept(visitor, _center, halflengths);
}

void DynamicOctree:: clear()
{
  TRACE();
  _rootCell.clear();
  assertion(_rootCell.content().empty());
}

}} // namespace precice, spacetree


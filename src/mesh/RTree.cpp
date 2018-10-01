#include "RTree.hpp"

namespace precice {
namespace mesh {

// Initialize static member
std::map<int, rtree::PtrRTree> precice::mesh::rtree::trees;
// Initialize static member
std::map<int, PtrPrimitiveRTree> precice::mesh::rtree::primitive_trees_;

rtree::PtrRTree rtree::getVertexRTree(PtrMesh mesh)
{
  RTreeParameters params;
  VertexIndexGetter ind(mesh->vertices());
    
  auto result = trees.emplace(std::piecewise_construct,
                              std::forward_as_tuple(mesh->getID()),
                              std::forward_as_tuple(std::make_shared<VertexRTree>(params, ind)));
    
  PtrRTree tree = std::get<0>(result)->second;
    
  if (std::get<1>(result)) // insertion took place, fill tree
    for (size_t i = 0; i < mesh->vertices().size(); ++i)
      tree->insert(i);
    
  return tree;
}

PtrPrimitiveRTree rtree::getPrimitiveRTree(PtrMesh mesh)
{
  assertion(mesh, "Empty meshes are not allowed.");
  auto iter = primitive_trees_.find(mesh->getID());
  if (iter != primitive_trees_.end()) {
    return iter->second;
  } else {
    auto treeptr = std::make_shared<PrimitiveRTree>(indexMesh(*mesh));
    primitive_trees_.insert(std::make_pair(
        mesh->getID(),
        treeptr));
    return treeptr;
  }
}

void rtree::clear(Mesh &mesh)
{
  trees.erase(mesh.getID());
  primitive_trees_.erase(mesh.getID());
}


Box3d getEnclosingBox(Vertex const & middlePoint, double sphereRadius)
{
  namespace bg = boost::geometry;
  auto & coords = middlePoint.getCoords();

  Box3d box;
  bg::set<bg::min_corner, 0>(box, bg::get<0>(coords) - sphereRadius);
  bg::set<bg::min_corner, 1>(box, bg::get<1>(coords) - sphereRadius);
  bg::set<bg::min_corner, 2>(box, bg::get<2>(coords) - sphereRadius);

  bg::set<bg::max_corner, 0>(box, bg::get<0>(coords) + sphereRadius);
  bg::set<bg::max_corner, 1>(box, bg::get<1>(coords) + sphereRadius);
  bg::set<bg::max_corner, 2>(box, bg::get<2>(coords) + sphereRadius);
  
  return box;
}

}}


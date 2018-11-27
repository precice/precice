#include "impl/RTree.hpp"

#include "RTree.hpp"

namespace precice {
namespace mesh {

// Initialize static member
std::map<int, rtree::PtrRTree> precice::mesh::rtree::trees;
// Initialize static member
std::map<int, PtrPrimitiveRTree> precice::mesh::rtree::_primitive_trees;

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
  auto iter = _primitive_trees.find(mesh->getID());
  if (iter != _primitive_trees.end()) {
    return iter->second;
  }
  auto treeptr = std::make_shared<PrimitiveRTree>(indexMesh(*mesh));
  _primitive_trees.emplace(std::piecewise_construct,
          std::forward_as_tuple(mesh->getID()),
          std::forward_as_tuple(treeptr));
  return treeptr;
}

void rtree::clear(Mesh &mesh)
{
  trees.erase(mesh.getID());
  _primitive_trees.erase(mesh.getID());
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

PrimitiveRTree indexMesh(const Mesh &mesh)
{
  using namespace impl;

  AABBGenerator  gen{mesh};
  PrimitiveRTree tree;
  indexPrimitive(tree, gen, mesh.vertices());
  indexPrimitive(tree, gen, mesh.edges());
  indexPrimitive(tree, gen, mesh.triangles());
  indexPrimitive(tree, gen, mesh.quads());
  return tree;
}

std::ostream &operator<<(std::ostream &out, Primitive val)
{
  switch (val) {
  case (Primitive::Vertex):
    out << "Vertex";
    break;
  case (Primitive::Edge):
    out << "Edge";
    break;
  case (Primitive::Triangle):
    out << "Triangle";
    break;
  case (Primitive::Quad):
    out << "Quad";
    break;
  }
  return out;
}

std::ostream& operator<<(std::ostream& out, PrimitiveIndex val) {
    return out << val.type << ":" << val.index;
}

bool operator==(const PrimitiveIndex& lhs, const PrimitiveIndex& rhs)
{
    return lhs.type == rhs.type && lhs.index == rhs.index;
}

/// Standard non-equality test for PrimitiveIndex
bool operator!=(const PrimitiveIndex& lhs, const PrimitiveIndex& rhs)
{
    return !(lhs == rhs);
}

}}


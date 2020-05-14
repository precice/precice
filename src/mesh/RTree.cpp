#include "mesh/impl/RTree.hpp"
#include "mesh/impl/RTreeAdapter.hpp"

#include <boost/range/irange.hpp>
#include "mesh/RTree.hpp"

namespace precice {
namespace mesh {

namespace bg = boost::geometry;

// Initialize static member
std::map<int, rtree::MeshIndices> precice::mesh::rtree::_cached_trees;
std::map<int, PtrPrimitiveRTree>  precice::mesh::rtree::_primitive_trees;

rtree::MeshIndices &rtree::cacheEntry(int meshID)
{
  auto result = _cached_trees.emplace(std::make_pair(meshID, rtree::MeshIndices{}));
  return result.first->second;
}

rtree::vertex_traits::Ptr rtree::getVertexRTree(const PtrMesh &mesh, int toPatchID)
{
  PRECICE_ASSERT(mesh);
  auto &cache = cacheEntry(mesh->getID());
  //int patchID = mesh->getTotalPatches();      //This actually needs to be the patch number, but totalPatches will suffice for intial testing
  if (cache.vertices) {
    return cache.vertices;
  }

  // Generating the rtree is expensive, so passing everything in the ctor is
  // the best we can do. Even passing an index range instead of calling
  // tree->insert repeatedly is about 10x faster.
  RTreeParameters            params;
  // If adding a patch class, mesh->patch()->vertices();

  for (int i = 0; i < mesh->vertices().size(); i++){
    int vertexPatchID = mesh->vertices()[i].getPatchID();
  }
  // If not using a patchVertexContainer, then check if patchID matches, then add the vertex if it does.
  // This slows down this process, but NN with boost remains O(NlogN). This is a larger savings.
  // This can also be used for NP mapping. RBF does not matter. It runs through both. Maybe include this for RBF
  // and use boost to find points N nearest points or within a radius? Speed up building RBF function.
  

  vertex_traits::IndexGetter ind(mesh->vertices());
  auto                       tree = std::make_shared<vertex_traits::RTree>(
      boost::irange<std::size_t>(0lu, mesh->vertices().size()), params, ind);

  cache.vertices = tree;
  return tree;
}

rtree::edge_traits::Ptr rtree::getEdgeRTree(const PtrMesh &mesh)
{
  PRECICE_ASSERT(mesh);
  auto &cache = cacheEntry(mesh->getID());
  if (cache.edges) {
    return cache.edges;
  }

  // Generating the rtree is expensive, so passing everything in the ctor is
  // the best we can do. Even passing an index range instead of calling
  // tree->insert repeatedly is about 10x faster.
  RTreeParameters          params;
  edge_traits::IndexGetter ind(mesh->edges());
  auto                     tree = std::make_shared<edge_traits::RTree>(
      boost::irange<std::size_t>(0lu, mesh->edges().size()), params, ind);

  cache.edges = tree;
  return tree;
}

rtree::triangle_traits::Ptr rtree::getTriangleRTree(const PtrMesh &mesh)
{
  PRECICE_ASSERT(mesh);
  auto &cache = cacheEntry(mesh->getID());
  if (cache.triangles) {
    return cache.triangles;
  }

  // We first generate the values for the triangle rtree.
  // The resulting vector is a random access range, which can be passed to the
  // constructor of the rtree for more efficient indexing.
  std::vector<triangle_traits::IndexType> elements;
  elements.reserve(mesh->triangles().size());
  for (size_t i = 0; i < mesh->triangles().size(); ++i) {
    auto box = bg::return_envelope<RTreeBox>(mesh->triangles()[i]);
    elements.emplace_back(std::move(box), i);
  }

  // Generating the rtree is expensive, so passing everything in the ctor is
  // the best we can do.
  RTreeParameters              params;
  triangle_traits::IndexGetter ind;
  auto                         tree = std::make_shared<triangle_traits::RTree>(elements, params, ind);
  cache.triangles                   = tree;
  return tree;
}

PtrPrimitiveRTree rtree::getPrimitiveRTree(const PtrMesh &mesh)
{
  PRECICE_ASSERT(mesh, "Empty meshes are not allowed.");
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
  _cached_trees.erase(mesh.getID());
  _primitive_trees.erase(mesh.getID());
}

void rtree::clear()
{
  _cached_trees.clear();
  _primitive_trees.clear();
}

Box3d getEnclosingBox(Vertex const &middlePoint, double sphereRadius)
{
  namespace bg = boost::geometry;
  auto &coords = middlePoint.getCoords();

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

std::ostream &operator<<(std::ostream &out, PrimitiveIndex val)
{
  return out << val.type << ":" << val.index;
}

bool operator==(const PrimitiveIndex &lhs, const PrimitiveIndex &rhs)
{
  return lhs.type == rhs.type && lhs.index == rhs.index;
}

/// Standard non-equality test for PrimitiveIndex
bool operator!=(const PrimitiveIndex &lhs, const PrimitiveIndex &rhs)
{
  return !(lhs == rhs);
}

} // namespace mesh
} // namespace precice

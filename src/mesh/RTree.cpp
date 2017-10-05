#include "RTree.hpp"

namespace precice {
namespace mesh {

// Initialize static member
std::map<int, rtree::PtrRTree> precice::mesh::rtree::trees;

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

void rtree::clear()
{
  trees.clear();
}

void rtree::clear(int id)
{
  trees.erase(id);
}

}}

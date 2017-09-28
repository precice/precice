#pragma once

#include <map>
#include "mesh/impl/RTreeAdapter.hpp"
#include "mesh/Mesh.hpp"
#include <boost/geometry.hpp>

namespace precice {
namespace mesh {

class rtree {
public:
  using VertexIndexGetter = impl::PtrVectorIndexable<Mesh::VertexContainer>;
  using RTreeParameters   = boost::geometry::index::rstar<16>;
  using VertexRTree       = boost::geometry::index::rtree<Mesh::VertexContainer::container::size_type,
                                                          RTreeParameters,
                                                          VertexIndexGetter>;
  using PtrRTree = std::shared_ptr<VertexRTree>;

  /// Returns the pointer to boost::geometry::rtree for the given mesh
  /*
   * Creates and fills the tree, if it wasn't requested before, otherwise it returns the cached tree.
   */
  static PtrRTree getVertexRTree(PtrMesh mesh)
  {
    RTreeParameters params;
    VertexIndexGetter ind(mesh->vertices());
    
    auto result = trees.emplace(std::piecewise_construct,
                                std::forward_as_tuple(mesh->getID()),
                                std::forward_as_tuple(std::make_shared<VertexRTree>(params, ind)));
    
    PtrRTree rtree = std::get<0>(result)->second;
    
    if (std::get<1>(result)) // insertion took place, fill tree
      for (size_t i = 0; i < mesh->vertices().size(); ++i)
        rtree->insert(i);
    
    return rtree;
  }

  /// Clears the cache of trees, e.g. when the mesh has changed
  static void clear()
  {
    trees.clear();
  }

private:
  static std::map<int, PtrRTree> trees;
};

}}


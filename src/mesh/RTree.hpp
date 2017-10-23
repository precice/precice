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
  static PtrRTree getVertexRTree(PtrMesh mesh);
  
  /// Clears the cache of all trees, e.g. when the mesh has changed
  static void clear();

  /// Only removed the tree of that specific mesh
  static void clear(int id);

private:
  static std::map<int, PtrRTree> trees;
};


using Box3d = boost::geometry::model::box<boost::geometry::model::point<double, 3, boost::geometry::cs::cartesian>>;

/// Returns a boost::geometry box that encloses a sphere of given radius around a middle point
Box3d getEnclosingBox(Vertex const & middlePoint, double sphereRadius);

}}

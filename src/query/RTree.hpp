#pragma once

#include <Eigen/Core>
#include <boost/geometry.hpp>
#include <iosfwd>
#include <map>
#include <memory>
#include <type_traits>
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Triangle.hpp"
#include "query/impl/RTreeAdapter.hpp"

namespace precice {

namespace testing {
namespace accessors {
struct rtree;
}
} // namespace testing

namespace query {

class rtree {
public:
  using vertex_traits   = impl::RTreeTraits<mesh::Vertex>;
  using edge_traits     = impl::RTreeTraits<mesh::Edge>;
  using triangle_traits = impl::RTreeTraits<mesh::Triangle>;

  /// Returns the pointer to boost::geometry::rtree for the given mesh vertices
  /*
   * Creates and fills the tree, if it wasn't requested before, otherwise it returns the cached tree.
   */
  static vertex_traits::Ptr getVertexRTree(const mesh::PtrMesh &mesh);

  static edge_traits::Ptr getEdgeRTree(const mesh::PtrMesh &mesh);

  static triangle_traits::Ptr getTriangleRTree(const mesh::PtrMesh &mesh);

  /// Only clear the trees of that specific mesh
  static void clear(mesh::Mesh &mesh);

  /// Clear the complete cache
  static void clear();

  friend struct testing::accessors::rtree;

private:
  struct MeshIndices {
    vertex_traits::Ptr   vertices;
    edge_traits::Ptr     edges;
    triangle_traits::Ptr triangles;
  };

  static MeshIndices &cacheEntry(int MeshID);

  using RTreeCache = std::map<int, MeshIndices>;
  static RTreeCache _cached_trees; ///< Cache for all index trees
};

using Box3d = boost::geometry::model::box<boost::geometry::model::point<double, 3, boost::geometry::cs::cartesian>>;

/// Returns a boost::geometry box that encloses a sphere of given radius around a middle point
Box3d getEnclosingBox(mesh::Vertex const &middlePoint, double sphereRadius);

} // namespace query

namespace testing {
namespace accessors {
struct rtree {
  static query::rtree::RTreeCache &getCache()
  {
    return query::rtree::_cached_trees;
  }
};
} // namespace accessors
} // namespace testing

} // namespace precice

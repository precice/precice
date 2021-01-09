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
#include "utils/Statistics.hpp"

namespace precice {
namespace mesh {
class Vertex;
class Edge;
class Triangle;
} // namespace mesh
namespace testing {
namespace accessors {
class rtree;
} // namespace accessors
} // namespace testing
} // namespace precice

namespace precice {
namespace query {

using RTreeBox        = boost::geometry::model::box<Eigen::VectorXd>;
using Box3d           = boost::geometry::model::box<boost::geometry::model::point<double, 3, boost::geometry::cs::cartesian>>;
using RTreeParameters = boost::geometry::index::rstar<16>;

/// Type trait to extract information based on the type of a Primitive
template <class T>
struct PrimitiveTraits;

template <>
struct PrimitiveTraits<mesh::Vertex> {
  using MeshContainer = mesh::Mesh::VertexContainer;
};

template <>
struct PrimitiveTraits<mesh::Edge> {
  using MeshContainer = mesh::Mesh::EdgeContainer;
};

template <>
struct PrimitiveTraits<mesh::Triangle> {
  using MeshContainer = mesh::Mesh::TriangleContainer;
};

template <class Primitive>
struct IsDirectIndexable : impl::IsDirectIndexableHelper<Primitive>::type {
};

/// The type traits of a rtree based on a Primitive
template <class Primitive>
struct RTreeTraits {
  using MeshContainer      = typename PrimitiveTraits<Primitive>::MeshContainer;
  using MeshContainerIndex = typename MeshContainer::size_type;

  using IndexType = typename std::conditional<
      IsDirectIndexable<Primitive>::value,
      MeshContainerIndex,
      std::pair<RTreeBox, MeshContainerIndex>>::type;

  using IndexGetter = typename std::conditional<
      IsDirectIndexable<Primitive>::value,
      impl::VectorIndexable<MeshContainer>,
      boost::geometry::index::indexable<IndexType>>::type;

  using RTree = boost::geometry::index::rtree<IndexType, RTreeParameters, IndexGetter>;
  using Ptr   = std::shared_ptr<RTree>;
};

struct MatchType {
  double distance;
  int    index;

  MatchType() = default;
  MatchType(double d, int i)
      : distance(d), index(i){};
  constexpr bool operator<(MatchType const &other) const
  {
    return distance < other.distance;
  };
};

class rtree {
public:
  using vertex_traits   = RTreeTraits<mesh::Vertex>;
  using edge_traits     = RTreeTraits<mesh::Edge>;
  using triangle_traits = RTreeTraits<mesh::Triangle>;

  /// Returns the pointer to boost::geometry::rtree for the given mesh vertices
  /*
   * Creates and fills the tree, if it wasn't requested before, otherwise it returns the cached tree.
   */
  static vertex_traits::Ptr getVertexRTree(const mesh::PtrMesh &mesh);

  static edge_traits::Ptr getEdgeRTree(const mesh::PtrMesh &mesh);

  static triangle_traits::Ptr getTriangleRTree(const mesh::PtrMesh &mesh);

  /// Get the closest vertex/edge/triangle to a vertex, return index and distance
  static std::vector<MatchType> getClosest(const mesh::Vertex &source, const vertex_traits::Ptr &tree, const mesh::Mesh::VertexContainer &targetContainer, int n = 1);
  static std::vector<MatchType> getClosest(const mesh::Vertex &source, const edge_traits::Ptr &tree, const mesh::Mesh::EdgeContainer &targetContainer, int n = 1);
  static std::vector<MatchType> getClosest(const mesh::Vertex &source, const triangle_traits::Ptr &tree, const mesh::Mesh::TriangleContainer &targetContainer, int n = 1);

  /// Get all the vertices inside the box and surrounded by support radius
  static std::vector<size_t> getVerticesInsideBox(const Box3d &searchBox, const vertex_traits::Ptr &tree, const mesh::Mesh::VertexContainer &vertices, const mesh::Vertex &centerVertex, double supportRadius);

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

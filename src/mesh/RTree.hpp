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
#include "mesh/impl/RTreeAdapter.hpp"

namespace precice {

namespace testing {
namespace accessors {
struct rtree;
}
} // namespace testing

namespace mesh {
class Edge;
class Triangle;
class Vertex;
namespace impl {
template <typename Container>
class VectorIndexable;
} // namespace impl

/// The RTree box type
using RTreeBox = boost::geometry::model::box<Eigen::VectorXd>;

/// Type trait to extract information based on the type of a Primitive
template <class T>
struct PrimitiveTraits;

template <>
struct PrimitiveTraits<Vertex> {
  using MeshContainer = Mesh::VertexContainer;
};

template <>
struct PrimitiveTraits<Edge> {
  using MeshContainer = Mesh::EdgeContainer;
};

template <>
struct PrimitiveTraits<Triangle> {
  using MeshContainer = Mesh::TriangleContainer;
};

/// The general rtree parameter type used in precice
using RTreeParameters = boost::geometry::index::rstar<16>;

namespace impl {
template <typename Primitive>
class IsDirectIndexableHelper {
private:
  template <typename T, typename = typename std::enable_if<
                            std::is_same<
                                typename boost::geometry::traits::tag<T>::type,
                                boost::geometry::point_tag>::value,
                            std::nullptr_t>::type>
  static std::true_type test(char *);
  template <typename T, typename = typename std::enable_if<
                            std::is_same<
                                typename boost::geometry::traits::tag<T>::type,
                                boost::geometry::segment_tag>::value,
                            std::nullptr_t>::type>
  static std::true_type test(int *);
  template <typename T, typename = typename std::enable_if<
                            std::is_same<
                                typename boost::geometry::traits::tag<T>::type,
                                boost::geometry::box_tag>::value,
                            std::nullptr_t>::type>
  static std::true_type test(void *);

  template <typename T>
  static std::false_type test(...);

public:
  using type = decltype(test<Primitive>(nullptr));
};
} // namespace impl

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

class rtree {
public:
  using vertex_traits   = RTreeTraits<Vertex>;
  using edge_traits     = RTreeTraits<Edge>;
  using triangle_traits = RTreeTraits<Triangle>;

  /// Returns the pointer to boost::geometry::rtree for the given mesh vertices
  /*
   * Creates and fills the tree, if it wasn't requested before, otherwise it returns the cached tree.
   */
  static vertex_traits::Ptr getVertexRTree(const PtrMesh &mesh);

  static edge_traits::Ptr getEdgeRTree(const PtrMesh &mesh);

  static triangle_traits::Ptr getTriangleRTree(const PtrMesh &mesh);

  /// Only clear the trees of that specific mesh
  static void clear(Mesh &mesh);

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
Box3d getEnclosingBox(Vertex const &middlePoint, double sphereRadius);

} // namespace mesh

namespace testing {
namespace accessors {
struct rtree {
  static mesh::rtree::RTreeCache &getCache()
  {
    return mesh::rtree::_cached_trees;
  }
};
} // namespace accessors
} // namespace testing

} // namespace precice

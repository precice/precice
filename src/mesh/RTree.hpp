#pragma once

#include <boost/geometry.hpp>
#include <map>
#include <memory>
#include "mesh/Mesh.hpp"
#include "mesh/Quad.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/impl/RTreeAdapter.hpp"

// Forward declaration to friend the boost test struct
namespace MeshTests {
namespace RTree {
struct CacheClearing;
}
} // namespace MeshTests

namespace precice {
namespace mesh {

/** The enumeration of various primitive types.
 * \see PrimitiveIndex
 * \see as_primitve_enum
 */
enum class Primitive {
  Vertex,
  Edge,
  Triangle,
  Quad
};

/// A standard print operator for Primitive
std::ostream &operator<<(std::ostream &out, Primitive val);

/// The type traits to return the enum value of a primitive type.
template <class T>
struct as_primitive_enum {
};
template <>
struct as_primitive_enum<mesh::Vertex> {
  static constexpr Primitive value = Primitive::Vertex;
};
template <>
struct as_primitive_enum<mesh::Edge> {
  static constexpr Primitive value = Primitive::Edge;
};
template <>
struct as_primitive_enum<mesh::Triangle> {
  static constexpr Primitive value = Primitive::Triangle;
};
template <>
struct as_primitive_enum<mesh::Quad> {
  static constexpr Primitive value = Primitive::Quad;
};

/** Binds an Index and Primitive into a type
 * \see AABBGenerator
 * \see indexPrimitives
 */
struct PrimitiveIndex {
  Primitive type;
  size_t    index;
};

/// Standard equality test for PrimitiveIndex
bool operator==(const PrimitiveIndex &lhs, const PrimitiveIndex &rhs);

/// Standard non-equality test for PrimitiveIndex
bool operator!=(const PrimitiveIndex &lhs, const PrimitiveIndex &rhs);

/// A standard print operator for PrimitiveIndex
std::ostream &operator<<(std::ostream &out, PrimitiveIndex val);

/// The axis aligned bounding box based on the Vertex Type
using AABB = boost::geometry::model::box<Eigen::VectorXd>;

/// The rtree capable of indexing primitives of an entire Mesh
using PrimitiveRTree = boost::geometry::index::rtree<std::pair<AABB, PrimitiveIndex>, boost::geometry::index::rstar<16>>;

/// The shared_ptr convenience type for PrimitiveRTree
using PtrPrimitiveRTree = std::shared_ptr<PrimitiveRTree>;

/** Indexes a given mesh and returns a PrimitiveRTree holding the index
 *
 * This indexes the vertices, edges, triangles, and quads of a given Mesh and retrurns the index tree
 *
 * \param mesh the mesh to index
 *
 * \return the tree containing the indexed primitives descibed above
 */
PrimitiveRTree indexMesh(const Mesh &mesh);

/// The RTree box type
using RTreeBox = boost::geometry::model::box<Eigen::VectorXd>;

/// Type trait to extract information based on the type of a Primitive
template <class T>
struct PrimitiveTraits;

template <>
struct PrimitiveTraits<pm::Vertex> {
  using MeshContainer = Mesh::VertexContainer;
};

template <>
struct PrimitiveTraits<pm::Edge> {
  using MeshContainer = Mesh::EdgeContainer;
};

template <>
struct PrimitiveTraits<pm::Triangle> {
  using MeshContainer = Mesh::TriangleContainer;
};

template <>
struct PrimitiveTraits<pm::Quad> {
  using MeshContainer = Mesh::VertexContainer;
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
  static vertex_traits::Ptr getVertexRTree(const PtrMesh &mesh, int toPatchID);

  static edge_traits::Ptr getEdgeRTree(const PtrMesh &mesh);

  static triangle_traits::Ptr getTriangleRTree(const PtrMesh &mesh);

  /// Returns the pointer to boost::geometry::rtree for the given mesh primitives
  /*
   * Creates and fills the tree, if it wasn't requested before, otherwise it returns the cached tree.
   */
  static PtrPrimitiveRTree getPrimitiveRTree(const PtrMesh &mesh);

  /// Only clear the trees of that specific mesh
  static void clear(Mesh &mesh);

  /// Clear the complete cache
  static void clear();

  friend struct MeshTests::RTree::CacheClearing;

private:
  struct MeshIndices {
    vertex_traits::Ptr   vertices;
    edge_traits::Ptr     edges;
    triangle_traits::Ptr triangles;
  };

  static MeshIndices &cacheEntry(int MeshID);

  static std::map<int, MeshIndices>       _cached_trees;    ///< Cache for all index trees
  static std::map<int, PtrPrimitiveRTree> _primitive_trees; ///< Cache for the primitive trees
};

using Box3d = boost::geometry::model::box<boost::geometry::model::point<double, 3, boost::geometry::cs::cartesian>>;

/// Returns a boost::geometry box that encloses a sphere of given radius around a middle point
Box3d getEnclosingBox(Vertex const &middlePoint, double sphereRadius);

} // namespace mesh
} // namespace precice

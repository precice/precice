#pragma once

#include <map>
#include <memory>
#include "mesh/impl/RTreeAdapter.hpp"
#include "mesh/Mesh.hpp"
#include <boost/geometry.hpp>
#include "mesh/Triangle.hpp"
#include "mesh/Quad.hpp"

// Forward declaration to friend the boost test struct
namespace MeshTests {
namespace RTree {
struct CacheClearing;
}}

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
std::ostream& operator<<(std::ostream& out, Primitive val);

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
bool operator==(const PrimitiveIndex& lhs, const PrimitiveIndex& rhs);

/// Standard non-equality test for PrimitiveIndex
bool operator!=(const PrimitiveIndex& lhs, const PrimitiveIndex& rhs);

/// A standard print operator for PrimitiveIndex
std::ostream& operator<<(std::ostream& out, PrimitiveIndex val);

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

class rtree {
public:
  using VertexIndexGetter = impl::PtrVectorIndexable<Mesh::VertexContainer>;
  using RTreeParameters   = boost::geometry::index::rstar<16>;
  using VertexRTree       = boost::geometry::index::rtree<Mesh::VertexContainer::container::size_type,
                                                          RTreeParameters,
                                                          VertexIndexGetter>;
  using PtrRTree = std::shared_ptr<VertexRTree>;

  /// Returns the pointer to boost::geometry::rtree for the given mesh vertices
  /*
   * Creates and fills the tree, if it wasn't requested before, otherwise it returns the cached tree.
   */
  static PtrRTree getVertexRTree(PtrMesh mesh);
  
  /// Returns the pointer to boost::geometry::rtree for the given mesh primitives
  /*
   * Creates and fills the tree, if it wasn't requested before, otherwise it returns the cached tree.
   */
  static PtrPrimitiveRTree getPrimitiveRTree(PtrMesh mesh);

  /// Only clear the trees of that specific mesh
  static void clear(Mesh & mesh);

  friend struct MeshTests::RTree::CacheClearing;
  
private:
  static std::map<int, PtrPrimitiveRTree> _primitive_trees; ///< Cache for the primitive trees
  static std::map<int, PtrRTree>          _vertex_trees; ///< Cache for the vertex trees
};


using Box3d = boost::geometry::model::box<boost::geometry::model::point<double, 3, boost::geometry::cs::cartesian>>;

/// Returns a boost::geometry box that encloses a sphere of given radius around a middle point
Box3d getEnclosingBox(Vertex const & middlePoint, double sphereRadius);

}}

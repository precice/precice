#pragma once

#include <Eigen/Core>
#include <boost/geometry.hpp>
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Tetrahedron.hpp"
#include "mesh/Vertex.hpp"
#include "utils/assertion.hpp"

namespace precice::mesh {
class Triangle;
} // namespace precice::mesh

namespace pm = precice::mesh;
namespace bg = boost::geometry;

namespace boost::geometry::traits {

/// Adapts Eigen::VectorXd to boost.geometry
/*
 * This adapts every VectorXd to a 3d point. For non-existing dimensions, zero is returned.
 */
template <>
struct tag<Eigen::VectorXd> {
  using type = point_tag;
};
template <>
struct coordinate_type<Eigen::VectorXd> {
  using type = double;
};
template <>
struct coordinate_system<Eigen::VectorXd> {
  using type = cs::cartesian;
};
template <>
struct dimension<Eigen::VectorXd> : boost::mpl::int_<3> {
};

template <size_t Dimension>
struct access<Eigen::VectorXd, Dimension> {
  static double get(Eigen::VectorXd const &p)
  {
    if (Dimension >= static_cast<size_t>(p.rows())) {
      return 0;
    } else {
      return p(Dimension);
    }
  }

  static void set(Eigen::VectorXd &p, double const &value)
  {
    if (Dimension >= static_cast<size_t>(p.rows())) {
      Eigen::VectorXd tmp(Dimension);
      std::copy_n(p.data(), Dimension - 1, tmp.data());
      p = std::move(tmp);
    }
    p(Dimension) = value;
  }
};

/// Adapts Vertex::RawCoords to boost.geometry
template <>
struct tag<pm::Vertex::RawCoords> {
  using type = point_tag;
};
template <>
struct coordinate_type<pm::Vertex::RawCoords> {
  using type = double;
};
template <>
struct coordinate_system<pm::Vertex::RawCoords> {
  using type = cs::cartesian;
};
template <>
struct dimension<pm::Vertex::RawCoords> : boost::mpl::int_<3> {
};

template <size_t Dimension>
struct access<pm::Vertex::RawCoords, Dimension> {
  static double get(pm::Vertex::RawCoords const &p)
  {
    return p[Dimension];
  }

  static void set(pm::Vertex::RawCoords &p, double const &value)
  {
    p[Dimension] = value;
  }
};

BOOST_CONCEPT_ASSERT((bg::concepts::Point<pm::Vertex::RawCoords>) );

/// Provides the necessary template specialisations to adapt precice's Vertex to boost.geometry
/*
 * This adapts every Vertex to a 3d point. For non-existing dimensions, zero is returned.
 */
template <>
struct tag<pm::Vertex> {
  using type = point_tag;
};
template <>
struct coordinate_type<pm::Vertex> {
  using type = double;
};
template <>
struct coordinate_system<pm::Vertex> {
  using type = cs::cartesian;
};
template <>
struct dimension<pm::Vertex> : boost::mpl::int_<3> {
};

template <size_t Dimension>
struct access<pm::Vertex, Dimension> {
  static double get(pm::Vertex const &p)
  {
    return p.rawCoords()[Dimension];
  }

  static void set(pm::Vertex &p, double const &value)
  {
    PRECICE_UNREACHABLE("Boost.Geometry is not allowed to write to mesh::Vertex.");
  }
};

BOOST_CONCEPT_ASSERT((concepts::Point<pm::Vertex>) );

/** @brief Provides the necessary template specialisations to adapt precice's Edge to boost.geometry
 *
 * This adapts every Edge to the segment concept of boost.geometry.
 * Include impl/RangeAdapter.hpp for full support.
 */
template <>
struct tag<pm::Edge> {
  using type = segment_tag;
};
template <>
struct point_type<pm::Edge> {
  using type = pm::Vertex::RawCoords;
};

template <size_t Index, size_t Dimension>
struct indexed_access<pm::Edge, Index, Dimension> {
  static_assert((Index <= 1), "Valid Indices are {0, 1}");
  static_assert((Dimension <= 2), "Valid Dimensions are {0, 1, 2}");

  static double get(pm::Edge const &e)
  {
    return access<point_type<pm::Edge>::type, Dimension>::get(e.vertex(Index).rawCoords());
  }

  static void set(pm::Edge &e, double const &value)
  {
    PRECICE_UNREACHABLE("Boost.Geometry is not allowed to write to mesh::Edge.");
  }
};

/** @brief Provides the necessary template specialisations to adapt precice's Triangle to boost.geometry
 *
 * This adapts every Triangle to the ring concept (filled planar polygon) of boost.geometry.
 * Include impl/RangeAdapter.hpp for full support.
 */
template <>
struct tag<pm::Triangle> {
  using type = ring_tag;
};
template <>
struct point_order<pm::Triangle> {
  static const order_selector value = clockwise;
};
template <>
struct closure<pm::Triangle> {
  static const closure_selector value = open;
};

// BOOST_CONCEPT_ASSERT( (bg::concepts::Ring<pm::Triangle>));

} // namespace boost::geometry::traits

namespace precice::query {

/// The RTree box type
using RTreeBox = boost::geometry::model::box<pm::Vertex::RawCoords>;

inline RTreeBox makeBox(const pm::Vertex::RawCoords &min, const pm::Vertex::RawCoords &max)
{
  return {min, max};
}

inline pm::Vertex::RawCoords eigenToRaw(const Eigen::VectorXd &v)
{
  const auto size = v.size();
  PRECICE_ASSERT(size == 2 || size == 3, size);
  pm::Vertex::RawCoords r{0.0, 0.0, 0.0};
  std::copy_n(v.data(), size, r.data());
  return r;
}

inline Eigen::VectorXd rawToEigen(const pm::Vertex::RawCoords &v)
{
  Eigen::VectorXd r(3);
  std::copy(v.begin(), v.end(), r.data());
  return r;
}

inline RTreeBox makeBox(const Eigen::VectorXd &min, const Eigen::VectorXd &max)
{
  return {eigenToRaw(min), eigenToRaw(max)};
}

// Overload for a tetrahedron
inline RTreeBox makeBox(const precice::mesh::Tetrahedron &tetra)
{

  precice::mesh::BoundingBox box(tetra.getDimensions());
  for (int i = 0; i < 4; ++i) {
    box.expandBy(tetra.vertex(i));
  }

  // Convert to Boost type
  return makeBox(box.minCorner(), box.maxCorner());
}

namespace impl {

/// The general rtree parameter type used in precice
using RTreeParameters = boost::geometry::index::rstar<16>;

/// Type trait to extract information based on the type of a Primitive
template <class T>
struct PrimitiveTraits;

template <>
struct PrimitiveTraits<pm::Vertex> {
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

template <>
struct PrimitiveTraits<mesh::Tetrahedron> {
  using MeshContainer = mesh::Mesh::TetraContainer;
};

/// Makes a utils::PtrVector indexable and thus be usable in boost::geometry::rtree
template <typename Container>
class PtrVectorIndexable {
  using size_type = typename Container::container::size_type;
  using cref      = const typename Container::value_type &;
  Container const &container;

public:
  using result_type = cref;

  explicit PtrVectorIndexable(Container const &c)
      : container(c)
  {
  }

  result_type operator()(size_type i) const
  {
    return container[i];
  }
};

/// Makes a std::vector indexable and thus be usable in boost::geometry::rtree
template <typename Container>
class VectorIndexable {
  using size_type = typename Container::size_type;
  using cref      = const typename Container::value_type &;
  Container const &container;

public:
  using result_type = cref;

  explicit VectorIndexable(Container const &c)
      : container(c)
  {
  }

  result_type operator()(size_type i) const
  {
    return container[i];
  }
};

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

} // namespace impl
} // namespace precice::query

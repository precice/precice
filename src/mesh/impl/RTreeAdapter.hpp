#pragma once

#include <Eigen/Core>
#include <boost/geometry.hpp>
#include "mesh/Edge.hpp"
#include "mesh/Vertex.hpp"

namespace precice {
namespace mesh {
class Triangle;
} // namespace mesh
} // namespace precice

namespace pm = precice::mesh;
namespace bg = boost::geometry;

namespace boost {
namespace geometry {
namespace traits {

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
    if (Dimension >= static_cast<size_t>(p.rows()))
      return 0;

    return p[Dimension];
  }

  static void set(Eigen::VectorXd &p, double const &value)
  {
    // This handles default initialized VectorXd
    if (p.size() == 0) {
      p = Eigen::VectorXd::Zero(3);
    }
    p[Dimension] = value;
  }
};

BOOST_CONCEPT_ASSERT((bg::concepts::Point<Eigen::VectorXd>) );

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
    if (Dimension >= static_cast<size_t>(p.getDimensions()))
      return 0;

    return p.getCoords()[Dimension];
  }

  static void set(pm::Vertex &p, double const &value)
  {
    Eigen::VectorXd vec = p.getCoords();
    vec[Dimension]      = value;
    p.setCoords(vec);
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
  using type = Eigen::VectorXd;
};

template <size_t Index, size_t Dimension>
struct indexed_access<pm::Edge, Index, Dimension> {
  static_assert((Index <= 1), "Valid Indices are {0, 1}");
  static_assert((Dimension <= 2), "Valid Dimensions are {0, 1, 2}");

  static double get(pm::Edge const &e)
  {
    return access<Eigen::VectorXd, Dimension>::get(e.vertex(Index).getCoords());
  }

  static void set(pm::Edge &e, double const &value)
  {
    Eigen::VectorXd v = e.vertex(Index).getCoords();
    access<Eigen::VectorXd, Dimension>::set(v, value);
    e.vertex(Index).setCoords(std::move(v));
  }
};

/** @brief Provides the necessary template specialisations to adapt precice's Triangle to boost.geometry
*
* This adapts every Triangle to the ring concept (filled planar polygone) of boost.geometry.
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

} // namespace traits
} // namespace geometry
} // namespace boost

namespace precice {
namespace mesh {
namespace impl {

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

} // namespace impl
} // namespace mesh
} // namespace precice

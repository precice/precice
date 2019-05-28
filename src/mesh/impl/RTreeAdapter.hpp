#pragma once

#include <boost/geometry.hpp>
#include <Eigen/Core>
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"

namespace precice {
namespace mesh {
class Triangle;
class Quad;
} // namespace mesh
} // namespace precice

using precice::mesh::Edge;
using precice::mesh::Quad;
using precice::mesh::Triangle;
using precice::mesh::Vertex;

namespace boost {
namespace geometry {
namespace traits {

/// Provides the necessary template specialisations to adapt precice's Vertex to boost.geometry
/*
* This adapts every Vertex to a 3d point. For non-existing dimensions, zero is returned.
*/
template<> struct tag<Vertex>               { using type = point_tag; };
template<> struct coordinate_type<Vertex>   { using type = double; };
template<> struct coordinate_system<Vertex> { using type = cs::cartesian; };
template<> struct dimension<Vertex> : boost::mpl::int_<3> {};

template<size_t Dimension>
struct access<Vertex, Dimension>
{
  static double get(Vertex const& p)
  {
    if (Dimension > static_cast<size_t>(p.getDimensions())-1)
      return 0;
   
    return p.getCoords()[Dimension];
  }
  
  static void set(Vertex& p, double const& value)
  {
    Eigen::VectorXd vec = p.getCoords();
    vec[Dimension] = value;
    p.setCoords(vec);
  }
};

BOOST_CONCEPT_ASSERT( (concepts::Point<Vertex>));

/** @brief Provides the necessary template specialisations to adapt precice's Edge to boost.geometry
*
* This adapts every Edge to the segment concept of boost.geometry.
* Include impl/RangeAdapter.hpp for full support.
*/
template <>
struct tag<Edge> {
  using type = segment_tag;
};
template <>
struct point_type<Edge> {
  using type = Eigen::VectorXd;
};

template <size_t Index, size_t Dimension>
struct indexed_access<Edge, Index, Dimension> {
  static_assert((Index <= 1), "Valid Indices are {0, 1}");
  static_assert((Dimension <= 2), "Valid Dimensions are {0, 1, 2}");

  static double get(Edge const &e)
  {
    return access<Eigen::VectorXd, Dimension>::get(e.vertex(Index).getCoords());
  }

  static void set(Edge &e, double const &value)
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
struct tag<Triangle> {
  using type = ring_tag;
};
template <>
struct closure<Triangle> {
  static const closure_selector value = closed;
};

/** @brief Provides the necessary template specialisations to adapt precice's Quad to boost.geometry
*
* This adapts every Quad to the ring concept (filled planar polygone) of boost.geometry.
*/
template <>
struct tag<Quad> {
  using type = ring_tag;
};
template <>
struct closure<Quad> {
  static const closure_selector value = closed;
};

/// Adapts Eigen::VectorXd to boost.geometry
/*
 * This adapts every VectorXd to a 3d point. For non-existing dimensions, zero is returned.
 */
template<> struct tag<Eigen::VectorXd>               { using type = point_tag; };
template<> struct coordinate_type<Eigen::VectorXd>   { using type = double; };
template<> struct coordinate_system<Eigen::VectorXd> { using type = cs::cartesian; };
template<> struct dimension<Eigen::VectorXd> : boost::mpl::int_<3> {};

template<size_t Dimension>
struct access<Eigen::VectorXd, Dimension>
{
  static double get(Eigen::VectorXd const& p)
  {
    if (Dimension > static_cast<size_t>(p.rows())-1)
      return 0;
   
    return p[Dimension];
  }
  
  static void set(Eigen::VectorXd& p, double const& value)
  {
    // This handles default initialized VectorXd
    if (p.size() == 0) {
        p.resize(3);
    }
    p[Dimension] = value;
  }
};

BOOST_CONCEPT_ASSERT( (concepts::Point<Eigen::VectorXd>));

/// Adapts precice's Mesh::BoundingBox to boost.geometry
/*
 * Mesh::BoundingBox should be fulfilling the boost.geometry Box concept
 */
using BoundingBox = std::vector<std::pair<double, double>>;

template <>
struct tag<BoundingBox>
{
  using type = box_tag;
};

namespace bg = ::boost::geometry;
template <>
struct point_type<BoundingBox>
{
  using point_t = bg::model::point<double, 3, bg::cs::cartesian>; //fake point type.
  using type = point_t; //BoundingBox does not consist of this point type, actually.
};

template <std::size_t Dimension>
struct indexed_access<BoundingBox, min_corner, Dimension>
{
  static inline double get(const BoundingBox& bb)
  {
    if (Dimension >= bb.size())
        return std::numeric_limits<double>::lowest();
    return bb[Dimension].first;
  }
  static inline void set(BoundingBox& bb, double value)
  {
    if (Dimension >= bb.size())
        return;
    bb[Dimension].first = value;
  }
};

template <std::size_t Dimension>
struct indexed_access<BoundingBox, max_corner, Dimension>
{
  static inline double get(const BoundingBox& bb)
  {
    if (Dimension >= bb.size())
        return std::numeric_limits<double>::max();
    return bb[Dimension].second;
  }
  static inline void set(BoundingBox& bb, const double& value)
  {
    if (Dimension >= bb.size())
        return;
    bb[Dimension].second = value;
  }
};

BOOST_CONCEPT_ASSERT( (bg::concepts::Box<BoundingBox>) );

}}}

namespace precice {
namespace mesh {
namespace impl {

/// Makes a utils::PtrVector indexable and thus be usable in boost::geometry::rtree
template <typename Container>
class PtrVectorIndexable
{
  using size_type = typename Container::container::size_type;
  using cref = const typename Container::value_type&;
  Container const& container;

public:
  using result_type = cref;

  explicit PtrVectorIndexable(Container const& c) : container(c)
  {}

  result_type operator()(size_type i) const
  {
    return container[i];
  }
};


}}}

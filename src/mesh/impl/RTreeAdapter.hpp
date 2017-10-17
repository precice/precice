#pragma once

#include <boost/geometry.hpp>
#include <Eigen/Core>
#include "mesh/Vertex.hpp"

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
    p[Dimension] = value;
  }
};

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

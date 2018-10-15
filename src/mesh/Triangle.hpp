#pragma once

#include <Eigen/Core>
#include <array>
#include <boost/noncopyable.hpp>
#include "mesh/Edge.hpp"
#include "mesh/PropertyContainer.hpp"
#include "mesh/RangeAccessor.hpp"
#include "utils/assertion.hpp"

namespace precice
{
namespace mesh
{
class Vertex;
}
} // namespace precice

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice
{
namespace mesh
{

/// Triangle of a mesh, defined by three edges (and vertices).
class Triangle : public PropertyContainer, private boost::noncopyable
{
public:
  /// Type of the read-only const random-access iterator over Vertex coords
  /**
   * This index-based iterator iterates over the vertices of this Quad.
   * The returned value is the forwarded result of Vertex::getCoords.
   * It is thus a read-only random-access iterator.
   */
  using const_iterator = IndexRangeIterator<const Triangle, const Eigen::VectorXd>;

  /// Type of the read-only random access vertex iterator
  using iterator = const_iterator;


  /// Constructor, the order of edges defines the outer normal direction.
  Triangle(
      Edge &edgeOne,
      Edge &edgeTwo,
      Edge &edgeThree,
      int   id);

  /// Destructor, empty.
  virtual ~Triangle() {}

  /// Returns dimensionalty of space the triangle is embedded in.
  int getDimensions() const;

  /**
   * @brief Returns triangle vertex with index 0, 1 or 2.
   *
   * Vertex 0 is the first vertex of edge 0. Vertex 1 the second vertex of
   * edge 0. Vertex 2 is either the first or second vertex of edge 1, which
   * is determined on construction of the triangle.
   */
  Vertex &vertex(int i);

  /**
   * @brief Returns const triangle vertex with index 0, 1 or 2.
   *
   * Vertex 0 is the first vertex of edge 0. Vertex 1 the second vertex of
   * edge 0. Vertex 2 is either the first or second vertex of edge 1, which
   * is determined on construction of the triangle.
   */
  const Vertex &vertex(int i) const;

  /// Returns triangle edge with index 0, 1 or 2.
  Edge &edge(int i);

  /// Returns const triangle edge with index 0, 1 or 2.
  const Edge &edge(int i) const;

  ///@name Iterators
  ///@{

  /// Returns a read-only random-access iterator to the begin (0) of the vertex range [0,1,2]
  iterator begin();

  /// Returns a read-only random-access iterator to the end (3) of the vertex range [0,1,2]
  iterator end();

  /// Returns a read-only random-access iterator to the begin (0) of the vertex range [0,1,2]
  const_iterator begin() const;

  /// Returns a read-only random access iterator to the end (3) of the vertex range [0,1,2]
  const_iterator end() const;

  /// Returns a read-only random-access iterator to the begin (0) of the vertex range [0,1,2]
  const_iterator cbegin() const;

  /// Returns a read-only random access iterator to the end (3) of the vertex range [0,1,2]
  const_iterator cend() const;

  ///@}

  /// Sets the outer normal of the triangle.
  template <typename VECTOR_T>
  void setNormal(const VECTOR_T &normal);

  /// Sets the barycenter of the triangle.
  template <typename VECTOR_T>
  void setCenter(const VECTOR_T &center);

  /// Sets the radius of the circle enclosing the triangle.
  void setEnclosingRadius(double radius);

  /// Returns a among triangles globally unique ID.
  int getID() const;

  /**
   * @brief Returns the outer normal of the triangle.
   *
   * @pre The normal has to be computed and set from outside before.
   */
  const Eigen::VectorXd &getNormal() const;

  /**
   * @brief Returns the barycenter of the triangle.
   *
   * @pre The center has to be computed and set from outside before.
   */
  const Eigen::VectorXd &getCenter() const;

  /**
   * @brief Returns the radius of the circle enclosing the triangle.
   *
   * @pre The radius has to be computed and set from outside before.
   */
  double getEnclosingRadius() const;

private:
  /// Edges defining the triangle.
  std::array<Edge *, 3> _edges;

  /// Decider for choosing unique vertices from _edges.
  std::array<int, 3> _vertexMap;
  //  bool _vertexDeciderFirst;

  /// ID of the edge.
  int _id;

  /// Normal vector of the triangle.
  Eigen::VectorXd _normal;

  /// Center point of the triangle.
  Eigen::VectorXd _center;

  /// Minimal radius of circle enclosing the triangle.
  double _enclosingRadius = 0;
};

// --------------------------------------------------------- HEADER DEFINITIONS

inline Vertex &Triangle::vertex(int i)
{
  assertion((i >= 0) && (i < 3), i);
  return edge(i).vertex(_vertexMap[i]);
}

inline const Vertex &Triangle::vertex(int i) const
{
  assertion((i >= 0) && (i < 3), i);
  return edge(i).vertex(_vertexMap[i]);
}

inline Edge &Triangle::edge(int i)
{
  return *_edges[i];
}

inline const Edge &Triangle::edge(int i) const
{
  return *_edges[i];
}

inline Triangle::iterator Triangle::begin()
{
  return {this, 0};
}

inline Triangle::iterator Triangle::end()
{
  return {this, 3};
}

inline Triangle::const_iterator Triangle::begin() const
{
  return {this, 0};
}

inline Triangle::const_iterator Triangle::end() const
{
  return {this, 3};
}

inline Triangle::const_iterator Triangle::cbegin() const
{
  return begin();
}

inline Triangle::const_iterator Triangle::cend() const
{
  return end();
}

template <typename VECTOR_T>
void Triangle::setNormal(
    const VECTOR_T &normal)
{
  assertion(normal.size() == getDimensions(), normal.size(), getDimensions());
  _normal = normal;
}

template <typename VECTOR_T>
void Triangle::setCenter(
    const VECTOR_T &center)
{
  assertion(center.size() == getDimensions(), center.size(), getDimensions());
  _center = center;
}

inline int Triangle::getID() const
{
  return _id;
}

} // namespace mesh
} // namespace precice

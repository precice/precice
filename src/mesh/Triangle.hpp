#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <iostream>

#include "math/differences.hpp"
#include "mesh/Edge.hpp"
#include "mesh/RangeAccessor.hpp"
#include "precice/types.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace mesh {
class Vertex;
}
} // namespace precice

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace mesh {

/// Triangle of a mesh, defined by three edges (and vertices).
class Triangle {
public:
  /// Type of the read-only const random-access iterator over Vertex coords
  /**
   * This index-based iterator iterates over the vertices of this Triangle.
   * The returned value is the forwarded result of Vertex::getCoords.
   * It is thus a read-only random-access iterator.
   */
  using const_iterator = IndexRangeIterator<const Triangle, const Vertex::RawCoords>;

  /// Type of the read-only random access vertex iterator
  using iterator = const_iterator;

  /// Fix for the Boost.Test versions 1.65.1 - 1.67
  using value_type = Vertex::RawCoords;

  /// Constructor, the order of edges defines the outer normal direction.
  Triangle(
      Edge &     edgeOne,
      Edge &     edgeTwo,
      Edge &     edgeThree,
      TriangleID id);

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

  /// Computes the normal of the triangle.
  Eigen::VectorXd computeNormal() const;

  /// Returns a among triangles globally unique ID.
  TriangleID getID() const;

  /// Returns the surface area of the triangle
  double getArea() const;

  /// Returns the barycenter of the triangle.
  const Eigen::VectorXd getCenter() const;

  /// Returns the radius of the circle enclosing the triangle.
  double getEnclosingRadius() const;

  /**
   * @brief Compares two Triangles for equality
   *
   * Two Triangles are equal if their normal vector is equal AND
   * if the three edges are equal, whereas the order of edges is NOT important.
   */
  bool operator==(const Triangle &other) const;

  /// Not equal, implemented in terms of equal.
  bool operator!=(const Triangle &other) const;

private:
  /// Edges defining the triangle.
  std::array<Edge *, 3> _edges;

  /// Decider for choosing unique vertices from _edges.
  std::array<bool, 3> _vertexMap;

  /// ID of the triangle.
  TriangleID _id;
};

// --------------------------------------------------------- HEADER DEFINITIONS

inline Vertex &Triangle::vertex(int i)
{
  PRECICE_ASSERT((i >= 0) && (i < 3), i);
  return edge(i).vertex(_vertexMap[i]);
}

inline const Vertex &Triangle::vertex(int i) const
{
  PRECICE_ASSERT((i >= 0) && (i < 3), i);
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

inline TriangleID Triangle::getID() const
{
  return _id;
}

std::ostream &operator<<(std::ostream &os, const Triangle &t);

} // namespace mesh
} // namespace precice

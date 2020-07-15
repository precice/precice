#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <iostream>
#include "mesh/Edge.hpp"
#include "mesh/RangeAccessor.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace mesh {
class Vertex;

/// Quadrilateral (or Quadrangle) geometric primitive.
class Quad {
public:
  /// Type of the read-only const random-access iterator over Vertex coords
  /**
   * This index-based iterator iterates over the vertices of this Quad.
   * The returned value is the forwarded result of Vertex::getCoords.
   * It is thus a read-only random-access iterator.
   */
  using const_iterator = IndexRangeIterator<const Quad, const Eigen::VectorXd>;

  /// Type of the random access vertex iterator
  using iterator = const_iterator; //IndexRangeIterator<Quad, Eigen::Vector3d>;

  /// Fix for the Boost.Test versions 1.65 - 1.67.1
  using value_type = Eigen::VectorXd;

  /// Constructor, the order of edges defines the outer normal direction.
  Quad(
      Edge &edgeOne,
      Edge &edgeTwo,
      Edge &edgeThree,
      Edge &edgeFour,
      int   id);

  /// Returns dimensionalty of space the quad is embedded in.
  int getDimensions() const;

  /**
   * @brief Returns quad vertex with index 0, 1, 2, or 3.
   *
   * Vertex 0 is the first vertex of edge 0. Vertex 1 the second vertex of
   * edge 0. Vertex 2 is either the first or second vertex of edge 1, which
   * is determined on construction of the triangle.
   */
  Vertex &vertex(int i);

  /**
   * @brief Returns const quad vertex with index 0, 1, 2, or 3.
   *
   * Vertex 0 is the first vertex of edge 0. Vertex 1 the second vertex of
   * edge 0. Vertex 2 is either the first or second vertex of edge 1, which
   * is determined on construction of the triangle.
   */
  const Vertex &vertex(int i) const;

  /// Returns quad edge with index 0, 1, 2, or 3.
  Edge &edge(int i);

  /// Returns const quad edge with index 0, 1, 2, or 3.
  const Edge &edge(int i) const;

  ///@name Iterators
  ///@{

  /// Returns a read-only random-access iterator to the begin (0) of the vertex range [0,1,2,3]
  iterator begin();

  /// Returns a read-only random-access iterator to the end (4) of the vertex range [0,1,2,3]
  iterator end();

  /// Returns a read-only random-access iterator to the begin (0) of the vertex range [0,1,2,3]
  const_iterator begin() const;

  /// Returns a read-only random-access iterator to the end (4) of the vertex range [0,1,2,3]
  const_iterator end() const;

  /// Returns a read-only random-access iterator to the begin (0) of the vertex range [0,1,2,3]
  const_iterator cbegin() const;

  /// Returns a read-only random-access iterator to the end (4) of the vertex range [0,1,2,3]
  const_iterator cend() const;

  ///@}

  /// Sets the outer normal of the quad.
  template <typename VECTOR_T>
  void setNormal(const VECTOR_T &normal);

  /// Computes and sets the normal of the triangle, returns the area-weighted normal.
  const Eigen::VectorXd computeNormal(bool flip = false);

  /// Returns a among quads globally unique ID.
  int getID() const;

  /**
   * @brief Returns the outer normal of the quad.
   *
   * @pre The normal has to be computed and set from outside before.
   */
  const Eigen::VectorXd &getNormal() const;

  /// Returns the barycenter of the quad.
  const Eigen::VectorXd getCenter() const;

  /// Returns the radius of the circle enclosing the quad.
  double getEnclosingRadius() const;

  /**
   * @brief Compares two Quads for equality
   *
   * Two Quads are equal if their normal vector is equal AND
   * if the four edges are equal, whereas the order of edges is NOT important.
   */
  bool operator==(const Quad &other) const;

  /// Not equal, implemented in terms of equal.
  bool operator!=(const Quad &other) const;

private:
  /// Edges defining the quad.
  std::array<Edge *, 4> _edges;

  /// Decider for choosing unique vertices from _edges.
  std::array<int, 4> _vertexMap;

  /// ID of the edge.
  int _id;

  /// Normal vector of the quad.
  Eigen::VectorXd _normal;
};

// --------------------------------------------------------- HEADER DEFINITIONS

inline Vertex &Quad::vertex(int i)
{
  PRECICE_ASSERT((i >= 0) && (i < 4), i);
  return edge(i).vertex(_vertexMap[i]);
}

inline const Vertex &Quad::vertex(int i) const
{
  PRECICE_ASSERT((i >= 0) && (i < 4), i);
  return edge(i).vertex(_vertexMap[i]);
}

inline Quad::iterator Quad::begin()
{
  return {this, 0};
}

inline Quad::iterator Quad::end()
{
  return {this, 4};
}

inline Quad::const_iterator Quad::begin() const
{
  return {this, 0};
}

inline Quad::const_iterator Quad::end() const
{
  return {this, 4};
}

inline Quad::const_iterator Quad::cbegin() const
{
  return begin();
}

inline Quad::const_iterator Quad::cend() const
{
  return end();
}

inline Edge &Quad::edge(int i)
{
  return *_edges[i];
}

inline const Edge &Quad::edge(int i) const
{
  return *_edges[i];
}

template <typename VECTOR_T>
void Quad::setNormal(const VECTOR_T &normal)
{
  PRECICE_ASSERT(normal.size() == getDimensions(), normal.size(), getDimensions());
  _normal = normal;
}

inline int Quad::getID() const
{
  return _id;
}

std::ostream &operator<<(std::ostream &os, const Quad &q);

} // namespace mesh
} // namespace precice

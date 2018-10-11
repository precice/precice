#pragma once

#include <iostream>
#include <algorithm>
#include "boost/noncopyable.hpp"
#include "mesh/Edge.hpp"
#include "mesh/PropertyContainer.hpp"
#include "math/differences.hpp"

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

/// Quadrilateral (or Quadrangle) geometric primitive.
class Quad : public PropertyContainer, private boost::noncopyable
{
public:
  /// Constructor, the order of edges defines the outer normal direction.
  Quad(
      Edge &edgeOne,
      Edge &edgeTwo,
      Edge &edgeThree,
      Edge &edgeFour,
      int   id);

  /// Destructor, empty.
  virtual ~Quad() {}

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

  /// Sets the outer normal of the quad.
  template <typename VECTOR_T>
  void setNormal(const VECTOR_T &normal);

  /// Sets the center of the quad.
  template <typename VECTOR_T>
  void setCenter(const VECTOR_T &center);

  /// Sets the radius of the circle enclosing the quad.
  void setEnclosingRadius(double radius);

  /// Returns a among quads globally unique ID.
  int getID() const;

  /**
   * @brief Returns the outer normal of the quad.
   *
   * Prerequesits: The normal has to be computed and set from outside before.
   */
  const Eigen::VectorXd &getNormal() const;

  /**
   * @brief Returns the barycenter of the quad.
   *
   * Prerequesits: The center has to be computed and set from outside before.
   */
  const Eigen::VectorXd &getCenter() const;

  /**
   * @brief Returns the radius of the circle enclosing the quad.
   *
   * Prerequesits: The radius has to be computed and set from outside before.
   */
  double getEnclosingRadius() const;

  friend std::ostream& operator<<(std::ostream& os, const Quad& q);

  inline bool operator==(const Quad& other) const;

  inline bool operator!=(const Quad& other) const;

private:
  /// Edges defining the quad.
  std::array<Edge *, 4> _edges;

  /// Decider for choosing unique vertices from _edges.
  std::array<int, 4> _vertexMap;

  /// ID of the edge.
  int _id;

  /// Normal vector of the quad.
  Eigen::VectorXd _normal;

  /// Center point of the quad.
  Eigen::VectorXd _center;

  /// Minimal radius of circle enclosing the quad.
  double _enclosingRadius = 0;
};

// --------------------------------------------------------- HEADER DEFINITIONS

inline Vertex &Quad::vertex(int i)
{
  assertion((i >= 0) && (i < 4), i);
  return edge(i).vertex(_vertexMap[i]);
}

inline const Vertex &Quad::vertex(int i) const
{
  assertion((i >= 0) && (i < 4), i);
  return edge(i).vertex(_vertexMap[i]);
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
  assertion(normal.size() == getDimensions(), normal.size(), getDimensions());
  _normal = normal;
}

template <typename VECTOR_T>
void Quad::setCenter(const VECTOR_T &center)
{
  assertion(center.size() == getDimensions(), center.size(), getDimensions());
  _center = center;
}

inline int Quad::getID() const
{
  return _id;
}

inline std::ostream& operator<<(std::ostream& os, const Quad& q)
{
  return os << "Quad " << q._id << " defined by:\n"
      << "\t" << *q._edges[0] << "\t" << *q._edges[1] 
      << "\t" << *q._edges[2] << "\t" << *q._edges[3];
}

inline bool Quad::operator==(const Quad& other) const
{
    return math::equals(_normal, other._normal) &&
        std::is_permutation(_edges.begin(), _edges.end(), other._edges.begin(), 
                [](const Edge* e1, const Edge* e2){return *e1 == *e2;});
}

inline bool Quad::operator!=(const Quad& other) const
{
  return !(*this == other);
}
} // namespace mesh
} // namespace precice

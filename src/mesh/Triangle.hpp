#pragma once

#include <Eigen/Core>
#include <boost/noncopyable.hpp>
#include <iostream>
#include <algorithm>

#include "mesh/Edge.hpp"
#include "mesh/PropertyContainer.hpp"
#include "utils/assertion.hpp"
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

/// Triangle of a mesh, defined by three edges (and vertices).
class Triangle : public PropertyContainer, private boost::noncopyable
{
public:
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

  friend std::ostream& operator<<(std::ostream& os, const Triangle& t);

  inline bool operator==(const Triangle& other) const;

  inline bool operator!=(const Triangle& other) const;

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

inline std::ostream& operator<<(std::ostream& os, const Triangle& t)
{
  return os << "Triangle " << t._id << " defined by:\n" 
      << "\t" << *t._edges[0] << "\t" << *t._edges[1] 
      << "\t" << *t._edges[2];
}

inline bool Triangle::operator==(const Triangle& other) const
{
    return math::equals(_normal, other._normal) &&
        std::is_permutation(_edges.begin(), _edges.end(), other._edges.begin(), 
                [](const Edge* e1, const Edge* e2){return *e1 == *e2;});
}
inline bool Triangle::operator!=(const Triangle& other) const
{
    return !(*this == other);
}

} // namespace mesh
} // namespace precice

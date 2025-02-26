#pragma once

#include <Eigen/Core>
#include <array>
#include <iostream>

#include "math/differences.hpp"
#include "mesh/Vertex.hpp"
#include "precice/impl/Types.hpp"
#include "utils/assertion.hpp"

namespace precice::mesh {

/// Linear edge of a mesh, defined by two Vertex objects.
class Edge {
public:
  /// Amount of vertices
  static constexpr int vertexCount{2};

  /**
   * @brief Constructor.
   *
   * @param[in] vertexOne First Vertex object defining the edge.
   * @param[in] vertexTwo Second Vertex object defining the edge.
   */
  Edge(
      Vertex &vertexOne,
      Vertex &vertexTwo);

  /// Returns number of spatial dimensions (2 or 3) the edge is embedded to.
  int getDimensions() const;

  /// Returns the edge's vertex with index 0 or 1.
  Vertex &vertex(int i);

  /// Returns the edge's vertex as const object with index 0 or 1.
  const Vertex &vertex(int i) const;

  /// Returns the length of the edge
  double getLength() const;

  /// Returns the center of the edge.
  const Eigen::VectorXd getCenter() const;

  /// Returns the radius of the enclosing circle of the edge.
  double getEnclosingRadius() const;

  /// Checks whether both edges share a vertex.
  bool connectedTo(const Edge &other) const;

  /**
   * @brief Compares two Edges for equality
   *
   * Two Edges are equal if the two vertices are equal,
   * whereas the order of vertices is NOT important.
   */
  bool operator==(const Edge &other) const;

  /// Not equal, implemented in terms of equal.
  bool operator!=(const Edge &other) const;

  /// Weak ordering based on vertex ids
  bool operator<(const Edge &other) const
  {
    return std::make_tuple(_vertices[0]->getID(), _vertices[1]->getID()) <
           std::make_tuple(other._vertices[0]->getID(), other._vertices[1]->getID());
  }

private:
  /// Pointers to Vertex objects defining the edge, ordered by Vertex::getID().
  std::array<Vertex *, 2> _vertices;
};

// ------------------------------------------------------ HEADER IMPLEMENTATION

inline Vertex &Edge::vertex(
    int i)
{
  PRECICE_ASSERT((i == 0) || (i == 1), i);
  return *_vertices[i];
}

inline const Vertex &Edge::vertex(
    int i) const
{
  PRECICE_ASSERT((i == 0) || (i == 1), i);
  return *_vertices[i];
}

inline int Edge::getDimensions() const
{
  return _vertices[0]->getDimensions();
}

std::ostream &operator<<(std::ostream &stream, const Edge &edge);

} // namespace precice::mesh

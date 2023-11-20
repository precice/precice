#pragma once

#include <array>
#include <iostream>
#include <tuple>

#include "mesh/Vertex.hpp"
#include "precice/types.hpp"
#include "utils/assertion.hpp"

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace mesh {

/// Tetrahedron of a mesh, defined by 4 vertices
class Tetrahedron {
public:
  /// Amount of vertices
  static constexpr int vertexCount{4};

  /** Constructor based on 4 vertices
   *
   * The vertices will be sorted by Vertex::getID().
   * This allows to weakly order tetrahedra.
   */
  Tetrahedron(
      Vertex &vertexOne,
      Vertex &vertexTwo,
      Vertex &vertexThree,
      Vertex &vertexFour);

  /// Returns dimensionalty of space the Tetrahedron is embedded in.
  int getDimensions() const;

  /**
   * @brief Returns tetrahedron vertex with index 0, 1, 2 or 3.
   */
  Vertex &vertex(int i);

  /**
   * @brief Returns const tetrahedron vertex with index 0, 1, 2 or 3.
   */
  const Vertex &vertex(int i) const;

  /// Returns the unsigned volume of the tetrahedron
  double getVolume() const;

  /// Returns the barycenter of the tetrahedron.
  const Eigen::VectorXd getCenter() const;

  /// Returns the radius of the sphere enclosing the tetrahedron.
  double getEnclosingRadius() const;

  /**
   * @brief Compares two Tetrahedrons for equality
   *
   * Two Tetrahedrons are equal if their vertices are the same, up to permutations
   */
  bool operator==(const Tetrahedron &other) const;

  /// Not equal, implemented in terms of equal.
  bool operator!=(const Tetrahedron &other) const;

  /// Weak ordering based on vertex ids
  bool operator<(const Tetrahedron &other) const
  {
    return std::make_tuple(_vertices[0]->getID(), _vertices[1]->getID(), _vertices[2]->getID(), _vertices[3]->getID()) <
           std::make_tuple(other._vertices[0]->getID(), other._vertices[1]->getID(), other._vertices[2]->getID(), other._vertices[3]->getID());
  }

private:
  /// Vertices defining the Tetrahedron.
  std::array<Vertex *, 4> _vertices;
};

// --------------------------------------------------------- HEADER DEFINITIONS

inline Vertex &Tetrahedron::vertex(int i)
{
  PRECICE_ASSERT((i >= 0) && (i < 4), i);
  return *_vertices[i];
}

inline const Vertex &Tetrahedron::vertex(int i) const
{
  PRECICE_ASSERT((i >= 0) && (i < 4), i);
  return *_vertices[i];
}

std::ostream &operator<<(std::ostream &os, const Tetrahedron &t);

} // namespace mesh
} // namespace precice

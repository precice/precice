#pragma once

#include <Eigen/Core>
#include <array>
#include <iostream>
#include <utility>

#include "math/differences.hpp"
#include "precice/types.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace mesh {

/// Vertex of a mesh.
class Vertex {
public:
  //( Used as the raw representation of the coordinates
  using RawCoords = std::array<double, 3>;

  /// Constructor for vertex
  template <typename VECTOR_T>
  Vertex(
      const VECTOR_T &coordinates,
      VertexID        id);

  /// Returns spatial dimensionality of vertex.
  int getDimensions() const;

  /// Sets the coordinates of the vertex.
  template <typename VECTOR_T>
  void setCoords(const VECTOR_T &coordinates);

  /// Returns the unique (among vertices of one mesh on one processor) ID of the vertex.
  VertexID getID() const;

  /// Returns the coordinates of the vertex.
  Eigen::VectorXd getCoords() const;

  /// Direct access to the coordinates
  const RawCoords &rawCoords() const;

  /// Globally unique index
  int getGlobalIndex() const;

  void setGlobalIndex(int globalIndex);

  bool isOwner() const;

  void setOwner(bool owner);

  bool isTagged() const;

  void tag();

  inline bool operator==(const Vertex &rhs) const;

  inline bool operator!=(const Vertex &rhs) const;

private:
  /// Coordinates of the vertex
  std::array<double, 3> _coords;

  /// Dimension of the coordinates. 3D or 2D
  short _dim;

  /// Unique (among vertices in one mesh) ID of the vertex.
  VertexID _id;

  /// global (unique) index for parallel simulations
  int _globalIndex = -1;

  /// true if this processors is the owner of the vertex (for parallel simulations)
  bool _owner = true;

  /// true if this vertex is tagged for partition
  bool _tagged = false;
};

// ------------------------------------------------------ HEADER IMPLEMENTATION

template <typename VECTOR_T>
Vertex::Vertex(
    const VECTOR_T &coordinates,
    int             id)
    : _dim(coordinates.size()),
      _id(id)
{
  PRECICE_ASSERT(_dim == 2 || _dim == 3, _dim);
  _coords[0] = coordinates[0];
  _coords[1] = coordinates[1];
  _coords[2] = (_dim == 3) ? coordinates[2] : 0.0;
}

template <typename VECTOR_T>
void Vertex::setCoords(
    const VECTOR_T &coordinates)
{
  PRECICE_ASSERT(coordinates.size() == _dim, coordinates.size(), _dim);
  _coords[0] = coordinates[0];
  _coords[1] = coordinates[1];
  _coords[2] = (_dim == 3) ? coordinates[2] : 0.0;
}

inline VertexID Vertex::getID() const
{
  return _id;
}

inline Eigen::VectorXd Vertex::getCoords() const
{
  Eigen::VectorXd v(_dim);
  std::copy_n(_coords.data(), _dim, v.data());
  return v;
}

inline const Vertex::RawCoords &Vertex::rawCoords() const
{
  return _coords;
}

inline bool Vertex::operator==(const Vertex &rhs) const
{
  return math::equals(getCoords(), rhs.getCoords());
}

inline bool Vertex::operator!=(const Vertex &rhs) const
{
  return !(*this == rhs);
}

/// Make Vertex printable
std::ostream &operator<<(std::ostream &os, const Vertex &v);

} // namespace mesh
} // namespace precice

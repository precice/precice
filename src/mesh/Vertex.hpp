#pragma once

#include <Eigen/Core>
#include <iostream>

#include "math/differences.hpp"

namespace precice {
namespace mesh {

/// Vertex of a mesh.
class Vertex {
public:
  /// Constructor for vertex
  template <typename VECTOR_T>
  Vertex(
      const VECTOR_T &coordinates,
      int             id);

  /// Returns spatial dimensionality of vertex
  int getDimensions() const;

  /// Sets the coordinates of the vertex.
  template <typename VECTOR_T>
  void setCoords(const VECTOR_T &coordinates);

  /// Sets the coordinates of the vertex, rvalue variant
  template <typename VECTOR_T>
  void setCoords(VECTOR_T &&coordinates);

  /// Sets the normal of the vertex.
  template <typename VECTOR_T>
  void setNormal(const VECTOR_T &normal);

  /// Sets the normal of the vertex, rvalue variant
  template <typename VECTOR_T>
  void setNormal(VECTOR_T &&normal);

  /// Returns the unique (among vertices of one mesh on one processor) ID of the vertex.
  int getID() const;

  /// Returns the patchID that a vertex belongs to
  int getPatchID() const;

  /// Sets the patchID of a vertex
  void setPatchID(int patchid);

  /// Sets the name of the patch of that vertex
  void setPatchName(std::string patchName);

  /// Returns the name of the patch of the vertex
  void getPatchName() const;

  /// Returns the coordinates of the vertex.
  const Eigen::VectorXd &getCoords() const;

  /// Returns the normal of the vertex.
  const Eigen::VectorXd &getNormal() const;

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
  /// Unique (among vertices in one mesh) ID of the vertex.
  int _id;

  /// Patch id that a vertex belongs to
  int _patchid = 0;

  std::string _patchname;

  /// Coordinates of the vertex.
  Eigen::VectorXd _coords;

  /// Normal of the vertex.
  Eigen::VectorXd _normal;

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
    : _id(id),
      _coords(coordinates),
      _normal(Eigen::VectorXd::Constant(_coords.size(), 0.0))
{
}

template <typename VECTOR_T>
void Vertex::setCoords(
    const VECTOR_T &coordinates)
{
  PRECICE_ASSERT(coordinates.size() == _coords.size(), coordinates.size(), _coords.size());
  _coords = coordinates;
}

template <typename VECTOR_T>
void Vertex::setCoords(
    VECTOR_T &&coordinates)
{
  PRECICE_ASSERT(coordinates.size() == _coords.size(), coordinates.size(), _coords.size());
  _coords = std::forward<VECTOR_T>(coordinates);
}

template <typename VECTOR_T>
void Vertex::setNormal(
    const VECTOR_T &normal)
{
  PRECICE_ASSERT(normal.size() == _normal.size(), normal.size(), _normal.size());
  _normal = normal;
}

template <typename VECTOR_T>
void Vertex::setNormal(
    VECTOR_T &&normal)
{
  PRECICE_ASSERT(normal.size() == _normal.size(), normal.size(), _normal.size());
  _normal = std::forward<VECTOR_T>(normal);
}

inline int Vertex::getID() const
{
  return _id;
}

inline int Vertex::getPatchID() const
{
 return _patchid;
}

inline int Vertex::getPatchName() const
{
 return _patchname;
}

inline const Eigen::VectorXd &Vertex::getCoords() const
{
  return _coords;
}

inline bool Vertex::operator==(const Vertex &rhs) const
{
  return math::equals(_coords, rhs._coords);
}

inline bool Vertex::operator!=(const Vertex &rhs) const
{
  return !(*this == rhs);
}

/// Make Vertex printable
std::ostream &operator<<(std::ostream &os, const Vertex &v);

} // namespace mesh
} // namespace precice

#pragma once

#include "mesh/PropertyContainer.hpp"
#include "boost/noncopyable.hpp"
#include <Eigen/Core>

namespace precice {
  namespace mesh {
    class Mesh;
  }
}

namespace precice {
namespace mesh {

/// Vertex of a mesh.
class Vertex : public PropertyContainer, private boost::noncopyable
{
public:

  /**
   * @brief Constructor for vertex, parent mesh is not assigned.
   */
  template<typename VECTOR_T>
  Vertex (
    const VECTOR_T& coordinates,
    int             id );

  /**
   * @brief Constructor for vertex, parent mesh is assigned.
   */
  template<typename VECTOR_T>
  Vertex (
    const VECTOR_T& coordinates,
    int             id,
    Mesh&           mesh );

  /**
   * @brief Destructor, empty.
   */
  virtual ~Vertex() {}

  /**
   * @brief Returns spatial dimenionality of vertex.
   */
  int getDimensions() const;

  /**
   * @brief Sets the coordinates of the vertex.
   */
  template<typename VECTOR_T>
  void setCoords ( const VECTOR_T& coordinates );

  /**
   * @brief Sets the normal of the vertex.
   */
  template<typename VECTOR_T>
  void setNormal ( const VECTOR_T& normal );

  /**
   * @brief Returns the unique (among vertices of one mesh) ID of the vertex.
   */
  int getID() const;

  /// Returns the coordinates of the vertex.
  const Eigen::VectorXd& getCoords() const;

  /// Returns the normal of the vertex.
  const Eigen::VectorXd& getNormal() const;

  /**
   * @brief Returns (possibly NULL) pointer to parent const Mesh object.
   */
  const Mesh* mesh() const;

  /**
   * @brief Returns possibly NULL pointer to parent Mesh object.
   */
  Mesh* mesh();

  int getGlobalIndex() const;

  void setGlobalIndex(int globalIndex);

  bool isOwner() const;

  void setOwner(bool owner);

  bool isTagged() const;

  void tag();

private:

  // @brief Unique (among vertices in one mesh) ID of the vertex.
  int _id;

  /// Coordinates of the vertex.
  Eigen::VectorXd _coords;

  /// Normal of the vertex.
  Eigen::VectorXd _normal;

  // @brief global (unique) index for parallel simulations
  int _globalIndex;

  // @brief true if this processors is the owner of the vertex (for parallel simulations)
  bool _owner;

  // @brief true if this vertex is tagged for partition
  bool _tagged;

  // @brief Pointer to parent mesh, possibly NULL.
  Mesh * _mesh;
};

// ------------------------------------------------------ HEADER IMPLEMENTATION

template<typename VECTOR_T>
Vertex:: Vertex
(
  const VECTOR_T& coordinates,
  int             id )
:
  PropertyContainer (),
  _id ( id ),
  _coords ( coordinates ),
  _normal ( Eigen::VectorXd::Constant(_coords.size(), 0.0) ),
  _globalIndex(-1),
  _owner(true),
  _tagged(false),
  _mesh ( NULL )
{}

template<typename VECTOR_T>
Vertex:: Vertex (
  const VECTOR_T& coordinates,
  int             id,
  Mesh&           mesh )
:
  PropertyContainer (),
  _id ( id ),
  _coords ( coordinates ),
  _normal ( Eigen::VectorXd::Constant(_coords.size(), 0.0) ),
  _globalIndex(-1),
  _owner(true),
  _tagged(false),
  _mesh ( & mesh )
{}

template<typename VECTOR_T>
void Vertex:: setCoords
(
  const VECTOR_T& coordinates )
{
  assertion ( coordinates.size() == _coords.size(), coordinates.size(), _coords.size() );
  _coords = coordinates;
}

template<typename VECTOR_T>
void Vertex:: setNormal
(
  const VECTOR_T& normal )
{
  assertion ( normal.size() == _normal.size(), normal.size(), _normal.size() );
  _normal = normal;
}

inline int Vertex:: getID() const
{
  return _id;
}

inline const Eigen::VectorXd& Vertex::getCoords() const
{
  return _coords;
}

}} // namespace precice, mesh


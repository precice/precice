// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_MESH_VERTEX_HPP_
#define PRECICE_MESH_VERTEX_HPP_

#include "mesh/PropertyContainer.hpp"
#include "utils/Dimensions.hpp"
#include "boost/noncopyable.hpp"
#include <map>

namespace precice {
  namespace mesh {
    class Mesh;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace mesh {

/**
 * @brief Vertex of a mesh.
 */
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

  /**
   * @brief Returns the coordinates of the vertex.
   */
  const utils::DynVector& getCoords() const;

  /**
   * @brief Returns the normal of the vertex.
   */
  const utils::DynVector& getNormal() const;

  /**
   * @brief Returns (possibly NULL) pointer to parent const Mesh object.
   */
  const Mesh* mesh() const;

  /**
   * @brief Returns possibly NULL pointer to parent Mesh object.
   */
  Mesh* mesh();

  int getGlobalIndex();

  void setGlobalIndex(int globalIndex);

  bool isOwner();

  void setOwner(bool owner);

private:

  // @brief Unique (among vertices in one mesh) ID of the vertex.
  int _id;

  // @brief Coordinates of the vertex.
  utils::DynVector _coords;

  // @brief Normal of the vertex.
  utils::DynVector _normal;

  // @brief global (unique) index for parallel simulations
  int _globalIndex;

  // @brief true if this processors is the owner of the vertex (for parallel simulations)
  bool _owner;

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
  _normal ( _coords.size(), 0.0 ),
  _globalIndex(-1),
  _owner(true),
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
  _normal ( _coords.size(), 0.0 ),
  _globalIndex(-1),
  _owner(true),
  _mesh ( & mesh )
{}

template<typename VECTOR_T>
void Vertex:: setCoords
(
  const VECTOR_T& coordinates )
{
  assertion2 ( coordinates.size() == _coords.size(), coordinates.size(), _coords.size() );
  _coords = coordinates;
}

template<typename VECTOR_T>
void Vertex:: setNormal
(
  const VECTOR_T& normal )
{
  assertion2 ( normal.size() == _normal.size(), normal.size(), _normal.size() );
  _normal = normal;
}

inline int Vertex:: getID() const
{
  return _id;
}

inline const utils::DynVector& Vertex:: getCoords() const
{
  return _coords;
}

}} // namespace precice, mesh

#endif /* PRECICE_MESH_VERTEX_HPP_ */

// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_MESH_EDGE_HPP_
#define PRECICE_MESH_EDGE_HPP_

#include "mesh/PropertyContainer.hpp"
#include "mesh/Vertex.hpp"
#include "boost/array.hpp"
#include "boost/noncopyable.hpp"

namespace precice {
namespace mesh {

struct ConstEdgeIteratorTypes;
struct EdgeIteratorTypes;
template<typename Types> class EdgeIterator;

/**
 * @brief Linear edge of a mesh, defined by two Vertex objects.
 */
class Edge : public PropertyContainer, private boost::noncopyable
{
public:

  /**
   * @brief Constructor.
   *
   * @param vertexOne [IN] First Vertex object defining the edge.
   * @param vertexOne [IN] Second Vertex object defining the edge.
   * @param id [IN] Unique (among edges in one mesh) ID.
   */
  Edge (
    Vertex& vertexOne,
    Vertex& vertexTwo,
    int     id );

  /**
   * @brief Destructor, empty.
   */
  virtual ~Edge () {}

  /**
   * @brief Returns number of spatial dimensions (2 or 3) the edge is embedded to.
   */
  int getDimensions() const;

  /**
   * @brief Returns the edge's vertex with index 0 or 1.
   */
  Vertex& vertex ( int i );

  /**
   * @brief Returns the edge's vertex as const object with index 0 or 1.
   */
  const Vertex& vertex ( int i ) const;

  /**
   * @brief Sets the normal of the edge.
   */
  template<typename VECTOR_T>
  void setNormal ( const VECTOR_T& normal );

  /**
   * @brief Sets the center of the edge.
   */
  template<typename VECTOR_T>
  void setCenter ( const VECTOR_T& center );

  /**
   * @brief Sets the radius of the circle enclosing the edge.
   */
  void setEnclosingRadius ( double radius );

  /**
   * @brief Returns the (among edges) unique ID of the edge.
   */
  int getID () const;

  /**
   * @brief Returns the normal of the edge.
   */
  const utils::DynVector& getNormal () const;

  /**
   * @brief Returns the center of the edge.
   */
  const utils::DynVector& getCenter () const;

  /**
   * @brief Returns the radius of the enclosing circle of the edge.
   */
  double getEnclosingRadius () const;

private:

  // @brief Pointers to Vertex objects defining the edge.
  boost::array<Vertex*,2> _vertices;

  // @brief Unique (among edges) ID of the edge.
  int _id;

  // @brief Normal of the edge.
  utils::DynVector _normal;

  // @brief Center of the edge.
  utils::DynVector _center;

  // @brief Radius of the enclosing circle.
  double _enclosingRadius;
};

// ------------------------------------------------------ HEADER IMPLEMENTATION

inline Vertex& Edge:: vertex
(
  int i )
{
  assertion1 ( (i == 0) || (i == 1), i );
  return *_vertices[i];
}

inline const Vertex& Edge:: vertex
(
  int i ) const
{
  assertion1 ( (i==0) || (i==1), i );
  return *_vertices[i];
}

inline int Edge:: getDimensions() const
{
  return _vertices[0]->getDimensions();
}

template<typename VECTOR_T>
void Edge:: setNormal
(
  const VECTOR_T& normal )
{
  assertion2 ( normal.size() == _vertices[0]->getDimensions(), normal,
               _vertices[0]->getDimensions() );
  _normal = normal;
}

template<typename VECTOR_T>
void Edge:: setCenter
(
  const VECTOR_T& center )
{
  assertion2 ( center.size() == _vertices[0]->getDimensions(), center,
               _vertices[0]->getDimensions() );
  _center = center;
}

inline const utils::DynVector& Edge:: getNormal () const
{
  return _normal;
}

}} // namespace precice, mesh

#endif /* PRECICE_MESH_EDGE_HPP_ */

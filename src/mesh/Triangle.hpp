// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_MESH_TRIANGLE_HPP_
#define PRECICE_MESH_TRIANGLE_HPP_

#include "mesh/PropertyContainer.hpp"
#include "mesh/Edge.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Helpers.hpp"
#include "boost/noncopyable.hpp"
#include "boost/array.hpp"

namespace precice {
  namespace mesh {
    class Vertex;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace mesh {

/**
 * @brief Triangle of a mesh, defined by three edges (and vertices).
 */
class Triangle : public PropertyContainer, private boost::noncopyable
{
public:

  /**
   * @brief Constructor, the order of edges defines the outer normal direction.
   */
  Triangle (
    Edge& edgeOne,
    Edge& edgeTwo,
    Edge& edgeThree,
    int   id );

  /**
   * @brief Destructor, empty.
   */
  virtual ~Triangle() {}

  /**
   * @brief Returns dimensionalty of space the triangle is embedded in.
   */
  int getDimensions() const;

  /**
   * @brief Returns triangle vertex with index 0, 1 or 2.
   *
   * Vertex 0 is the first vertex of edge 0. Vertex 1 the second vertex of
   * edge 0. Vertex 2 is either the first or second vertex of edge 1, which
   * is determined on construction of the triangle.
   */
  Vertex& vertex ( int i );

  /**
   * @brief Returns const triangle vertex with index 0, 1 or 2.
   *
   * Vertex 0 is the first vertex of edge 0. Vertex 1 the second vertex of
   * edge 0. Vertex 2 is either the first or second vertex of edge 1, which
   * is determined on construction of the triangle.
   */
  const Vertex& vertex ( int i ) const;

  /**
   * @brief Returns triangle edge with index 0, 1 or 2.
   */
  Edge& edge ( int i );

  /**
   * @brief Returns const triangle edge with index 0, 1 or 2.
   */
  const Edge& edge ( int i ) const;

  /**
   * @brief Sets the outer normal of the triangle.
   */
  template<typename VECTOR_T>
  void setNormal ( const VECTOR_T& normal );

  /**
   * @brief Sets the barycenter of the triangle.
   */
  template<typename VECTOR_T>
  void setCenter ( const VECTOR_T& center );

  /**
   * @brief Sets the radius of the circle enclosing the triangle.
   */
  void setEnclosingRadius ( double radius );

  /**
   * @brief Returns a among triangles globally unique ID.
   */
  int getID() const;

  /**
   * @brief Returns the outer normal of the triangle.
   *
   * Prerequesits: The normal has to be computed and set from outside before.
   */
  const utils::DynVector& getNormal () const;

  /**
   * @brief Returns the barycenter of the triangle.
   *
   * Prerequesits: The center has to be computed and set from outside before.
   */
  const utils::DynVector& getCenter () const;

  /**
   * @brief Returns the radius of the circle enclosing the triangle.
   *
   * Prerequesits: The radius has to be computed and set from outside before.
   */
  double getEnclosingRadius () const;

private:

  // @brief Edges defining the triangle.
  boost::array<Edge*,3> _edges;

  // @brief Decider for choosing unique vertices from _edges.
  boost::array<int,3> _vertexMap;
//  bool _vertexDeciderFirst;

  // @brief ID of the edge.
  int _id;

  // @brief Normal vector of the triangle.
  utils::DynVector _normal;

  // @brief Center point of the triangle.
  utils::DynVector _center;

  // @brief Minimal radius of circle enclosing the triangle.
  double _enclosingRadius;
};



// --------------------------------------------------------- HEADER DEFINITIONS

inline Vertex& Triangle:: vertex ( int i )
{
   assertion1 ( (i >= 0) && (i < 3), i );
   return edge(i).vertex(_vertexMap[i]);
};

inline const Vertex& Triangle:: vertex ( int i ) const
{
   assertion1 ( (i >= 0) && (i < 3), i );
   return edge(i).vertex(_vertexMap[i]);
};

inline Edge& Triangle:: edge ( int i )
{
   assertion1 ( (i >= 0) || (i < 3), i );
   return *_edges[i];
};

inline const Edge& Triangle:: edge ( int i ) const
{
   assertion1 ( (i >= 0) || (i < 3), 1 );
   return *_edges[i];
};

template<typename VECTOR_T>
void Triangle:: setNormal
(
  const VECTOR_T& normal )
{
  assertion2 ( normal.size() == getDimensions(), normal.size(), getDimensions() );
  _normal = normal;
}

template<typename VECTOR_T>
void Triangle:: setCenter
(
  const VECTOR_T& center )
{
  assertion2 ( center.size() == getDimensions(), center.size(), getDimensions() );
  _center = center;
}

inline int Triangle:: getID() const
{
  return _id;
}

}} // namespace precice, mesh

#endif /* PRECICE_MESH_TRIANGLE_HPP_ */

#ifndef PRECICE_MESH_QUAD_HPP_
#define PRECICE_MESH_QUAD_HPP_

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
 * @brief Quadrilateral (or Quadrangle) geometric primitive.
 */
class Quad : public PropertyContainer, private boost::noncopyable
{
public:

  /**
   * @brief Constructor, the order of edges defines the outer normal direction.
   */
  Quad (
    Edge& edgeOne,
    Edge& edgeTwo,
    Edge& edgeThree,
    Edge& edgeFour,
    int   id );

  /**
   * @brief Destructor, empty.
   */
  virtual ~Quad() {}

  /**
   * @brief Returns dimensionalty of space the quad is embedded in.
   */
  int getDimensions() const;

  /**
   * @brief Returns quad vertex with index 0, 1, 2, or 3.
   *
   * Vertex 0 is the first vertex of edge 0. Vertex 1 the second vertex of
   * edge 0. Vertex 2 is either the first or second vertex of edge 1, which
   * is determined on construction of the triangle.
   */
  Vertex& vertex( int i );

  /**
   * @brief Returns const quad vertex with index 0, 1, 2, or 3.
   *
   * Vertex 0 is the first vertex of edge 0. Vertex 1 the second vertex of
   * edge 0. Vertex 2 is either the first or second vertex of edge 1, which
   * is determined on construction of the triangle.
   */
  const Vertex& vertex( int i ) const;

  /**
   * @brief Returns quad edge with index 0, 1, 2, or 3.
   */
  Edge& edge( int i );

  /**
   * @brief Returns const quad edge with index 0, 1, 2, or 3.
   */
  const Edge& edge( int i ) const;

  /**
   * @brief Sets the outer normal of the quad.
   */
  template<typename VECTOR_T>
  void setNormal( const VECTOR_T& normal );

  /**
   * @brief Sets the center of the quad.
   */
  template<typename VECTOR_T>
  void setCenter( const VECTOR_T& center );

  /**
   * @brief Sets the radius of the circle enclosing the quad.
   */
  void setEnclosingRadius( double radius );

  /**
   * @brief Returns a among quads globally unique ID.
   */
  int getID() const;

  /**
   * @brief Returns the outer normal of the quad.
   *
   * Prerequesits: The normal has to be computed and set from outside before.
   */
  const utils::DynVector& getNormal() const;

  /**
   * @brief Returns the barycenter of the quad.
   *
   * Prerequesits: The center has to be computed and set from outside before.
   */
  const utils::DynVector& getCenter() const;

  /**
   * @brief Returns the radius of the circle enclosing the quad.
   *
   * Prerequesits: The radius has to be computed and set from outside before.
   */
  double getEnclosingRadius() const;

private:

  // @brief Edges defining the quad.
  boost::array<Edge*,4> _edges;

  // @brief Decider for choosing unique vertices from _edges.
  boost::array<int,4> _vertexMap;

  // @brief ID of the edge.
  int _id;

  // @brief Normal vector of the quad.
  utils::DynVector _normal;

  // @brief Center point of the quad.
  utils::DynVector _center;

  // @brief Minimal radius of circle enclosing the quad.
  double _enclosingRadius;
};

// --------------------------------------------------------- HEADER DEFINITIONS

inline Vertex& Quad:: vertex
(
  int i )
{
  assertion((i >= 0) && (i < 4), i);
  return edge(i).vertex(_vertexMap[i]);
};

inline const Vertex& Quad:: vertex
(
  int i ) const
{
  assertion((i >= 0) && (i < 4), i);
  return edge(i).vertex(_vertexMap[i]);
};

inline Edge& Quad:: edge
(
  int i )
{
  return *_edges[i];
};

inline const Edge& Quad:: edge
(
  int i ) const
{
  return *_edges[i];
};

template<typename VECTOR_T>
void Quad:: setNormal
(
  const VECTOR_T& normal )
{
  assertion(normal.size() == getDimensions(), normal.size(), getDimensions());
  _normal = normal;
}

template<typename VECTOR_T>
void Quad:: setCenter
(
  const VECTOR_T& center )
{
  assertion(center.size() == getDimensions(), center.size(), getDimensions());
  _center = center;
}

inline int Quad:: getID() const
{
  return _id;
}

}} // namespace precice, mesh

#endif /* PRECICE_MESH_QUAD_HPP_ */

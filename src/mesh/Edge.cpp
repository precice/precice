#include "Edge.hpp"

namespace precice {
namespace mesh {

Edge:: Edge
(
  Vertex& vertexOne,
  Vertex& vertexTwo,
  int     id )
:
  PropertyContainer (),
  _vertices( {&vertexOne, &vertexTwo} ),
  _id ( id ),
  _normal ( Eigen::VectorXd::Constant(vertexOne.getDimensions(), 0.0) ),
  _center ( Eigen::VectorXd::Constant(vertexOne.getDimensions(), 0.0) ),
  _enclosingRadius ( 0.0 )
{
  assertion ( vertexOne.getDimensions() == vertexTwo.getDimensions(),
              vertexOne.getDimensions(), vertexTwo.getDimensions() );
}

void Edge:: setEnclosingRadius
(
  double radius )
{
  _enclosingRadius = radius;
}

int Edge:: getID () const
{
  return _id;
}

const Eigen::VectorXd& Edge::getCenter () const
{
  return _center;
}

double Edge:: getEnclosingRadius () const
{
  return _enclosingRadius;
}

}} // namespace precice, mesh

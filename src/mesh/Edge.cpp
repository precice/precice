#include "Edge.hpp"
#include "utils/ManageUniqueIDs.hpp"
#include "boost/assign.hpp"

namespace precice {
namespace mesh {

Edge:: Edge
(
  Vertex& vertexOne,
  Vertex& vertexTwo,
  int     id )
:
  PropertyContainer (),
  _vertices ( boost::assign::list_of(&vertexOne)(&vertexTwo).to_array(_vertices) ),
  _id ( id ),
  _normal ( vertexOne.getDimensions(), 0.0 ),
  _center ( vertexOne.getDimensions(), 0.0 ),
  _enclosingRadius ( 0.0 )
{
  assertion ( vertexOne.getDimensions() == vertexTwo.getDimensions(),
               vertexOne.getDimensions(), vertexTwo.getDimensions() );
}

//int Edge:: getDimensions() const
//{
//  return _vertices[0]->getDimensions();
//}

//Vertex& Edge:: vertex
//(
//  int i )
//{
//  assertion ( (i == 0) || (i == 1), i );
//  return *_vertices[i];
//}
//
//const Vertex& Edge:: vertex
//(
//  int i ) const
//{
//  assertion ( (i==0) || (i==1), i );
//  return *_vertices[i];
//}

//void Edge:: setNormal
//(
//  const utils::DynVector& normal )
//{
//  assertion ( normal.size() == _vertices[0]->getDimensions(), normal,
//               _vertices[0]->getDimensions() );
//  _normal = normal;
//}

//void Edge:: setCenter
//(
//  const utils::DynVector& center )
//{
//  assertion ( center.size() == _vertices[0]->getDimensions(), center,
//               _vertices[0]->getDimensions() );
//  _center = center;
//}

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

const utils::DynVector& Edge:: getCenter () const
{
  return _center;
}

double Edge:: getEnclosingRadius () const
{
  return _enclosingRadius;
}

}} // namespace precice, mesh

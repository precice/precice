#include "Merge.hpp"
#include "mesh/PropertyContainer.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Vertex.hpp"

namespace precice {
namespace mesh {

Merge:: Merge()
:
  //PROPERTY_ID ( PropertyContainer::getFreePropertyID() ),
  _merged()
{}

Merge:: ~Merge()
{
//  for ( Vertex& vertex : _merged.vertices() ){
//    vertex.deleteProperty ( PROPERTY_ID );
//  }
//  for ( Edge& edge : _merged.edges() ){
//    edge.deleteProperty ( PROPERTY_ID );
//  }
//  for ( Triangle& triangle : _merged.triangles() ){
//    triangle.deleteProperty ( PROPERTY_ID );
//  }
}

Group& Merge:: content()
{
  return _merged;
}

const Group& Merge:: content() const
{
  return _merged;
}

//int Merge:: getPropertyIndex()
//{
//  return PROPERTY_ID;
//}
//
//bool Merge:: isUnique
//(
//  PropertyContainer& cont )
//{
//  if ( cont.hasProperty(PROPERTY_ID) ){
//    return false;
//  }
//  cont.setProperty ( PROPERTY_ID, 1 );
//  return true;
//}

}} // namespace precice, mesh


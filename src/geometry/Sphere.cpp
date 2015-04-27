// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "Sphere.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include <list>
#include <tuple>
#include "tarch/la/Scalar.h"

namespace precice {
namespace geometry {

tarch::logging::Log Sphere:: _log ( "precice::geometry::Sphere" );

Sphere:: Sphere
(
  const utils::DynVector& offset,
  double                  discretizationWidth,
  double                  radius )
:
  Geometry ( offset ),
  _discretizationWidth ( discretizationWidth ),
  _radius ( radius )
{}

void Sphere:: specializedCreate
(
  mesh::Mesh& seed )
{
  preciceTrace ( "specializedCreate()" );
  using namespace mesh;
  int dimensions = seed.getDimensions();
  assertion1 ( (dimensions == 2) || (dimensions == 3), dimensions );
  if ( dimensions == 2 ){
    utils::Vector2D start(0.0);
    start(0) += _radius;
    Vertex * oldUpper = &seed.createVertex ( -1.0 * start );
    Vertex * oldLower = oldUpper;

    int numLatitudinal = getNumberLatitudinalElements (_discretizationWidth)+1;
    double angleInc = tarch::la::PI / (double) (numLatitudinal-1.0);

    double latitude, currentRadius;
    for (int i=1; i < numLatitudinal-1; i++) {
      latitude = _radius * ( 1.0 - std::cos(i * angleInc) );
      currentRadius = _radius * std::sin (i * angleInc);
      Vertex* newUpper = & seed.createVertex(
          utils::Vector2D(-_radius + latitude, currentRadius));
      Vertex* newLower = & seed.createVertex(
          utils::Vector2D(-_radius + latitude, -1.0 * currentRadius));
      seed.createEdge ( *newUpper, *oldUpper );
      seed.createEdge ( *oldLower, *newLower );
      oldUpper = newUpper;
      oldLower = newLower;
    }
    Vertex& last = seed.createVertex (start);
    seed.createEdge (last, *oldUpper);
    seed.createEdge (*oldLower, last);
  }
  else { // dimensions == 3
    typedef std::tuple<Vertex*, Vertex*, Vertex*> TriangleVertices;
    std::list<TriangleVertices> triangles;

    // Construct initial icosahedron:
    // ------------------------------

    // An icosahderon consists of 12 vertices, 30 edges, and 20 identical
    // equilateral triangles (see http://en.wikipedia.org/wiki/Icosahedron).
    // The icosahedron is constructed from the corner vertices of 3 rectangles.
    double l = 4.0 * _radius / std::sqrt(10.0 + 2.0 * std::sqrt(5));
    //precicePrint ( "l = " << l );
    double goldenRatio = ( 1.0 + std::sqrt(5) ) / 2.0;
    double h = l * goldenRatio;

    // First rectangle lies in xy-plane:
    // ^ y     2----------3 -
    // |       |          | l
    // --> x   0----------1 -
    //         |---- h ---|
    h /= 2.0;
    l /= 2.0;
    Vertex * xy0 = & seed.createVertex ( utils::Vector3D(-h, -l, 0.0) );
    Vertex * xy1 = & seed.createVertex ( utils::Vector3D(h, -l, 0.0) );
    Vertex * xy2 = & seed.createVertex ( utils::Vector3D(-h, l, 0.0) );
    Vertex * xy3 = & seed.createVertex ( utils::Vector3D(h, l, 0.0) );

    // Second rectangle lies in zy-plane:
    // ^ z     2----------3 -
    // |       |          | l
    // --> y   0----------1 -
    //         |---- h ---|
    Vertex * yz0 = & seed.createVertex ( utils::Vector3D(0.0, -h, -l) );
    Vertex * yz1 = & seed.createVertex ( utils::Vector3D(0.0, h, -l) );
    Vertex * yz2 = & seed.createVertex ( utils::Vector3D(0.0, -h, l) );
    Vertex * yz3 = & seed.createVertex ( utils::Vector3D(0.0, h, l) );

    // Second rectangle lies in zx-plane:
    // ^ x     2----------3 -
    // |       |          | l
    // --> z   0----------1 -
    //         |---- h ---|
    Vertex * zx0 = & seed.createVertex ( utils::Vector3D(-l, 0.0, -h) );
    Vertex * zx1 = & seed.createVertex ( utils::Vector3D(-l, 0.0, h) );
    Vertex * zx2 = & seed.createVertex ( utils::Vector3D(l, 0.0, -h) );
    Vertex * zx3 = & seed.createVertex ( utils::Vector3D(l, 0.0, h) );

    // Add triangles of edge xy0-2
    triangles.push_back ( std::make_tuple(xy0, xy2, zx0) );
    triangles.push_back ( std::make_tuple(xy2, xy0, zx1) );

    // Add triangles of edge xy1-3
    triangles.push_back ( std::make_tuple(xy3, xy1, zx2) );
    triangles.push_back ( std::make_tuple(xy1, xy3, zx3) );

    // Add triangles of edge yz0-2
    triangles.push_back ( std::make_tuple(yz0, yz2, xy0) );
    triangles.push_back ( std::make_tuple(yz2, yz0, xy1) );

    // Add triangles of edge yz1-3
    triangles.push_back ( std::make_tuple(yz3, yz1, xy2) );
    triangles.push_back ( std::make_tuple(yz1, yz3, xy3) );

    // Add triangles of edge zx0-2
    triangles.push_back ( std::make_tuple(zx0, zx2, yz0) );
    triangles.push_back ( std::make_tuple(zx2, zx0, yz1) );

    // Add triangles of edge zx1-3
    triangles.push_back ( std::make_tuple(zx3, zx1, yz2) );
    triangles.push_back ( std::make_tuple(zx1, zx3, yz3) );

    // Create 3-rectangle triangles
    triangles.push_back ( std::make_tuple(xy0, zx0, yz0) );
    triangles.push_back ( std::make_tuple(xy0, yz2, zx1) );

    triangles.push_back ( std::make_tuple(xy2, yz1, zx0) );
    triangles.push_back ( std::make_tuple(xy2, zx1, yz3) );

    triangles.push_back ( std::make_tuple(xy1, yz0, zx2) );
    triangles.push_back ( std::make_tuple(xy1, zx3, yz2) );

    triangles.push_back ( std::make_tuple(xy3, zx2, yz1) );
    triangles.push_back ( std::make_tuple(xy3, yz3, zx3) );


    // Refine initial icosahedron to obey discretization width:
    // --------------------------------------------------------

    double sidelength = l * 2.0;
    while ( sidelength > _discretizationWidth ) {
      std::map<std::pair<int,int>,Vertex*> dividedEdges;
      std::list<TriangleVertices>::iterator iter = triangles.begin();
      sidelength = 0.0;
      while ( iter != triangles.end() ) {
        Vertex * v0 = std::get<0> ( *iter );
        Vertex * v1 = std::get<1> ( *iter );
        Vertex * v2 = std::get<2> ( *iter );

        Vertex * v01 = getVertex ( *v0, *v1, dividedEdges, seed );
        Vertex * v12 = getVertex ( *v1, *v2, dividedEdges, seed );
        Vertex * v20 = getVertex ( *v2, *v0, dividedEdges, seed );
        assertion ( v01 != NULL );
        assertion ( v12 != NULL );
        assertion ( v20 != NULL );
        triangles.insert ( iter, std::make_tuple(v0, v01, v20) );
        triangles.insert ( iter, std::make_tuple(v01, v1, v12 ) );
        triangles.insert ( iter, std::make_tuple(v20, v12, v2 ) );
        triangles.insert ( iter, std::make_tuple(v20, v01, v12) );
        iter = triangles.erase ( iter );

        /*
        double length01 = la::norm ( v1->getCoords() - v0->getCoords() );
        double length12 = la::norm ( v1->getCoords() - v2->getCoords() );
        double length20 = la::norm ( v2->getCoords() - v0->getCoords() );
        bool refine = false;
        if ( length01 > _discretizationWidth ) {
          refine = true;
          sidelength = length01 * 0.5;
        }
        else if ( length12 > _discretizationWidth ) {
          refine = true;
          sidelength = length12 * 0.5;
        }
        else if ( length20 > _discretizationWidth ) {
          refine = true;
          sidelength = length20 * 0.5;
        }
        if ( refine ) {
          Vertex * v01 = getVertex ( *v0, *v1, dividedEdges, seed );
          Vertex * v12 = getVertex ( *v1, *v2, dividedEdges, seed );
          Vertex * v20 = getVertex ( *v2, *v0, dividedEdges, seed );
          assertion ( v01 != NULL );
          assertion ( v12 != NULL );
          assertion ( v20 != NULL );
          triangles.insert ( iter, std::make_tuple(v0, v01, v20) );
          triangles.insert ( iter, std::make_tuple(v01, v1, v12 ) );
          triangles.insert ( iter, std::make_tuple(v20, v12, v2 ) );
          triangles.insert ( iter, std::make_tuple(v20, v01, v12) );
          iter = triangles.erase ( iter );
        }
        */

        /*
        Vertex * v01 = getVertex ( *v0, *v1, dividedEdges, seed );
        Vertex * v12 = getVertex ( *v1, *v2, dividedEdges, seed );
        Vertex * v20 = getVertex ( *v2, *v0, dividedEdges, seed );
        Vector coords ( v0->getCoords() );
        coords += v1->getCoords();
        coords += v2->getCoords();
        coords /= 3.0;
        coords *= _radius / la::norm(coords);
        Vertex * vMid = & seed.createVertex ( coords );
        assertion ( v01 != NULL );
        assertion ( v12 != NULL );
        assertion ( v20 != NULL );
        triangles.insert ( iter, std::make_tuple(v0, v01, v20) );
        triangles.insert ( iter, std::make_tuple(v20, v01, vMid) );
        triangles.insert ( iter, std::make_tuple(vMid, v01, v12) );
        triangles.insert ( iter, std::make_tuple(v01, v1, v12 ) );
        triangles.insert ( iter, std::make_tuple(v20, vMid, v12 ) );
        triangles.insert ( iter, std::make_tuple(v20, v12, v2 ) );
        iter = triangles.erase ( iter );
        */

        /*
        Vertex * v01 = getVertex ( *v0, *v1, dividedEdges, seed );
        Vertex * v12 = getVertex ( *v1, *v2, dividedEdges, seed );
        Vertex * v20 = getVertex ( *v2, *v0, dividedEdges, seed );
        Vector coords ( v0->getCoords() );
        coords += v1->getCoords();
        coords += v2->getCoords();
        coords /= 3.0;
        coords *= _radius / la::norm(coords);
        Vertex * vMid = & seed.createVertex ( coords );
        assertion ( v01 != NULL );
        assertion ( v12 != NULL );
        assertion ( v20 != NULL );
        triangles.insert ( iter, std::make_tuple(v0, vMid, v20) );
        triangles.insert ( iter, std::make_tuple(v0, v01, vMid) );
        triangles.insert ( iter, std::make_tuple(vMid, v01, v1) );
        triangles.insert ( iter, std::make_tuple(vMid, v1, v12 ) );
        triangles.insert ( iter, std::make_tuple(v2, vMid, v12 ) );
        triangles.insert ( iter, std::make_tuple(v20, vMid, v2 ) );
        iter = triangles.erase ( iter );
        */

        /*
        Vertex * v01 = getVertex ( *v0, *v1, dividedEdges, seed );
        Vertex * v12 = getVertex ( *v1, *v2, dividedEdges, seed );
        Vertex * v20 = getVertex ( *v2, *v0, dividedEdges, seed );
        Vector coords ( v0->getCoords() );
        coords += v1->getCoords();
        coords += v2->getCoords();
        coords /= 3.0;
        coords *= _radius / la::norm(coords);
        Vertex * vMid = & seed.createVertex ( coords );
        assertion ( v01 != NULL );
        assertion ( v12 != NULL );
        assertion ( v20 != NULL );
        triangles.insert ( iter, std::make_tuple(v0, vMid, v20) );
        triangles.insert ( iter, std::make_tuple(v0, v01, vMid) );
        triangles.insert ( iter, std::make_tuple(vMid, v01, v1) );
        triangles.insert ( iter, std::make_tuple(vMid, v1, v12 ) );
        triangles.insert ( iter, std::make_tuple(v2, vMid, v12 ) );
        triangles.insert ( iter, std::make_tuple(v20, vMid, v2 ) );
        iter = triangles.erase ( iter );
        */

        /*
        Vector coords ( v0->getCoords() );
        coords += v1->getCoords();
        coords += v2->getCoords();
        coords /= 3.0;
        coords *= _radius / la::norm(coords);
        Vertex * vMid = & seed.createVertex ( coords );
        triangles.insert ( iter, std::make_tuple(v0, v1, vMid) );
        triangles.insert ( iter, std::make_tuple(v2, vMid, v1) );
        triangles.insert ( iter, std::make_tuple(v0, vMid, v2) );
        iter = triangles.erase ( iter );
        */

         //      assertion ( triangles.size() < 1000 ); // TODO remove
        double length = tarch::la::norm2 ( v20->getCoords() - v01->getCoords() );
        sidelength = length > sidelength ? length : sidelength;
      }
    }

    // Create edges and triangles of final icosahedron:
    // ------------------------------------------------

    std::map<std::pair<int,int>, Edge*> edges;
    for (TriangleVertices & vertices : triangles) {
      Vertex * v0 = std::get<0> ( vertices );
      Vertex * v1 = std::get<1> ( vertices );
      Vertex * v2 = std::get<2> ( vertices );
      assertion ( v0 != NULL );
      assertion ( v1 != NULL );
      assertion ( v2 != NULL );
      Edge * e0 = getEdge ( *v0, *v1, edges, seed );
      Edge * e1 = getEdge ( *v1, *v2, edges, seed );
      Edge * e2 = getEdge ( *v2, *v0, edges, seed );
      assertion ( e0 != NULL );
      seed.createTriangle ( *e0, *e1, *e2 );
    }
  }
}

mesh::Vertex* Sphere:: getVertex
(
  mesh::Vertex&                               v0,
  mesh::Vertex&                               v1,
  std::map<std::pair<int,int>,mesh::Vertex*>& dividedEdges,
  mesh::Mesh&                                 seed )
{
  std::map<std::pair<int,int>,mesh::Vertex*>::iterator iter;
  std::pair<int,int> index = std::make_pair ( v0.getID(), v1.getID() );
  iter = dividedEdges.find ( index );
  if ( iter != dividedEdges.end() ) {
    return iter->second;
  }
  index = std::make_pair ( v1.getID(), v0.getID() );
  iter = dividedEdges.find ( index );
  if ( iter != dividedEdges.end() ) {
    return iter->second;
  }
  mesh::Vertex* vertex = NULL;
  if ( seed.getDimensions() == 2 ){ // 2D
    utils::Vector2D coords ( v0.getCoords() );
    coords += v1.getCoords();
    coords /= 2.0;
    coords *= _radius / tarch::la::norm2(coords);
    vertex = &seed.createVertex ( coords );
  }
  else { // 3D
    utils::Vector3D coords ( v0.getCoords() );
    coords += v1.getCoords();
    coords /= 2.0;
    coords *= _radius / tarch::la::norm2(coords);
    vertex = &seed.createVertex ( coords );
  }
  dividedEdges[index] = vertex;
  return vertex;
}

mesh::Edge* Sphere:: getEdge
(
  mesh::Vertex&                             v0,
  mesh::Vertex&                             v1,
  std::map<std::pair<int,int>,mesh::Edge*>& edges,
  mesh::Mesh&                               seed )
{
  std::map<std::pair<int,int>,mesh::Edge*>::iterator iter;
  std::pair<int,int> index = std::make_pair ( v0.getID(),v1.getID() );
  iter = edges.find ( index );
  if ( iter != edges.end() ) {
    return iter->second;
  }
  index = std::make_pair ( v1.getID(),v0.getID() );
  iter = edges.find ( index );
  if ( iter != edges.end() ) {
    return iter->second;
  }
  mesh::Edge* edge = & seed.createEdge ( v0, v1 );
  edges[index] = edge;
  return edge;
}

int Sphere:: getNumberLongitudinalElements
(
  double discretizationWidth ) const
{
  return getOffset().size() == 2 // dimensions == 2
         ? 2
         : (2.0 * _radius * tarch::la::PI / discretizationWidth) + 0.5;
}

int Sphere:: getNumberLatitudinalElements
(
  const double discretizationWidth ) const
{
  return static_cast<int>((1.6 * _radius * tarch::la::PI / discretizationWidth) + 0.5);
}

}} // namespace precice, geometry

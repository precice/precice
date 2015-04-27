// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "FindClosestTest.hpp"
#include "query/FindClosest.hpp"
#include "query/ExportVTKNeighbors.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/PropertyContainer.hpp"
#include "io/ExportVTK.hpp"
#include "utils/Parallel.hpp"
#include <vector>
#include <memory>

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::query::tests::FindClosestTest)

namespace precice {
namespace query {
namespace tests {

tarch::logging::Log FindClosestTest:: _log ( "precice::query::tests::FindClosestTest" );

FindClosestTest:: FindClosestTest ()
:
  TestCase ("query::FindClosestTest")
{}

void FindClosestTest:: run ()
{
  PRECICE_MASTER_ONLY {
    testMethod ( testFindClosestDistanceToVertices );
    testMethod ( testSignOfShortestDistance );
    testMethod ( testIndependenceOfSignOfShortestDistance );
    testMethod ( testFindClosestDistanceToEdges );
    testMethod ( testFindClosestDistanceToEdges3D );
    testMethod ( testFindClosestDistanceToTriangles );
    testMethod ( testFindClosestDistanceToTrianglesAndVertices );
    testMethod ( testMultipleGeometryIDs );
    testMethod ( testWeigthsOfVertices );
  }
}

void FindClosestTest:: testFindClosestDistanceToVertices ()
{
   preciceTrace ( "testFindClosestDistanceToVertices()" );
   for ( int dim=2; dim <= 3; dim++ ){
     preciceDebug ( "Dimension = " << dim );
     mesh::Mesh mesh ( "RootMesh", dim, false );
     mesh.createVertex ( utils::DynVector(dim, 0.0) );
     utils::DynVector queryCoords0 ( dim, 0.0 );
     queryCoords0[0] = 1.0;
     FindClosest find ( queryCoords0 );
     find ( mesh );
     query::ClosestElement closest = find.getClosest();
     double distance = closest.distance;
     validateNumericalEquals (distance, 1.0);
     if ( dim == 3 ){
       mesh::Mesh mesh3D ( "Mesh3D", dim, false );
       mesh3D.createVertex ( utils::DynVector(3, 1.0) );
       FindClosest find2 ( utils::DynVector(3, -1.0) );
       find2 ( mesh3D );
       distance = find2.getClosest().distance;
       validateNumericalEquals ( distance, tarch::la::norm2(utils::DynVector(3,2.0)) );
     }
   }
}

void FindClosestTest:: testSignOfShortestDistance ()
{
  preciceTrace ( "testSignOfShortestDistance()" );
  for ( int dim=2; dim <= 3; dim++ ){
    preciceDebug ( "Dimension = " << dim );
    mesh::Mesh mesh ( "Mesh", dim, false );
    mesh::Vertex & vertex = mesh.createVertex ( utils::DynVector(dim, 0.0) );
    utils::DynVector normal ( dim, 0.0 );
    normal[0] = 1.0;
    vertex.setNormal ( normal );

    // Check point that lies outside of geometry
    utils::DynVector queryCoords ( dim, 0.0 );
    queryCoords(0) = 1.0;
    FindClosest find ( queryCoords );
    find ( mesh );
    double distance = find.getClosest().distance;
    validateEquals (distance, 1.0);

    // Check point that lies inside of geometry
    normal[0] = -1.0;
    vertex.setNormal ( normal );
    find.reset ();
    find ( mesh );
    validateEquals (find.getClosest().distance, -1.0);
  }
}

void FindClosestTest:: testIndependenceOfSignOfShortestDistance ()
{
  preciceTrace ( "testIndependenceOfSignOfShortestDistance()" );
  using utils::DynVector;
  for ( int dim=2; dim <= 3; dim++ ){
    preciceDebug ( "Dimension = " << dim );
    mesh::Mesh mesh ( "Mesh", dim, false );
    mesh::Vertex& vertex = mesh.createVertex ( DynVector(dim, 1.0) );
    vertex.setNormal ( DynVector(dim, 1.0) );
    mesh::Vertex& vertex2 = mesh.createVertex ( DynVector(dim, -2.0) );
    vertex2.setNormal ( DynVector(dim, 1.0) );

    FindClosest find ( DynVector(dim, 0.0) );
    find ( mesh );
    double distance = find.getClosest().distance;
    validateEquals ( distance, -1.0 * tarch::la::norm2(DynVector(dim,1.0)));

    // Invert normal of futher away vertex, should have no influence
    DynVector normal ( dim, 0.0 );
    normal[1] = -1.0;
    vertex2.setNormal ( normal );
    find.reset ();
    find ( mesh );
    distance = find.getClosest().distance;
    validateNumericalEquals (  distance, -1.0 * tarch::la::norm2(DynVector(dim,1.0)) );

    // Invert normal of closer vertex, should invert sign of distance
    vertex.setNormal ( normal );
    find.reset ();
    find ( mesh );
    distance = find.getClosest().distance;
    validateNumericalEquals ( distance, tarch::la::norm2(DynVector(dim,1.0)) );
  }
}

void FindClosestTest:: testFindClosestDistanceToEdges ()
{
  preciceTrace ( "testFindClosestDistanceToEdges()" );
  using utils::DynVector;
  for ( int dim=2; dim <= 3; dim++ ){
    preciceDebug ( "Dimensions = " << dim );
    // Create geometry consisting of two vertices and an edge
    mesh::Mesh mesh ( "Mesh", dim, false );
    mesh::Vertex& v1 = mesh.createVertex ( DynVector(dim, -1.0) );
    mesh::Vertex& v2 = mesh.createVertex ( DynVector(dim, 1.0) );
    mesh::Edge& edge = mesh.createEdge ( v1, v2 );
    DynVector normal ( dim, -1.0 );
    normal[1] = 1.0;
    v1.setNormal ( normal );
    v2.setNormal ( normal );
    edge.setNormal ( normal );

    // Create query points
    // Query point 0 lies outside of the geometry
    DynVector queryCoords0 ( dim, -1.0 );
    queryCoords0[1] = 1.0;
    // Query point 1 lies inside of the geometry
    DynVector queryCoords1 ( dim, 1.0 );
    queryCoords1[1] = -1.0;
    // Query point 2 lies on the query edge
    DynVector queryCoords2 ( dim, 0.0 );
    // Query point 3 lies on the same line as the edge, but above
    DynVector queryCoords3 ( dim, 2.0 );
    // Query point 4 lies on the same line as the edge, but below
    DynVector queryCoords4 ( dim, -2.0 );

    // Create query objects
    FindClosest find0 ( queryCoords0 );
    FindClosest find1 ( queryCoords1 );
    FindClosest find2 ( queryCoords2 );
    FindClosest find3 ( queryCoords3 );
    FindClosest find4 ( queryCoords4 );

    // Perform queries
    find0 ( mesh );
    find1 ( mesh );
    find2 ( mesh );
    find3 ( mesh );
    find4 ( mesh );

    // Evaluate query results
    using tarch::la::norm2;
    DynVector expected ( dim, 0.0 );
    assign(expected) = 1.0;
    validateNumericalEquals (find0.getClosest().distance, norm2(expected));
    validateNumericalEquals (find1.getClosest().distance, -1.0 * norm2(expected));
    validateNumericalEquals (find2.getClosest().distance, 0.0);
    validateNumericalEquals (std::abs(find3.getClosest().distance), norm2(expected));
    validateNumericalEquals (std::abs(find4.getClosest().distance), norm2(expected));
  }
}

void FindClosestTest:: testFindClosestDistanceToEdges3D ()
{
   preciceTrace ( "testFindClosestDistanceToEdges3D()" );
   int dim = 3;
   using utils::Vector3D;
   // Create geometry consisting of two vertices and an edge
   mesh::Mesh mesh ( "Mesh", dim, false );
   mesh::Vertex& v1 = mesh.createVertex ( Vector3D(-1.0, -1.0, 0.0) );
   mesh::Vertex& v2 = mesh.createVertex ( Vector3D( 1.0,  1.0, 0.0) );
   mesh::Edge& edge = mesh.createEdge ( v1, v2 );
   utils::DynVector normal(dim);
   assignList(normal) = -1.0, 1.0, 0.0;
   v1.setNormal ( normal );
   v2.setNormal ( normal );
   edge.setNormal ( normal );

   io::ExportVTK exportMesh(true);
   exportMesh.doExport ( "query-FindClosestTest-testFindClosestDistanceToEdges3D",
                         mesh );

   // Create query points
   std::vector<Vector3D> queryPoints;
   // Query point 0 lies outside of the geometry
   queryPoints.push_back ( Vector3D(-0.5,  0.0,  0.0) );
   queryPoints.push_back ( Vector3D( 0.0,  0.5,  0.0) );
   queryPoints.push_back ( Vector3D(-0.5,  0.5,  0.0) );
   queryPoints.push_back ( Vector3D(-0.5,  0.5, -0.5) );
   queryPoints.push_back ( Vector3D(-0.5,  0.5,  0.5) );

   // Create query objects
   std::vector<query::FindClosest*> finds;
   for ( size_t i=0; i < queryPoints.size(); i++ ) {
      finds.push_back ( new FindClosest(queryPoints[i]) );
   }

   // Perform queries
   ExportVTKNeighbors exportNeighbors;
   for ( size_t i=0; i < finds.size(); i++ ) {
      validate ( (*finds[i])(mesh) );
      exportNeighbors.addNeighbors ( queryPoints[i], finds[i]->getClosest() );
   }

   exportNeighbors.exportNeighbors (
      "query-FindClosestTest-testFindClosestDistanceToEdges3D-neighbors" );

   // Evaluate query results
   validateNumericalEquals ( finds[0]->getClosest().distance, std::sqrt(1.0/8.0) );
   validateNumericalEquals ( finds[1]->getClosest().distance, std::sqrt(1.0/8.0) );
   validateNumericalEquals ( finds[2]->getClosest().distance,
                             tarch::la::norm2(Vector3D(0.5, -0.5, 0.0)) );
   validateNumericalEquals ( finds[3]->getClosest().distance,
                             tarch::la::norm2(Vector3D(0.5, -0.5, 0.5)) );
   validateNumericalEquals ( finds[4]->getClosest().distance,
                             tarch::la::norm2(Vector3D(0.5, -0.5, -0.5)) );

   // Clean up
   for ( size_t i=0; i < finds.size(); i++ ) {
      delete finds[i];
   }
}

void FindClosestTest:: testFindClosestDistanceToTriangles ()
{
  preciceTrace ( "testFindClosestDistanceToTriangles()" );
  using namespace tarch::la;
  using utils::Vector3D;

  // Create mesh to query
  mesh::Mesh mesh ( "Mesh", 3, true );
  mesh::Vertex& v0 = mesh.createVertex ( Vector3D(0.0, 0.0, 0.0) );
  mesh::Vertex& v1 = mesh.createVertex ( Vector3D(1.0, 1.0, 0.0) );
  mesh::Vertex& v2 = mesh.createVertex ( Vector3D(1.0, 1.0, 1.0) );
  mesh::Edge& e0 = mesh.createEdge ( v0, v1 );
  mesh::Edge& e1 = mesh.createEdge ( v1, v2 );
  mesh::Edge& e2 = mesh.createEdge ( v2, v0 );
  mesh.createTriangle ( e0, e1, e2 );
  mesh.computeState();

  // Prepare and issue queries
  std::vector<Vector3D> queries;
  queries.push_back ( Vector3D( 0.6,  0.6,  0.5) ); //  0: on triangle middle
  queries.push_back ( Vector3D( 0.0,  0.0,  0.0) ); //  1: on vertex0
  queries.push_back ( Vector3D( 1.0,  1.0,  0.0) ); //  2: on vertex1
  queries.push_back ( Vector3D( 1.0,  1.0,  1.0) ); //  3: on vertex2
  queries.push_back ( Vector3D( 0.5,  0.5,  0.0) ); //  4: on edge0
  queries.push_back ( Vector3D( 1.0,  1.0,  0.5) ); //  5: on edge1
  queries.push_back ( Vector3D( 0.5,  0.5,  0.5) ); //  6: on edge2
  queries.push_back ( Vector3D( 0.0,  1.0,  0.3) ); //  7: outside triangle
  queries.push_back ( Vector3D( 1.0,  0.0,  0.3) ); //  8: inside triangle
  queries.push_back ( Vector3D(-0.5, -0.5, -0.5) ); //  9: outside vertex0
  queries.push_back ( Vector3D( 1.5,  1.5, -0.5) ); // 10: outside vertex1
  queries.push_back ( Vector3D( 1.5,  1.5,  1.5) ); // 11: outside vertex2
  std::vector< std::shared_ptr<FindClosest> > findVisitors;
  for ( size_t i=0; i < queries.size(); i++ ) {
    std::shared_ptr<FindClosest> find ( new FindClosest(queries[i]) );
    findVisitors.push_back ( find );
    (*findVisitors[i])( mesh );
  }

  // Validate results
  for ( size_t i=0; i < 7; i++ ) {
    validateNumericalEquals ( findVisitors[i]->getClosest().distance, 0.0 );
  }
  Vector3D expect;
  assignList(expect) = 0.5, -0.5, 0.0;
  validate (equals(findVisitors[7]->getClosest().vectorToElement, expect));
  validateNumericalEquals (findVisitors[7]->getClosest().distance, norm2(expect));
  assignList(expect) = -0.5, 0.5, 0.0;
  validate ( equals(findVisitors[8]->getClosest().vectorToElement, expect) );
  validateNumericalEquals ( findVisitors[8]->getClosest().distance,
                            -1.0 * norm2(expect) );
  assignList(expect) = 0.5, 0.5, 0.5;
  validate ( equals(findVisitors[9]->getClosest().vectorToElement, expect) );
  validateNumericalEquals ( std::abs(findVisitors[9]->getClosest().distance),
                            norm2(expect) );
  assignList(expect) = -0.5, -0.5, 0.5;
  validate ( equals(findVisitors[10]->getClosest().vectorToElement, expect) );
  validateNumericalEquals ( std::abs(findVisitors[10]->getClosest().distance),
                            norm2(expect) );
  assignList(expect) = -0.5, -0.5, -0.5;
  validate ( equals(findVisitors[11]->getClosest().vectorToElement, expect) );
  validateNumericalEquals ( std::abs(findVisitors[11]->getClosest().distance),
                            norm2(expect) );
}


void FindClosestTest:: testFindClosestDistanceToTrianglesAndVertices ()
{
  preciceTrace ( "testFindClosestDistanceToTrianglesAndVertices()" );
  int dim = 2;
  using utils::Vector2D;
  mesh::Mesh mesh ( "Mesh", dim, false );
  mesh::Vertex& vertex1 = mesh.createVertex (Vector2D(0.0, 0.0));
  vertex1.setNormal (Vector2D(-0.5, 0.5));

  mesh::Vertex& vertex2 = mesh.createVertex (Vector2D(1.0, 0.0));
  vertex2.setNormal (Vector2D(0.5, 0.5));

  mesh::Edge& edge = mesh.createEdge (vertex1, vertex2);
  edge.setNormal (Vector2D(0.0, 1.0));

  query::FindClosest find (Vector2D(0.0, 0.0));
  find ( mesh );
  double distance = find.getClosest().distance;
  validateNumericalEquals (distance, 0.0);

  query::FindClosest find2 (Vector2D(0.5, 0.0));
  find2 ( mesh );
  distance = find2.getClosest().distance;
  validateNumericalEquals (distance, 0.0);

  query::FindClosest find3 (Vector2D(0.5, 0.1));
  find3 ( mesh );
  distance = find3.getClosest().distance;
  validateNumericalEquals (distance, 0.1);

  query::FindClosest find4 (Vector2D(0.0, 1.5));
  find4 ( mesh );
  distance = find4.getClosest().distance;
  validateNumericalEquals (distance, 1.5);

  query::FindClosest find5 (Vector2D(0.5, -1.0));
  find5 ( mesh );
  distance = find5.getClosest().distance;
  validateNumericalEquals (distance, -1.0);
//#   elif defined(Dim3)
//    mesh::Vertex & vertex1 = mesh.createVertex(Vector(0.0, 0.0, 0.0));
//    vertex1.setNormal (Vector(-0.5, -0.5, 0.5));
//    (edge.INDEX_GEOMETRY_ID, geometryID)
//
//    mesh::Vertex & vertex2 = mesh.createVertex(Vector(1.0, 0.0, 0.0));
//    vertex2.setNormal (Vector(0.5, -0.5, 0.5));
//
//    mesh::Vertex & vertex3 = mesh.createVertex(Vector(0.0, 1.0, 0.0));
//    vertex3.setNormal (Vector(-0.5, 0.5, 0.5));
//
//    int iVertex1 = vertex1.getIndex ();
//    int iVertex2 = vertex2.getIndex ();
//    int iVertex3 = vertex3.getIndex ();
//
//    mesh::Triangle & triangle = mesh.createTriangle (IntVector(iVertex1, iVertex2, iVertex3));
//    triangle.setNormal (Vector(0.0, 0.0, 1.0));
//
//    query::FindClosestOnMeshVisitor find (Vector(0.0, 0.0, 0.0));
//    container.accept (find);
//    double distance = find.getClosest().distance;
//    validateNumericalEquals (distance, 0.0);
//
//    query::FindClosestOnMeshVisitor find2 (Vector(0.5, 0.5, 0.0));
//    container.accept (find2);
//    distance = find2.getClosest().distance;
//    validateNumericalEquals (distance, 0.0);
//
//    query::FindClosestOnMeshVisitor find3 (Vector(0.5, 0.5, 0.1));
//    container.accept (find3);
//    distance = find3.getClosest().distance;
//    validateNumericalEquals (distance, 0.1);
//
//    query::FindClosestOnMeshVisitor find4 (Vector(0.0, 1.5, 0.0));
//    container.accept (find4);
//    distance = find4.getClosest().distance;
//    validateNumericalEquals (distance, 0.5);
//
//    query::FindClosestOnMeshVisitor find5 (Vector(1.0, 1.0, 0.0));
//    container.accept (find5);
//    distance = find5.getClosest().distance;
//    validateNumericalEquals (distance, -1.0);
}

void FindClosestTest:: testMultipleGeometryIDs ()
{
  preciceTrace ( "testMultipleGeometryIDs" );
  using utils::DynVector;
  int dim = 2;
  mesh::Mesh mesh ( "Mesh", dim, true );
  query::ExportVTKNeighbors exportNeighbors;
  std::vector<mesh::Vertex*> vertices(dim);
  for (int i=0; i < dim; i++){
    DynVector vertexCoords (dim, 0.0);
    vertexCoords[i] = 1.0;
    vertices[i] = & mesh.createVertex(vertexCoords);
  }
  mesh::Edge& face = mesh.createEdge (*vertices[0], *vertices[1]);
  face.addParent ( mesh.setSubID("face-2") );
  int idFace = mesh.getID ( "Mesh-face-2" );
  int idsVertices[2];
  mesh.computeState();
  for (int i=0; i < 2; i++){
    std::ostringstream stream;
    stream << "vertex-" << i;
    face.vertex(i).addParent ( mesh.setSubID(stream.str()) );
    idsVertices[i] = mesh.getID ( "Mesh-" + stream.str() );
  }

  // Perform queries
  std::vector<double> distances(dim);
  std::vector<std::vector<int> > geoIDs(dim);
  query::ClosestElement closest(dim);
  for (int i=0; i < dim; i++){
    DynVector query (dim, 0.0);
    query[i] = 2.0;
    FindClosest findClosest ( query );
    validate ( findClosest(mesh) );
    closest = findClosest.getClosest();
    exportNeighbors.addNeighbors ( query, closest );
    distances[i] = closest.distance;
    geoIDs[i] = closest.meshIDs;
  }
  DynVector query (dim, 0.0);
  FindClosest findClosest ( query );
  validate ( findClosest(mesh) );
  closest = findClosest.getClosest();
  exportNeighbors.addNeighbors ( query, closest );
  double faceDistance = closest.distance;
  std::vector<int> faceGeoIDs = closest.meshIDs;

  // Visualize queries
  io::ExportVTK exportVTK(true);
  exportVTK.doExport ("query-FindClosestTest-testMultipleGeometryIDs.vtk", mesh);
  exportNeighbors.exportNeighbors (
      "query-FindClosestTest-testMultipleGeometryIDs_neighb.vtk");

  // Validate queries
  for (int i=0; i < dim; i++){
    validateNumericalEquals (distances[i], -1.0);
    validateEquals (geoIDs[i][1], idsVertices[i]);
  }
  validateNumericalEquals (faceDistance, std::sqrt(1.0/(double)dim));
  validateNumericalEquals (faceGeoIDs[1], idFace);
}

void FindClosestTest:: testWeigthsOfVertices ()
{
  preciceTrace ( "testWeigthsOfVertices()" );
  int dim = 2;
  using namespace tarch::la;
  using utils::Vector2D;

  // Create geometry
  mesh::Mesh mesh ( "Mesh", dim, true );
  mesh.setProperty (mesh.INDEX_GEOMETRY_ID, 0);
  mesh::Vertex& vertex1 = mesh.createVertex (Vector2D(0.0, 0.0));
  mesh::Vertex& vertex2 = mesh.createVertex (Vector2D(1.0, 0.0));
  mesh.createEdge (vertex1, vertex2);
  mesh.computeState();

  // Query elements
  query::FindClosest findClosest (Vector2D(0.3, 1.0));
  findClosest ( mesh );
  const query::ClosestElement& closest = findClosest.getClosest ();

  // Validate results
  validateEquals ( closest.interpolationElements.size(), 2 );
  mesh::Vertex* const pointerVertex1 =
    dynamic_cast<mesh::Vertex* const> (closest.interpolationElements[0].element);
  mesh::Vertex* pointerVertex2 =
    dynamic_cast<mesh::Vertex*> (closest.interpolationElements[1].element);
  validate ( equals(pointerVertex1->getCoords(), Vector2D(0.0, 0.0)) );
  validate ( equals(pointerVertex2->getCoords(), Vector2D(1.0, 0.0)) );
  validateNumericalEquals ( closest.interpolationElements[0].weight, 0.7 );
  validateNumericalEquals ( closest.interpolationElements[1].weight, 0.3 );
}

}}} // namespace precice, query, tests

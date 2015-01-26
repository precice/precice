// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ModifyCoordinatesActionTest.hpp"
#include "action/ModifyCoordinatesAction.hpp"
#include "geometry/Cuboid.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "utils/Parallel.hpp"
#include "tarch/la/WrappedVector.h"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::action::tests::ModifyCoordinatesActionTest)

namespace precice {
namespace action {
namespace tests {

tarch::logging::Log ModifyCoordinatesActionTest::
   _log ( "precice::action::tests::ModifyCoordinatesActionTest" );

ModifyCoordinatesActionTest:: ModifyCoordinatesActionTest ()
:
  TestCase ( "precice::action::tests::ModifyCoordinatesActionTest" )
{}

void ModifyCoordinatesActionTest:: run ()
{
  PRECICE_MASTER_ONLY {
    testMethod ( testAddToCoordinates );
    testMethod ( testSubtractFromCoordinates );
  }
}

//void ModifyCoordinatesActionTest:: testSetAsDisplacement ()
//{
//   preciceDebug ( "testSetAsDisplacement()", "Entering" );
//
//   // Create geometryContext by faking a geometry but not using it to create
//   // the mesh. The mesh is created by hand, such that references to the vertices
//   // to be displaced are obtained.
//   action::MeshContext geoContext;
//   geoContext.geometry = geometry::PtrGeometry (
//      new geometry::Cuboid("test-cuboid", true, Vector(1.0), 1.0, Vector(1.0)) );
//   int dataID = 0;
//   mesh::Data data ( dataID, mesh::Data::TYPE_VECTOR );
//   geoContext.mesh = mesh::PtrMesh ( new mesh::Mesh() );
//   geoContext.mesh->setVertexData ( data );
//   geoContext.associatedData.push_back ( dataID );
//   mesh::Vertex & v0 = geoContext.mesh->createVertex ( Vector(0.0) );
//   mesh::Vertex & v1 = geoContext.mesh->createVertex ( Vector(1.0) );
//   mesh::Edge & edge = geoContext.mesh->createEdge ( v0, v1 );
//   mesh::ComputeMeshState computeState;
//   computeState.computeState ( *geoContext.mesh, true );
//
//   // Create ApplyDisplacementsMeshAction
//   action::ApplyDisplacementsMeshAction applyDisplacements (
//      data, action::ApplyDisplacementsMeshAction::SET_DISPLACEMENTS_MODE );
//   applyDisplacements.loadMeshContext ( & geoContext );
//
//   // Validate coordinates of mesh before modifying it
//   validate ( tarch::la::equals(v0.getCoords(), Vector(0.0)),
//              "testSetAsDisplacement" );
//   validate ( tarch::la::equals(v1.getCoords(), Vector(1.0)),
//              "testSetAsDisplacement" );
//   Vector normalizedNormal( Vector(0.5, -0.5) );
//   normalizedNormal = normalizedNormal /
//                      utils::LinearAlgebra::twoNorm(normalizedNormal);
//   validate ( tarch::la::equals(edge.getNormal(), normalizedNormal),
//              "testSetAsDisplacement" );
//
//   // Set displacements
//   v0.setProperty ( dataID, Vector(2.0) );
//   v1.setProperty ( dataID, Vector(-1.0) );
//
//   // Apply displacements to  node coordinates
//   applyDisplacements.performAction ();
//
//   // Validate coordinates of mesh after modifying it
//   validate ( tarch::la::equals(v0.getCoords(), Vector(3.0)),
//              "testSetAsDisplacement" );
//   validate ( tarch::la::equals(v1.getCoords(), Vector(0.0)),
//              "testSetAsDisplacement" );
//   normalizedNormal = Vector(-0.5, 0.5);
//   normalizedNormal = normalizedNormal /
//                      utils::LinearAlgebra::twoNorm(normalizedNormal);
//   validate ( tarch::la::equals(edge.getNormal(), normalizedNormal),
//              "testSetAsDisplacement" );
//
//   preciceDebug ( "testSetAsDisplacement()", "Leaving" );
//}

void ModifyCoordinatesActionTest:: testAddToCoordinates ()
{
   preciceTrace ( "testAddToCoordinates()" );
   using namespace mesh;
   using utils::Vector2D;
   // Create geometryContext by faking a geometry but not using it to create
   // the mesh. The mesh is created by hand, such that references to the vertices
   // to be displaced are obtained.
   PtrMesh mesh ( new Mesh("Mesh", 2, false) );
   PtrData data = mesh->createData ( "test-data", 2 );
   int dataID = data->getID ();
   Vertex& v0 = mesh->createVertex ( Vector2D(0.0) );
   Vertex& v1 = mesh->createVertex ( Vector2D(1.0) );
   Edge & edge = mesh->createEdge ( v0, v1 );
   mesh->computeState();
   mesh->allocateDataValues ();
   utils::DynVector& values = data->values ();

   // Create ApplyDisplacementsMeshAction
   action::ModifyCoordinatesAction modifyCoordinates (
      action::ModifyCoordinatesAction::ALWAYS_PRIOR, dataID, mesh,
      action::ModifyCoordinatesAction::ADD_TO_COORDINATES_MODE );

   // Validate coordinates of mesh before modifying it
   validate ( tarch::la::equals(v0.getCoords(), Vector2D(0.0)) );
   validate ( tarch::la::equals(v1.getCoords(), Vector2D(1.0)) );
   Vector2D normalizedNormal( Vector2D(0.5, -0.5) );
   normalizedNormal = normalizedNormal / tarch::la::norm2(normalizedNormal);
   validate ( tarch::la::equals(edge.getNormal(), normalizedNormal) );

   // Set displacements
   tarch::la::slice<2>(values,v0.getID()*2) = Vector2D(2.0);
   tarch::la::slice<2>(values,v1.getID()*2) = Vector2D(-2.0);

   // Apply displacements to  node coordinates
   modifyCoordinates.performAction(0.0, 0.0, 0.0, 0.0);

   // Validate coordinates of mesh after modifying it
   validate(tarch::la::equals(v0.getCoords(), Vector2D(2.0)));
   validate(tarch::la::equals(v1.getCoords(), Vector2D(-1.0)));
   normalizedNormal = Vector2D(-0.5, 0.5);
   normalizedNormal = normalizedNormal / tarch::la::norm2(normalizedNormal);
   validate(tarch::la::equals(edge.getNormal(), normalizedNormal));
}

void ModifyCoordinatesActionTest:: testSubtractFromCoordinates ()
{
   preciceTrace ( "testSubtractFromCoordinates()" );
   using namespace mesh;
   using utils::Vector2D;
   // Create geometryContext by faking a geometry but not using it to create
   // the mesh. The mesh is created by hand, such that references to the vertices
   // to be displaced are obtained.
   PtrMesh mesh ( new Mesh("Mesh", 2, false) );
   PtrData data = mesh->createData ( "test-data", 2 );
   int dataID = data->getID ();
   Vertex& v0 = mesh->createVertex ( Vector2D(0.0) );
   Vertex& v1 = mesh->createVertex ( Vector2D(1.0) );
   Edge& edge = mesh->createEdge ( v0, v1 );
   mesh->computeState();
   mesh->allocateDataValues ();
   utils::DynVector& values = data->values ();

   // Create ApplyDisplacementsMeshAction
   action::ModifyCoordinatesAction modifyCoordinates (
       action::ModifyCoordinatesAction::ALWAYS_PRIOR, dataID, mesh,
       action::ModifyCoordinatesAction::SUBTRACT_FROM_COORDINATES_MODE );
//   modifyCoordinates.loadMeshContext ( meshContext );

   // Validate coordinates of mesh before modifying it
   validate ( tarch::la::equals(v0.getCoords(), Vector2D(0.0)) );
   validate ( tarch::la::equals(v1.getCoords(), Vector2D(1.0)) );
   Vector2D normalizedNormal( Vector2D(0.5, -0.5) );
   normalizedNormal = normalizedNormal / tarch::la::norm2(normalizedNormal);
   validate ( tarch::la::equals(edge.getNormal(), normalizedNormal) );

   // Set displacements
   tarch::la::slice<2>(values,v0.getID()*2) = Vector2D(-2.0);
   tarch::la::slice<2>(values,v1.getID()*2) = Vector2D(2.0);

   // Apply displacements to  node coordinates
   modifyCoordinates.performAction(0.0, 0.0, 0.0, 0.0);

   // Validate coordinates of mesh after modifying it
   validate ( tarch::la::equals(v0.getCoords(), Vector2D(2.0)) );
   validate ( tarch::la::equals(v1.getCoords(), Vector2D(-1.0)) );
   normalizedNormal = Vector2D(-0.5, 0.5);
   normalizedNormal = normalizedNormal / tarch::la::norm2(normalizedNormal);
   validate ( tarch::la::equals(edge.getNormal(), normalizedNormal) );
}

}}} // namespace precice, action, tests

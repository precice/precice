// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "BalanceVertexPositionActionTest.hpp"
#include "action/BalanceVertexPositionAction.hpp"
#include "action/config/ActionConfiguration.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "utils/Globals.hpp"
#include "utils/Dimensions.hpp"
#include "io/ExportVTK.hpp"
#include "tarch/la/Scalar.h"
#include "geometry/Sphere.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::action::tests::BalanceVertexPositionActionTest)

namespace precice {
namespace action {
namespace tests {

tarch::logging::Log BalanceVertexPositionActionTest::
  _log ( "precice::action::tests::BalanceVertexPositionActionTest" );

BalanceVertexPositionActionTest:: BalanceVertexPositionActionTest ()
:
  TestCase ( "precice::action::tests::BalanceVertexPositionActionTest" )
{}

void BalanceVertexPositionActionTest:: run ()
{
  testMethod ( testSmoothCircle );
  testMethod ( testSmoothSphere );
  testMethod ( testSmoothHexahedron );
  testMethod ( testConfiguration );
}

void BalanceVertexPositionActionTest:: testSmoothCircle ()
{
  preciceTrace ( "testSmoothCircle()" );
  mesh::PtrMesh mesh ( new mesh::Mesh("Mesh", 2, false) );

  using namespace tarch::la;
  using utils::Vector2D;
  int nodes = 10;
  double circumference = 2.0 * PI; // radius = 1
  double h = circumference / (2.0 * nodes); // angle resolution
  mesh->createVertex ( Vector2D(1.0, 0.0) );
  for ( int i=1; i < nodes; i++ ){
    Vector2D coord ( std::cos((double)i * h), std::sin((double)i * h) );
    mesh::Vertex& v0 = mesh->vertices()[i-1];
    mesh::Vertex& v1 = mesh->createVertex ( coord );
    mesh->createEdge ( v0, v1 );
  }
  int size = (int)mesh->vertices().size();
  nodes = 20;
  h = circumference / (2.0 * nodes);
  for ( int i=0; i < nodes; i++ ){
    Vector2D coord ( std::cos(PI + (double)i * h), std::sin(PI + (double)i * h) );
    mesh::Vertex & v0 = mesh->vertices()[size-1 + i];
    mesh::Vertex & v1 = mesh->createVertex ( coord );
    mesh->createEdge ( v0, v1 );
  }
  size = (int)mesh->vertices().size();
  mesh->createEdge ( mesh->vertices()[size-1], mesh->vertices()[0] );
  mesh->computeState();
  io::ExportVTK exportVTK ( true );
  exportVTK.doExport ( "BalanceVertexPositionActionTest-testSmoothCircle-init.vtk",
                       *mesh );
  action::BalanceVertexPositionAction balance (
      action::Action::ALWAYS_PRIOR, mesh, 1e-10, 100 );

  double dummy = 0.0;
  balance.performAction(dummy, dummy, dummy, dummy);

  exportVTK.doExport("BalanceVertexPositionActionTest-testSmoothCircle-1.vtk",
                     *mesh);
}

void BalanceVertexPositionActionTest:: testSmoothSphere ()
{
  preciceTrace ( "testSmoothSphere()" );
  mesh::PtrMesh mesh ( new mesh::Mesh("Mesh", 3, false) );
  geometry::Sphere sphere ( utils::Vector3D(0.0), 0.1, 1.0 );
  sphere.create ( *mesh );

  io::ExportVTK exportVTK ( true );
  exportVTK.doExport ( "BalanceVertexPositionActionTest-testSmoothSphere-init.vtk",
                       *mesh );
  action::BalanceVertexPositionAction balance (
      action::Action::ALWAYS_PRIOR, mesh, 1e-10, 100 );
  double dummy = 0.0;
  balance.performAction(dummy, dummy, dummy, dummy);
  std::string name = "BalanceVertexPositionActionTest-testSmoothSphere-smooth.vtk";
  exportVTK.doExport(name, *mesh);
}

void BalanceVertexPositionActionTest:: testSmoothHexahedron()
{
  preciceTrace("testSmoothHexahedron()");
  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 3, false));

  using namespace tarch::la;
  using utils::Vector3D;
  using namespace mesh;

  Vertex& v000 = mesh->createVertex(Vector3D(-2.0, -1.0, -1.0));
  Vertex& v001 = mesh->createVertex(Vector3D(-2.0, -1.0,  1.0));
  Vertex& v010 = mesh->createVertex(Vector3D(-2.0,  1.0, -1.0));
  Vertex& v011 = mesh->createVertex(Vector3D(-2.0,  1.0,  1.0));
  Vertex& v100 = mesh->createVertex(Vector3D( 2.0, -1.0, -1.0));
  Vertex& v101 = mesh->createVertex(Vector3D( 2.0, -1.0,  1.0));
  Vertex& v110 = mesh->createVertex(Vector3D( 2.0,  1.0, -1.0));
  Vertex& v111 = mesh->createVertex(Vector3D( 2.0,  1.0,  1.0));

  Edge& e000to100 = mesh->createEdge(v000, v100);
  Edge& e010to110 = mesh->createEdge(v010, v110);
  Edge& e001to101 = mesh->createEdge(v001, v101);
  Edge& e011to111 = mesh->createEdge(v011, v111);

  Edge& e000to010 = mesh->createEdge(v000, v010);
  Edge& e100to110 = mesh->createEdge(v100, v110);
  Edge& e001to011 = mesh->createEdge(v001, v011);
  Edge& e101to111 = mesh->createEdge(v101, v111);

  Edge& e000to001 = mesh->createEdge(v000, v001);
  Edge& e100to101 = mesh->createEdge(v100, v101);
  Edge& e010to011 = mesh->createEdge(v010, v011);
  Edge& e110to111 = mesh->createEdge(v110, v111);
  Edge& e010to001 = mesh->createEdge(v010, v001);

  Edge& e100to111 = mesh->createEdge(v100, v111);
  Edge& e100to001 = mesh->createEdge(v100, v001);
  Edge& e110to011 = mesh->createEdge(v110, v011);
  Edge& e010to100 = mesh->createEdge(v010, v100);
  Edge& e001to111 = mesh->createEdge(v001, v111);

  mesh->createTriangle(e000to001, e010to001, e000to010); // x = 0
  mesh->createTriangle(e010to011, e010to001, e001to011);
  mesh->createTriangle(e100to101, e100to111, e101to111); // x = 1
  mesh->createTriangle(e100to110, e110to111, e100to111);

  mesh->createTriangle(e000to100, e100to001, e000to001); // y = 0
  mesh->createTriangle(e100to101, e001to101, e100to001);
  mesh->createTriangle(e010to110, e010to011, e110to011); // y = 1
  mesh->createTriangle(e110to111, e110to011, e011to111);

  mesh->createTriangle(e000to100, e000to010, e010to100); // z = 0
  mesh->createTriangle(e010to100, e010to110, e100to110);
  mesh->createTriangle(e001to101, e101to111, e001to111); // z = 1
  mesh->createTriangle(e001to011, e001to111, e011to111);

  mesh->computeState();
  io::ExportVTK exportVTK(true);
  exportVTK.doExport("BalanceVertexPositionActionTest-testSmoothHexahedron-init.vtk",
                     *mesh);
  action::BalanceVertexPositionAction balance(
      action::Action::ALWAYS_PRIOR, mesh, 1e-10, 100);

  double dummy = 0.0;
  balance.performAction(dummy, dummy, dummy, dummy);
  std::string name = "BalanceVertexPositionActionTest-testSmoothHexahedron-smooth.vtk";
  exportVTK.doExport(name, *mesh);
}

void BalanceVertexPositionActionTest:: testConfiguration()
{
  preciceTrace("testConfiguration()");
  {
    std::string filename = utils::Globals::getPathToSources() +
                           "/action/tests/BalanceVertexPositionActionTest-testConfiguration.xml";
    utils::XMLTag tag = utils::getRootTag();
    mesh::PtrDataConfiguration dataConfig(new mesh::DataConfiguration(tag));
    dataConfig->setDimensions(3);
    mesh::PtrMeshConfiguration meshConfig(new mesh::MeshConfiguration(tag,dataConfig));
    meshConfig->setDimensions(3);
    action::ActionConfiguration config(tag,meshConfig);
    utils::configure(tag, filename);
    validateEquals(config.actions().size(), 1);
    action::PtrAction action = config.actions().front();
    validate(action.get() != NULL);
  }
}

}}} // namespace precice, action, tests

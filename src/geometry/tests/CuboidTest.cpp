#include "CuboidTest.hpp"
#include "geometry/Cuboid.hpp"
#include "geometry/config/GeometryConfiguration.hpp"
#include "geometry/Geometry.hpp"
#include "geometry/SharedPointer.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "io/ExportVTK.hpp"
#include "utils/Parallel.hpp"
#include "utils/xml/XMLTag.hpp"
#include "utils/Globals.hpp"
#include <iostream>

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::geometry::tests::CuboidTest)

namespace precice {
namespace geometry {
namespace tests {

tarch::logging::Log CuboidTest:: _log ( "precice::geometry::tests::CuboidTest" );

CuboidTest:: CuboidTest ()
:
  TestCase ("geometry::CuboidTest")
{}

void CuboidTest:: run ()
{
  PRECICE_MASTER_ONLY {
    testMethod ( testCreation );
    testMethod ( testConfiguration );
    testMethod ( testSubIDs3D );
    testMethod ( testSubIDs2D );
  }
}

void CuboidTest:: testCreation ()
{
  preciceTrace ( "testCreation()" );

  for ( int dim=2; dim <= 3; dim++ ){
    utils::DynVector offset(dim, 0.0);
    utils::DynVector length(dim, 0.0);
    for (int i=0; i < dim; i++){
      offset(i) = static_cast<double>(i);
      length(i) = static_cast<double>(dim - i);
    }
    bool flipNormals = false;
    mesh::Mesh mesh("test-cuboid", dim, flipNormals);
    geometry::Cuboid cuboid(offset, 0.5, length);
    cuboid.create(mesh);
    io::ExportVTK exportVTK(true);
    std::ostringstream name;
    name << "geometry-CuboidTest-dim" << dim;
    exportVTK.doExport(name.str(), mesh);
  }
}

void CuboidTest:: testConfiguration()
{
  preciceTrace("testConfiguration()");
  using namespace boost;

  for (int dim=2; dim <= 3; dim++){
    std::string xmlFilename = utils::Globals::getPathToSources() + "/geometry/tests/";
    preciceDebug("dim = " << dim);
    if (dim == 2){
      xmlFilename += "cuboid2d.xml";
    }
    else {
      assertion ( dim == 3 );
      xmlFilename += "cuboid3d.xml";
    }
    utils::XMLTag tag = utils::getRootTag();
    mesh::PtrDataConfiguration dataConfig ( new mesh::DataConfiguration(tag) );
    dataConfig->setDimensions(dim);
    mesh::PtrMeshConfiguration meshConfig ( new mesh::MeshConfiguration(tag, dataConfig) );
    meshConfig->setDimensions(dim);
    GeometryConfiguration geoConfig ( tag, meshConfig );
    geoConfig.setDimensions(dim);

    utils::configure ( tag, xmlFilename );
    //validate ( dataConfig->isValid() );
    //validate ( meshConfig->isValid() );
    meshConfig->setMeshSubIDs();
    validateEquals ( meshConfig->meshes().size(), 1 );
    //validate ( geoConfig.isValid() );
    validateEquals ( geoConfig.geometries().size(), 1 );

    mesh::PtrMesh mesh = meshConfig->meshes()[0];
    std::map<std::string,int> nameIDs = mesh->getNameIDPairs ();
    if (dim == 2){
      validateEquals ( nameIDs.size(), 5 );
    }
    else {
      assertion ( dim == 3 );
      validateEquals ( nameIDs.size(), 7 );
    }
    geometry::PtrGeometry geo = geoConfig.geometries()[0];
    geo->create ( *mesh );
  }
}

void CuboidTest:: testSubIDs2D ()
{
  preciceTrace ( "testSubIDs2D" );
  using utils::Vector2D;
  mesh::Mesh::resetGeometryIDsGlobally ();
  bool flipNormals = false;
  mesh::Mesh mesh ( "test-cuboid", 2, flipNormals );
  utils::DynVector offset ( 2, 0.0 );
  double discretizationWidth = 1.0;
  utils::DynVector sidelengths ( 2, 1.0 );
  Cuboid cuboid ( offset, discretizationWidth, sidelengths );
  mesh.setSubID ( "side-1" );
  mesh.setSubID ( "side-2" );
  cuboid.create ( mesh );
  int id = mesh.getID ( "test-cuboid" );
  int idSide1 = mesh.getID ( "test-cuboid-side-1" );
  int idSide2 = mesh.getID ( "test-cuboid-side-2" );
  validateEquals ( mesh.vertices().size(), 4 );
  validateEquals ( mesh.edges().size(), 4 );
  validateEquals ( mesh.triangles().size(), 0 );

  // Test geometry IDs of vertices
  using tarch::la::equals;
  for (mesh::Vertex& vertex : mesh.vertices()) {
    std::vector<int> geometryIDs;
    vertex.getProperties ( vertex.INDEX_GEOMETRY_ID, geometryIDs );
    if ( equals(vertex.getCoords(), Vector2D(0.0, 0.0)) ) {
      validateEquals ( geometryIDs.size(), 2 );
      validate ( utils::contained(id, geometryIDs) );
      validate ( utils::contained(idSide2, geometryIDs) );
    }
    else if ( equals(vertex.getCoords(), Vector2D(0.0, 1.0)) ) {
      validateEquals ( geometryIDs.size(), 1 );
      validate ( utils::contained(id, geometryIDs) );
    }
    else if ( equals(vertex.getCoords(), Vector2D(1.0, 0.0)) ) {
      validateEquals ( geometryIDs.size(), 3 );
      validate ( utils::contained(id, geometryIDs) );
      validate ( utils::contained(idSide1, geometryIDs) );
      validate ( utils::contained(idSide2, geometryIDs) );
    }
    else if ( equals(vertex.getCoords(), Vector2D(1.0, 1.0)) ) {
      validateEquals ( geometryIDs.size(), 2 );
      validate ( utils::contained(id, geometryIDs) );
      validate ( utils::contained(idSide1, geometryIDs) );
    }
    else {
      preciceDebug ( "Wrong coords = " << vertex.getCoords() );
      validate ( false );
    }
  }

  // Test geometry IDs of edges
  for (mesh::Edge& edge : mesh.edges()) {
    std::vector<int> geometryIDs;
    edge.getProperties ( edge.INDEX_GEOMETRY_ID, geometryIDs );
    if ( equals(edge.getCenter(), Vector2D(0.5, 0.0)) ) {
      validateEquals ( geometryIDs.size(), 2 );
      validate ( utils::contained(idSide2, geometryIDs) );
    }
    else if ( equals(edge.getCenter(), Vector2D(0.0, 0.5)) ) {
      validateEquals ( geometryIDs.size(), 1 );
      validate ( utils::contained(id, geometryIDs) );
    }
    else if ( equals(edge.getCenter(), Vector2D(1.0, 0.5)) ) {
      validateEquals ( geometryIDs.size(), 2 );
      validate ( utils::contained(idSide1, geometryIDs) );
    }
    else if ( equals(edge.getCenter(), Vector2D(0.5, 1.0)) ) {
      validateEquals ( geometryIDs.size(), 1 );
      validate ( utils::contained(id, geometryIDs) );
    }
    else {
      preciceDebug ( "Wrong center = " << edge.getCenter() );
      validate ( false );
    }
  }
}

void CuboidTest:: testSubIDs3D ()
{
  preciceTrace ( "testSubIDs3D" );
  using utils::contained;
  using namespace tarch::la;
  using utils::Vector3D;
  mesh::Mesh::resetGeometryIDsGlobally ();

  bool flipNormals = false;
  mesh::Mesh mesh ( "test-cuboid", 3, flipNormals );
  utils::DynVector offset ( 3, 0.0 );
  double discretizationWidth = 1.0;
  utils::DynVector sidelengths ( 3, 1.0 );
  Cuboid cuboid ( offset, discretizationWidth, sidelengths );
  int baseID = mesh.getID ( "test-cuboid" );
  mesh.setSubID ( "side-1" );
  int idSide1 = mesh.getID ( "test-cuboid-side-1" );
  mesh.setSubID ( "side-2" );
  int idSide2 = mesh.getID ( "test-cuboid-side-2" );
  mesh.setSubID ( "side-5" );
  int idSide5 = mesh.getID ( "test-cuboid-side-5" );
  cuboid.create ( mesh );
  validateEquals ( mesh.vertices().size(), 8 );
  validateEquals ( mesh.edges().size(), 18 );
  validateEquals ( mesh.triangles().size(), 12 );

  // Test geometry IDs of vertices
  for (mesh::Vertex& vertex : mesh.vertices()) {
    std::vector<int> geometryIDs;
    vertex.getProperties ( vertex.INDEX_GEOMETRY_ID, geometryIDs );
    if ( equals(vertex.getCoords(), Vector3D(0.0, 0.0, 0.0)) ) {
      validateEquals ( geometryIDs.size(), 2 );
      validate ( contained(baseID, geometryIDs) );
      validate ( contained(idSide2, geometryIDs) );
    }
    else if ( equals(vertex.getCoords(),Vector3D(0.0, 0.0, 1.0)) ) {
      validateEquals ( geometryIDs.size(), 3 );
      validate ( contained(baseID, geometryIDs) );
      validate ( contained(idSide2, geometryIDs) );
      validate ( contained(idSide5, geometryIDs) );
    }
    else if ( equals(vertex.getCoords(),Vector3D(0.0, 1.0, 0.0)) ) {
      validateEquals ( geometryIDs.size(), 1 );
      validate ( contained(baseID, geometryIDs) );
    }
    else if ( equals(vertex.getCoords(),Vector3D(0.0, 1.0, 1.0)) ) {
      validateEquals ( geometryIDs.size(), 2 );
      validate ( contained(baseID, geometryIDs) );
      validate ( contained(idSide5, geometryIDs) );
    }
    else if ( equals(vertex.getCoords(),Vector3D(1.0, 0.0, 0.0)) ) {
      validateEquals ( geometryIDs.size(), 3 );
      validate ( contained(baseID, geometryIDs) );
      validate ( contained(idSide1, geometryIDs) );
      validate ( contained(idSide2, geometryIDs) );
    }
    else if ( equals(vertex.getCoords(),Vector3D(1.0, 0.0, 1.0)) ) {
      validateEquals ( geometryIDs.size(), 4 );
      validate ( contained(baseID, geometryIDs) );
      validate ( contained(idSide1, geometryIDs) );
      validate ( contained(idSide2, geometryIDs) );
      validate ( contained(idSide5, geometryIDs) );
    }
    else if ( equals(vertex.getCoords(),Vector3D(1.0, 1.0, 0.0)) ) {
      validateEquals ( geometryIDs.size(), 2 );
      validate ( contained(baseID, geometryIDs) );
      validate ( contained(idSide1, geometryIDs) );
    }
    else if ( equals(vertex.getCoords(),Vector3D(1.0, 1.0, 1.0)) ) {
      validateEquals ( geometryIDs.size(), 3 );
      validate ( contained(baseID, geometryIDs) );
      validate ( contained(idSide1, geometryIDs) );
      validate ( contained(idSide5, geometryIDs) );
    }
    else {
      preciceDebug ( "Wrong coords = " << vertex.getCoords() );
      validate ( false );
    }

    // Test geometryIDs of edges
    for (mesh::Edge & edge : mesh.edges()) {
      std::vector<int> geometryIDs;
      edge.getProperties ( edge.INDEX_GEOMETRY_ID, geometryIDs );
      if ( equals(edge.getCenter(),Vector3D(0.5, 0.0, 0.0)) ) {
        validateEquals ( geometryIDs.size(), 2 );
        validate ( contained(baseID, geometryIDs) );
        validate ( contained(idSide2, geometryIDs) );
      }
      else if ( equals(edge.getCenter(),Vector3D(0.0, 0.5, 0.0)) ) {
        validateEquals ( geometryIDs.size(), 1 );
        validate ( contained(baseID, geometryIDs) );
      }
      else if ( equals(edge.getCenter(),Vector3D(0.0, 0.0, 0.5)) ) {
        validateEquals ( geometryIDs.size(), 2 );
        validate ( contained(baseID, geometryIDs) );
        validate ( contained(idSide2, geometryIDs) );
      }
      else if ( equals(edge.getCenter(),Vector3D(0.0, 0.5, 0.5)) ) {
        validateEquals ( geometryIDs.size(), 1 );
        validate ( contained(baseID, geometryIDs) );
      }
      else if ( equals(edge.getCenter(),Vector3D(0.5, 0.5, 0.0)) ) {
        validateEquals ( geometryIDs.size(), 1 );
        validate ( contained(baseID, geometryIDs) );
      }
      else if ( equals(edge.getCenter(),Vector3D(0.5, 0.0, 0.5)) ) {
        validateEquals ( geometryIDs.size(), 2 );
        validate ( contained(baseID, geometryIDs) );
        validate ( contained(idSide2, geometryIDs) );
      }
      else if ( equals(edge.getCenter(),Vector3D(0.5, 1.0, 0.0)) ) {
        validateEquals ( geometryIDs.size(), 1 );
        validate ( contained(baseID, geometryIDs) );
      }
      else if ( equals(edge.getCenter(),Vector3D(0.5, 1.0, 1.0)) ) {
        validateEquals ( geometryIDs.size(), 2 );
        validate ( contained(baseID, geometryIDs) );
        validate ( contained(idSide5, geometryIDs) );
      }
      else if ( equals(edge.getCenter(),Vector3D(0.5, 1.0, 0.5)) ) {
        validateEquals ( geometryIDs.size(), 1 );
        validate ( contained(baseID, geometryIDs) );
      }
      else if ( equals(edge.getCenter(),Vector3D(0.0, 1.0, 0.5)) ) {
        validateEquals ( geometryIDs.size(), 1 );
        validate ( contained(baseID, geometryIDs) );
      }
      else if ( equals(edge.getCenter(),Vector3D(0.0, 0.5, 1.0)) ) {
        validateEquals ( geometryIDs.size(), 2 );
        validate ( contained(baseID, geometryIDs) );
        validate ( contained(idSide5, geometryIDs) );
      }
      else if ( equals(edge.getCenter(),Vector3D(0.5, 0.0, 1.0)) ) {
        validateEquals ( geometryIDs.size(), 3 );
        validate ( contained(baseID, geometryIDs) );
        validate ( contained(idSide2, geometryIDs) );
        validate ( contained(idSide5, geometryIDs) );
      }
      else if ( equals(edge.getCenter(),Vector3D(0.5, 0.5, 1.0)) ) {
        validateEquals ( geometryIDs.size(), 2 );
        validate ( contained(baseID, geometryIDs) );
        validate ( contained(idSide5, geometryIDs) );
      }
      else if ( equals(edge.getCenter(),Vector3D(1.0, 0.0, 0.5)) ) {
        validateEquals ( geometryIDs.size(), 3 );
        validate ( contained(baseID, geometryIDs) );
        validate ( contained(idSide2, geometryIDs) );
        validate ( contained(idSide1, geometryIDs) );
      }
      else if ( equals(edge.getCenter(),Vector3D(1.0, 0.5, 0.0)) ) {
        validateEquals ( geometryIDs.size(), 2 );
        validate ( contained(baseID, geometryIDs) );
        validate ( contained(idSide1, geometryIDs) );
      }
      else if ( equals(edge.getCenter(),Vector3D(1.0, 1.0, 0.5)) ) {
        validateEquals ( geometryIDs.size(), 2 );
        validate ( contained(baseID, geometryIDs) );
        validate ( contained(idSide1, geometryIDs) );
      }
      else if ( equals(edge.getCenter(),Vector3D(1.0, 0.5, 1.0)) ) {
        validateEquals ( geometryIDs.size(), 3 );
        validate ( contained(baseID, geometryIDs) );
        validate ( contained(idSide1, geometryIDs) );
        validate ( contained(idSide5, geometryIDs) );
      }
      else if ( equals(edge.getCenter(),Vector3D(1.0, 0.5, 0.5)) ) {
        validateEquals ( geometryIDs.size(), 2 );
        validate ( contained(baseID, geometryIDs) );
        validate ( contained(idSide1, geometryIDs) );
      }
      else {
        preciceDebug ( "Wrong coords = " << edge.getCenter() );
        validate ( false );
      }
    }

    // Test geometryIDs of edges
    for (mesh::Triangle & triangle : mesh.triangles()) {
      std::vector<int> geometryIDs;
      triangle.getProperties ( triangle.INDEX_GEOMETRY_ID, geometryIDs );
      if ( tarch::la::equals(triangle.getCenter()(0), 0.0) ) {
        validateEquals ( geometryIDs.size(), 1 );
        validate ( contained(baseID, geometryIDs) );
      }
      else if ( equals(triangle.getCenter()(0), 1.0) ) {
        validateEquals ( geometryIDs.size(), 2 );
        validate ( contained(baseID, geometryIDs) );
        validate ( contained(idSide1, geometryIDs) );
      }
      else if ( equals(triangle.getCenter()(1), 0.0) ) {
        validateEquals ( geometryIDs.size(), 2 );
        validate ( contained(baseID, geometryIDs) );
        validate ( contained(idSide2, geometryIDs) );
      }
      else if ( equals(triangle.getCenter()(1), 1.0) ) {
        validateEquals ( geometryIDs.size(), 1 );
        validate ( contained(baseID, geometryIDs) );
      }
      else if ( equals(triangle.getCenter()(2), 0.0) ) {
        validateEquals ( geometryIDs.size(), 1 );
        validate ( contained(baseID, geometryIDs) );
      }
      else if ( equals(triangle.getCenter()(2), 1.0) ) {
        validateEquals ( geometryIDs.size(), 2 );
        validate ( contained(baseID, geometryIDs) );
        validate ( contained(idSide5, geometryIDs) );
      }
      else {
        preciceDebug ( "Wrong triangle center = " << triangle.getCenter() );
        validate ( false );
      }
    }
  }
}

}}} // namespace precice, geometry, tests

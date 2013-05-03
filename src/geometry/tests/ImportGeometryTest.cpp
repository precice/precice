// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ImportGeometryTest.hpp"
#include "geometry/Geometry.hpp"
#include "geometry/SharedPointer.hpp"
#include "geometry/config/GeometryConfiguration.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "mesh/Mesh.hpp"
#include "com/config/CommunicationConfiguration.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"
#include "utils/Dimensions.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::geometry::tests::ImportGeometryTest)

namespace precice {
namespace geometry {
namespace tests {

tarch::logging::Log ImportGeometryTest::
  _log ( "precice::geometry::tests::ImportGeometryTest" );

ImportGeometryTest:: ImportGeometryTest ()
:
  TestCase ( "geometry::tests::ImportGeometryTest" )
{}

void ImportGeometryTest:: run ()
{
# ifndef PRECICE_NO_SPIRIT2
  PRECICE_MASTER_ONLY {
    testMethod ( testImportVRMLConfig );
  }
# endif // not PRECICE_NO_SPIRIT2
}

void ImportGeometryTest:: testImportVRMLConfig ()
{
  preciceTrace ( "testImportVRMLConfig" );
  using namespace boost;

  std::string xmlFilename ( utils::Globals::getPathToSources() +
    "/geometry/tests/import-vrml-config.xml" );
  utils::XMLTag tag = utils::getRootTag();
  mesh::PtrDataConfiguration dataConfig ( new mesh::DataConfiguration(tag) );
  dataConfig->setDimensions(3);
  mesh::PtrMeshConfiguration meshConfig ( new mesh::MeshConfiguration(tag, dataConfig) );
  meshConfig->setDimensions(3);
  GeometryConfiguration geoConfig ( tag, meshConfig );
  geoConfig.setDimensions(3);

  utils::configure ( tag, xmlFilename );
  //validate ( dataConfig->isValid() );
  //validate ( meshConfig->isValid() );
  meshConfig->setMeshSubIDs();
  validateEquals ( meshConfig->meshes().size(), 1 );
  //validate ( geoConfig.isValid() );
  validateEquals ( geoConfig.geometries().size(), 1 );

  PtrGeometry geo = geoConfig.getGeometry ( "VRMLGeometry" );
  mesh::PtrMesh mesh = meshConfig->meshes()[0];
  validateEquals ( mesh->getName(), std::string("VRMLGeometry") );
  validate ( tarch::la::equals(geo->getOffset(), utils::Vector3D(0.0)) );
  validateEquals ( mesh->data().size(), 1 );
}

}}} // namespace precice, geometry, tests

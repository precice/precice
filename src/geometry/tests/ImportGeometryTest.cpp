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

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::geometry::tests::ImportGeometryTest)

namespace precice {
namespace geometry {
namespace tests {

logging::Logger ImportGeometryTest:: _log("precice::geometry::tests::ImportGeometryTest");

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
  TRACE();
  std::string xmlFilename ( utils::Globals::getPathToSources() + "/geometry/tests/import-vrml-config.xml" );
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
  validate ( math::equals(geo->getOffset(), Eigen::Vector3d::Zero()) );
  validateEquals ( mesh->data().size(), 1 );
}

}}} // namespace precice, geometry, tests

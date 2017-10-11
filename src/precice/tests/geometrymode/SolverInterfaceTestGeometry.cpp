#include "SolverInterfaceTestGeometry.hpp"
#include "precice/Constants.hpp"
#include "precice/config/SolverInterfaceConfiguration.hpp"
#include "precice/impl/MeshContext.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"
#include "precice/config/Configuration.hpp"
#include "geometry/Cuboid.hpp"
#include "mesh/Mesh.hpp"
#include "query/FindClosest.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"
#include "io/ExportVTK.hpp"
#include <vector>
#include <set>
#include <algorithm>
#include <fstream>

#include "tarch/tests/TestCaseFactory.h"
//registerIntegrationTest(precice::tests::SolverInterfaceTestGeometry)

namespace precice {
namespace tests {

logging::Logger SolverInterfaceTestGeometry::
    _log ( "precice::tests::SolverInterfaceTestGeometry" );

SolverInterfaceTestGeometry:: SolverInterfaceTestGeometry ()
:
  TestCase ( "tests::SolverInterfaceTestGeometry" ),
  _pathToTests (),
  _geoID ( -1 )
{}

void SolverInterfaceTestGeometry:: setUp ()
{
  _pathToTests = utils::Globals::getPathToSources() + "/precice/tests/geometrymode/";
}

void SolverInterfaceTestGeometry:: run()
{
  TRACE();
  typedef utils::Parallel Par;
  Par::Communicator comm = Par::getRestrictedCommunicator({0});
  PRECICE_MASTER_ONLY {
    TRACE();
    Par::setGlobalCommunicator(comm);
    testMethod(testConfiguration);
    testMethod(testSearchQuery);
    testMethod(testVoxelQuery);
    testMethod(testDataActions);
    testMethod(testVoxelQueryMultipleGeometryIDs);
    testMethod(testVoxelQueryDFGChannel);
    testMethod(testVoxelQueryFSIChannel);
    testMethod(testVoxelQueryChannelFour);
    testMethod(testVoxelQueryEpsBox);
    testMethod(testConservativeStationaryDataMapping);
    testMethod(testMappingRBF);
    testMethod(testCustomGeometryCreation);
#   ifndef PRECICE_NO_PYTHON
    testMethod(testPinelli);
#   endif // not PRECICE_NO_PYTHON
#   ifndef PRECICE_NO_SPIRIT2
//    testMethod(testBug);   // excluded to speed up integration tests
//    testMethod(testBug2);  // excluded to speed up integration tests
#   endif // not PRECICE_NO_SPIRIT2
    testMethod(testBug3);
    testMethod(testBug4);
    testMethod(testBug5);
    testMethod(testUpdateSpacetree);
    testMethod(testMultipleMeshSpacetree);
    Par::setGlobalCommunicator(Par::getCommunicatorWorld());
  }
}

void SolverInterfaceTestGeometry:: configureSolverInterface
(
  const std::string& configFilename,
  SolverInterface&   interface )
{
  TRACE(configFilename);
  mesh::Mesh::resetGeometryIDsGlobally();
  mesh::Data::resetDataCount();
  impl::Participant::resetParticipantCount();
  config::Configuration config;
  utils::configure ( config.getXMLTag(), configFilename );
  //validate ( config.isValid() );
  interface._impl->configure ( config.getSolverInterfaceConfiguration() );
}

void SolverInterfaceTestGeometry:: testConfiguration()
{
  TRACE();
  mesh::Mesh::resetGeometryIDsGlobally ();
  DEBUG ( "Test 2D configuration");
  { // 2D
    SolverInterface geoInterface ( "TestAccessor", 0, 1 );
    configureSolverInterface (
        _pathToTests + "configuration2D.xml", geoInterface );
    validateEquals ( geoInterface.getDimensions(), 2 );
    geoInterface.initialize ();
    geoInterface.exportMesh ( "testConfiguration2D" );
  }
  DEBUG ( "Test 3D configuration");
  { // 3D
    SolverInterface geoInterface ( "TestAccessor", 0, 1 );
    configureSolverInterface (
        _pathToTests + "configuration3D.xml", geoInterface );
    validateEquals ( geoInterface.getDimensions(), 3 );
    geoInterface.initialize ();
    geoInterface.exportMesh ( "testConfiguration3D" );
  }
}

void SolverInterfaceTestGeometry:: testSearchQuery()
{
  TRACE();
  for ( int dim=2; dim <= 3; dim++ ) {
    SolverInterface geoInterface ( "TestAccessor", 0, 1 );
    if (dim == 2){
      configureSolverInterface ( _pathToTests + "2D.xml",
                                 geoInterface );
    }
    else {
      configureSolverInterface ( _pathToTests + "3D.xml",
                                 geoInterface );
    }
    
    geoInterface.initialize();

    int meshID = geoInterface.getMeshID("SolverMesh");
    Eigen::VectorXd pos = Eigen::VectorXd::Zero(dim);
    pos.setConstant(50);
    geoInterface.setMeshVertex(meshID, pos.data());

    _geoID = geoInterface.getMeshID("itest-cuboid");

    std::set<int> ids;
    pos.setConstant(0);
    ClosestMesh closest = geoInterface.inquireClosestMesh (pos.data(), ids);
    validateEquals ( closest.meshIDs().size(), 1 );
    validateEquals ( closest.meshIDs()[0], _geoID );
    validateEquals ( closest.position(), constants::positionOutsideOfGeometry() );

    pos.setConstant(4);
    closest = geoInterface.inquireClosestMesh (pos.data(), ids);
    validateEquals ( closest.meshIDs().size(), 1 );
    validateEquals ( closest.meshIDs()[0], _geoID );
    validateEquals ( closest.position(), constants::positionOutsideOfGeometry() );

    pos.setConstant(5);
    closest = geoInterface.inquireClosestMesh (pos.data(), ids);
    validateEquals ( closest.meshIDs().size(), 1 );
    //validateEquals ( closest.meshIDs()[0], _geoID );  // why does this not work
    validateEquals ( closest.position(), constants::positionOnGeometry() );

    pos.setConstant(6);
    closest = geoInterface.inquireClosestMesh ( pos.data(), ids );
    validateEquals ( closest.meshIDs().size(), 1 );
    validateEquals ( closest.meshIDs()[0], _geoID );
    validateEquals ( closest.position(), constants::positionInsideOfGeometry() );

    if (dim == 2){
      geoInterface.exportMesh (
          "SolverInterfaceTestGeometry-testSearchQuery-testCreateInterface-2D" );
    }
    else {
      geoInterface.exportMesh (
          "SolverInterfaceTestGeometry-testSearchQuery-testCreateInterface-3D" );
    }
  }
}

void SolverInterfaceTestGeometry:: testVoxelQuery()
{
  TRACE();
  
  for ( int dim=2; dim <= 3; dim++ ){
    SolverInterface geoInterface ( "TestAccessor", 0, 1 );
    if (dim == 2){
      configureSolverInterface ( _pathToTests + "2D.xml",
                                 geoInterface );
    }
    else {
      configureSolverInterface ( _pathToTests + "3D.xml",
                                 geoInterface );
    }

    geoInterface.initialize();

    int meshID = geoInterface.getMeshID("SolverMesh");
    Eigen::VectorXd posVertex = Eigen::VectorXd::Constant(dim,50.0);
    geoInterface.setMeshVertex(meshID, posVertex.data());


    _geoID = geoInterface.getMeshID ("itest-cuboid");

    // Voxel completely contained in center
    std::set<int> ids;
    bool include = false;
    Eigen::VectorXd center = Eigen::VectorXd::Zero(dim);
    Eigen::VectorXd h = Eigen::VectorXd::Constant(dim, 0.1);
    VoxelPosition pos = geoInterface.inquireVoxelPosition (
      center.data(), h.data(), include, ids );
    validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

    if (dim == 2){
      // Voxels in corners
      center << -4.0, -4.0;
      h << 1.0, 1.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids);
      validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
      center << 4.0, -4.0;
      h << 1.0, 1.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
      center << -4.0, 4.0;
      h << 1.0, 1.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
      center << 4.0, 4.0;
      h << 1.0, 1.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

      // Voxels on sides (outside of geometry)
      center << -4.0, 0.0;
      h << 1.0, 1.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
      center << 4.0, 0.0;
      h << 1.0, 1.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
      center << -4.0, 4.0;
      h << 1.0, 1.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
      center << 4.0, 4.0;
      h << 1.0, 1.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
      center << -4.0, 0.0;
      h << 1.0 + (math::NUMERICAL_ZERO_DIFFERENCE/2.0), 1.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

      // Voxels on sides (inside of geometry)
      center << -6.0, 0.0;
      h << 1.0, 1.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
      center << 6.0, 0.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
      center << -6.0, 6.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
      center << 6.0, 6.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionInsideOfGeometry() );

      center << 5.0, 5.0;
      h << 5.0, 5.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.meshIDs().size(), 1 );
      validateEquals ( pos.meshIDs()[0], _geoID );
      geoInterface.exportMesh ("SolverInterfaceTestGeometry-testVoxelQuery-2D");
    }
    else { // 3D
      // Voxels in corners (outside of geometry)
      center << -4.0, -4.0, -4.0;
      h << 1.0, 1.0, 1.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
      center << 4.0, -4.0, -4.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
      center << -4.0, 4.0, -4.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
      center << 4.0, 4.0, -4.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
      center << -4.0, -4.0, 4.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
      center << 4.0, -4.0, 4.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
      center << -4.0, 4.0, 4.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
      center << 4.0, 4.0, 4.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

      // Voxels in corners (inside of geometry)
      center << -6.0, -6.0, -6.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
      center << 6.0, -6.0, -6.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
      center << -6.0, 6.0, -6.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
      center << 6.0, 6.0, -6.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
      center << -6.0, -6.0, 6.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
      center << 6.0, -6.0, 6.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
      center << -6.0, 6.0, 6.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
      center << 6.0, 6.0, 6.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionInsideOfGeometry() );

      // Voxels on sides (outside of geometry)
      center << -4.0, 0.0, 0.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionOutsideOfGeometry());
      center << 4.0, 0.0, 0.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
      center << 0.0, -4.0, 0.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
      center << 0.0, 4.0, 0.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
      center << 0.0, 0.0, -4.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
      center << 0.0, 0.0, 4.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

      // Voxels on sides (inside of geometry)
      center << -6.0, 0.0, 0.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
      center << 6.0, 0.0, 0.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
      center << 0.0, -6.0, 0.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
      center << 0.0, 6.0, 0.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
      center << 0.0, 0.0, -6.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
      center << 0.0, 0.0, 6.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.position(), constants::positionInsideOfGeometry() );

      // Voxels on geometry
      center << 5.0, 5.0, 5.0;
      h << 5.0, 5.0, 5.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.meshIDs().size(), 1 );
      validateEquals ( pos.meshIDs()[0], _geoID );

      center << -5.0, -5.0, -5.0;
      pos = geoInterface.inquireVoxelPosition (center.data(), h.data(), include, ids );
      validateEquals ( pos.meshIDs().size(), 1 );
      validateEquals ( pos.meshIDs()[0], _geoID );
      geoInterface.exportMesh ("SolverInterfaceTestGeometry-testVoxelQuery-3D");
    }
  }
}

void SolverInterfaceTestGeometry:: testDataActions()
{
  TRACE();
  SolverInterface geo("Accessor", 0, 1);
  configureSolverInterface(_pathToTests + "testDataActions.xml", geo);
  impl::SolverInterfaceImpl* impl = geo._impl.get();
  int meshID = geo.getMeshID("Box");
  int dataID = geo.getDataID("VectorData", meshID);
  geo.initialize();
  mesh::PtrMesh mesh = impl->_accessor->meshContext(meshID).mesh;
  std::vector<Eigen::Vector3d> coords ( mesh->vertices().size() );
  for ( size_t i=0; i < mesh->vertices().size(); i++ ){
    coords[i] = mesh->vertices()[i].getCoords();
  }
  mesh::PtrData data = mesh->data ( dataID );
  data->values() = Eigen::VectorXd::Constant(data->values().size(), 1.0);
  //assign(data->values()) = 1.0;

  geo.advance ( 1.0 );

  for ( size_t i=0; i < mesh->vertices().size(); i++ ){
    validate ( math::equals(coords[i] + Eigen::VectorXd::Constant(coords[i].size(), 1.0),
                            mesh->vertices()[i].getCoords()) );
  }

  geo.advance ( 1.0 );

  for ( size_t i=0; i < mesh->vertices().size(); i++ ){
    validate ( math::equals(coords[i] + Eigen::VectorXd::Constant(coords[i].size(), 2.0),
                            mesh->vertices()[i].getCoords()) );
  }

  geo.finalize ();
}

void SolverInterfaceTestGeometry:: testVoxelQueryMultipleGeometryIDs()
{
  TRACE();
  using Eigen::Vector2d;
  SolverInterface geoInterface ( "TestAccessor", 0, 1 );
  configureSolverInterface (
      _pathToTests + "testMultipleGeometryIDs.xml",
      geoInterface );
  validateEquals ( geoInterface.getDimensions(), 2 );
  geoInterface.initialize();

  int idQuad0Base = geoInterface.getMeshID ( "Quad0" );
  int idQuad0Side1 = geoInterface.getMeshID ( "Quad0-side-1" );
  int idQuad0Side2 = geoInterface.getMeshID ( "Quad0-side-2" );
  int idQuad0Side3 = geoInterface.getMeshID ( "Quad0-side-3" );

  int idQuad1Base = geoInterface.getMeshID ( "Quad1" );
  int idQuad1Side2 = geoInterface.getMeshID ( "Quad1-side-2" );

  geoInterface.exportMesh ( "testVoxelQueryMultipleGeometryIDs" );

  std::set<int> ids;
  VoxelPosition pos = geoInterface.inquireVoxelPosition (
    Vector2d(0, 0).data(), Vector2d(1.0, 1.0).data(), false, ids );
  validateEquals ( pos.position(), constants::positionOnGeometry() );
  validateEquals ( pos.meshIDs().size(), 2 );
  if ( pos.meshIDs()[0] == idQuad0Base ) {
    validateEquals ( pos.meshIDs()[1], idQuad0Side2 );
  }
  else {
    validateEquals ( pos.meshIDs()[0], idQuad0Side2 );
    validateEquals ( pos.meshIDs()[1], idQuad0Base );
  }

  pos = geoInterface.inquireVoxelPosition ( Vector2d(5.0, 5.0).data(),
                                            Vector2d(1.5, 1.5).data(),
                                            false, ids );
  validateEquals ( pos.position(), constants::positionOnGeometry() );
  validateEquals ( pos.meshIDs().size(), 5 );

  for ( size_t i=0; i < 4; i++ ){
    if ( (pos.meshIDs()[i] == idQuad0Base)
      || (pos.meshIDs()[i] == idQuad0Side1)
      || (pos.meshIDs()[i] == idQuad0Side3)
      || (pos.meshIDs()[i] == idQuad1Base)
      || (pos.meshIDs()[i] == idQuad1Side2) )
    {
    }
    else {
      validate ( false );
    }
  }
}

void SolverInterfaceTestGeometry:: testVoxelQueryDFGChannel()
{
  TRACE();
  SolverInterface geoInterface ( "TestAccessor", 0, 1 );
  configureSolverInterface ( _pathToTests + "dfgchannel.xml", geoInterface );
  validateEquals ( geoInterface.getDimensions(), 2 );
  geoInterface.initialize();
  _geoID = geoInterface.getMeshID ("itest-dfg-cuboid");
  Eigen::Vector2d inquiryCenter (0.41/9.0, (46.0 * 0.41) / 9.0);
  Eigen::Vector2d inquiryHalflengths (0.41/9.0, 0.41/9.0);
  std::set<int> ids;
  VoxelPosition pos = geoInterface.inquireVoxelPosition (
    inquiryCenter.data(), inquiryHalflengths.data(), false, ids );
  validateEquals ( pos.position(), constants::positionInsideOfGeometry() );

  inquiryCenter << 0.0455556, 2.09556;
  inquiryHalflengths << 0.0455556, 0.0455556;
  pos = geoInterface.inquireVoxelPosition ( inquiryCenter.data(), inquiryHalflengths.data(),
    false, ids );
  validateEquals ( pos.position(), constants::positionInsideOfGeometry() );

  inquiryCenter << 0.41/9.0, (44.0 * 0.41) / 9.0;
  inquiryHalflengths << 0.41/9.0, 0.41/9.0;
  pos = geoInterface.inquireVoxelPosition ( inquiryCenter.data(), inquiryHalflengths.data(),
    false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  geoInterface.exportMesh ("SolverInterfaceTestGeometry-testVoxelQueryDFGChannel");
}

void SolverInterfaceTestGeometry:: testVoxelQueryFSIChannel()
{
  TRACE();
  SolverInterface geoInterface ( "TestAccessor", 0, 1 );
  configureSolverInterface (
      _pathToTests + "fsichannel.xml", geoInterface );
  validateEquals ( geoInterface.getDimensions(), 2 );
  geoInterface.initialize();
  geoInterface.exportMesh ("SolverInterfaceTestGeometry-testVoxelQueryFSIChannel");
  Eigen::Vector2d inquiryCenter (0.230027434842249633995, 1.87790123456790114531);
  Eigen::Vector2d inquiryHalflengths = Eigen::Vector2d::Constant(0.000562414);
  std::set<int> ids;
  VoxelPosition pos = geoInterface.inquireVoxelPosition (
    inquiryCenter.data(), inquiryHalflengths.data(), false, ids );
  validateEquals ( pos.position(), constants::positionOnGeometry() );
}

void SolverInterfaceTestGeometry:: testVoxelQueryChannelFour ()
{
  TRACE();
  SolverInterface geoInterface ( "TestAccessor", 0, 1 );
  configureSolverInterface ( _pathToTests + "four.xml", geoInterface );
  validateEquals ( geoInterface.getDimensions(), 2 );
  geoInterface.initialize();
  
  Eigen::Vector2d inquiryCenter (-4.44089e-16, -4.44089e-16);
  std::set<int> ids;
  ClosestMesh closest = geoInterface.inquireClosestMesh ( inquiryCenter.data(), ids );
  validateEquals ( closest.position(), constants::positionOnGeometry() );
}

void SolverInterfaceTestGeometry:: testVoxelQueryEpsBox()
{
  TRACE();
  SolverInterface geoInterface ( "TestAccessor", 0, 1 );
  configureSolverInterface ( _pathToTests + "eps-box.xml", geoInterface );
  validateEquals ( geoInterface.getDimensions(), 2 );
  geoInterface.initialize();
  Eigen::Vector2d center(1.1, 0.9);
  Eigen::Vector2d h = Eigen::Vector2d::Constant(0.1);
  std::set<int> ids;
  VoxelPosition pos = geoInterface.inquireVoxelPosition (
    center.data(), h.data(), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[0] = 1.1 - (math::NUMERICAL_ZERO_DIFFERENCE / 2.0);
  pos = geoInterface.inquireVoxelPosition ( center.data(), h.data(), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[1] = 1.0;
  pos = geoInterface.inquireVoxelPosition ( center.data(), h.data(), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[1] = 0.9;
  pos = geoInterface.inquireVoxelPosition ( center.data(), h.data(), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[1] = 0.5;
  pos = geoInterface.inquireVoxelPosition ( center.data(), h.data(), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[1] = 0.4;
  pos = geoInterface.inquireVoxelPosition ( center.data(), h.data(), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[1] = 0.0;
  pos = geoInterface.inquireVoxelPosition ( center.data(), h.data(), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[0] = 1.1 + (math::NUMERICAL_ZERO_DIFFERENCE / 2.0);
  pos = geoInterface.inquireVoxelPosition ( center.data(), h.data(), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[1] = 1.0;
  pos = geoInterface.inquireVoxelPosition ( center.data(), h.data(), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[1] = 0.9;
  pos = geoInterface.inquireVoxelPosition ( center.data(), h.data(), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[1] = 0.5;
  pos = geoInterface.inquireVoxelPosition ( center.data(), h.data(), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[1] = 0.4;
  pos = geoInterface.inquireVoxelPosition ( center.data(), h.data(), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[1] = 0.0;
  pos = geoInterface.inquireVoxelPosition ( center.data(), h.data(), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[0] = 1.1 - (5.0 * math::NUMERICAL_ZERO_DIFFERENCE);
  pos = geoInterface.inquireVoxelPosition ( center.data(), h.data(), false, ids );
  validateEquals ( pos.position(), constants::positionOnGeometry() );

  center[1] = 1.0;
  pos = geoInterface.inquireVoxelPosition ( center.data(), h.data(), false, ids );
  validateEquals ( pos.position(), constants::positionOnGeometry() );

  center[1] = 0.9;
  pos = geoInterface.inquireVoxelPosition ( center.data(), h.data(), false, ids );
  validateEquals ( pos.position(), constants::positionOnGeometry() );

  center[1] = 0.5;
  pos = geoInterface.inquireVoxelPosition ( center.data(), h.data(), false, ids );
  validateEquals ( pos.position(), constants::positionOnGeometry() );

  center[1] = 0.4;
  pos = geoInterface.inquireVoxelPosition ( center.data(), h.data(), false, ids );
  validateEquals ( pos.position(), constants::positionOnGeometry() );

  center[1] = 0.0;
  pos = geoInterface.inquireVoxelPosition ( center.data(), h.data(), false, ids );
  validateEquals ( pos.position(), constants::positionOnGeometry() );
}

void SolverInterfaceTestGeometry:: testConservativeStationaryDataMapping()
{

  TRACE();
  SolverInterface precice("Accessor", 0, 1);
  configureSolverInterface(_pathToTests + "stationary-mapping.xml", precice);
  validateEquals(precice.getDimensions(), 2);

  int meshID = precice.getMeshID("SolverMesh");
  int dataID = precice.getDataID("Forces", meshID);
  int indices[4];
  double values[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
  double x[2];
  x[0] = 0.0; x[1] = 0.0;
  indices[0] = precice.setMeshVertex(meshID, x);
  x[0] = 1.0; x[1] = 0.0;
  indices[1] = precice.setMeshVertex(meshID, x);
  x[0] = 1.0; x[1] = 1.0;
  indices[2] = precice.setMeshVertex(meshID, x);
  x[0] = 0.0; x[1] = 1.0;
  indices[3] = precice.setMeshVertex(meshID, x);
  precice.writeBlockVectorData(dataID, 4, indices, values);

  precice.initialize();
  DEBUG ( "preCICE initialized");
  precice.mapWriteDataFrom(meshID);
  // Validate results
  impl::PtrParticipant p = precice._impl->_accessor;
  DEBUG ( "Participant found");
  validate(p != nullptr);
  DEBUG ( "dataContexts: " << p->_dataContexts << " and dataID: " << dataID);
  validate(p->_dataContexts[dataID] != nullptr);
  mesh::PtrData data = p->_dataContexts[dataID]->toData;
  DEBUG ( "ToData found");
  validate(data.get() != nullptr);
  auto& writtenValues = data->values();

  validateEquals(writtenValues.size(), 8);
  validateNumericalEquals(writtenValues[0], 1.0);
  validateNumericalEquals(writtenValues[1], 2.0);
  validateNumericalEquals(writtenValues[2], 3.0);
  validateNumericalEquals(writtenValues[3], 4.0);
  validateNumericalEquals(writtenValues[4], 7.0); //order is different than above
  validateNumericalEquals(writtenValues[5], 8.0);
  validateNumericalEquals(writtenValues[6], 5.0);
  validateNumericalEquals(writtenValues[7], 6.0);
}


void SolverInterfaceTestGeometry:: testMappingRBF()
{
  TRACE();
  SolverInterface interface("TestAccessor", 0, 1);
  configureSolverInterface(_pathToTests + "mapping-rbf.xml", interface);
  validateEquals(interface.getDimensions(), 2);
  interface.initialize();
  using Eigen::Vector2d;
  
  // Set write data for consistent thin plate splines
  {
    int solverMeshID = interface.getMeshID("SolverMesh-ConsistentTPS");
    // set positions
    Vector2d position ( 0.0, 0.0 );
    int i0 = interface.setMeshVertex ( solverMeshID, position.data() );
    position << 1.0, 0.0;
    int i1 = interface.setMeshVertex ( solverMeshID, position.data() );
    position << 1.0, 1.0;
    int i2 = interface.setMeshVertex ( solverMeshID, position.data() );
    position << 0.0, 1.0;
    int i3 = interface.setMeshVertex ( solverMeshID, position.data() );
    int dataID = interface.getDataID ( "ConsistentTPS", solverMeshID );
    double data = 0.0;
    interface.writeScalarData ( dataID, i0, data );
    data = 2.0;
    interface.writeScalarData ( dataID, i1, data );
    data = 6.0;
    interface.writeScalarData ( dataID, i2, data );
    data = 2.0;
    interface.writeScalarData ( dataID, i3, data );
  }

  // Set  write data for conservative thin plate splines
  {
    int solverMeshID = interface.getMeshID("SolverMesh-ConservativeTPS");
    // set positions
    Vector2d position ( 0.0, 0.0 );
    int i0 = interface.setMeshVertex ( solverMeshID, position.data() );
    position << 1.0, 0.0;
    int i1 = interface.setMeshVertex ( solverMeshID, position.data() );
    position << 1.0, 1.0;
    int i2 = interface.setMeshVertex ( solverMeshID, position.data() );
    position << 0.0, 1.0;
    int i3 = interface.setMeshVertex ( solverMeshID, position.data() );
    int dataID = interface.getDataID("ConservativeTPS", solverMeshID);
    double data = 0.0;
    interface.writeScalarData(dataID, i0, data);
    data = 2.0;
    interface.writeScalarData(dataID, i1, data);
    data = 6.0;
    interface.writeScalarData(dataID, i2, data);
    data = 2.0;
    interface.writeScalarData(dataID, i3, data);
  }

  // Set write data for consistent volume splines
  {
    int solverMeshID = interface.getMeshID("SolverMesh-ConsistentVS");
    // set positions
    Vector2d position ( 0.0, 0.0 );
    int i0 = interface.setMeshVertex ( solverMeshID, position.data() );
    position << 1.0, 0.0;
    int i1 = interface.setMeshVertex ( solverMeshID, position.data() );
    position << 1.0, 1.0;
    int i2 = interface.setMeshVertex ( solverMeshID, position.data() );
    position << 0.0, 1.0;
    int i3 = interface.setMeshVertex ( solverMeshID, position.data() );
    int dataID = interface.getDataID ( "ConsistentVS", solverMeshID );
    Vector2d data ( 0.0, 0.0 );
    interface.writeVectorData ( dataID, i0, data.data() );
    data << 2.0, 2.0;
    interface.writeVectorData ( dataID, i1, data.data() );
    data << 6.0, 6.0;
    interface.writeVectorData ( dataID, i2, data.data() );
    data << 2.0, 2.0;
    interface.writeVectorData ( dataID, i3, data.data() );
  }

  // Set write data for conservative volume splines
  {
    int solverMeshID = interface.getMeshID("SolverMesh-ConservativeVS");
    // set positions
    Vector2d position ( 0.0, 0.0 );
    int i0 = interface.setMeshVertex ( solverMeshID, position.data() );
    position << 1.0, 0.0;
    int i1 = interface.setMeshVertex ( solverMeshID, position.data() );
    position << 1.0, 1.0;
    int i2 = interface.setMeshVertex ( solverMeshID, position.data() );
    position << 0.0, 1.0;
    int i3 = interface.setMeshVertex ( solverMeshID, position.data() );
    int dataID = interface.getDataID ( "ConservativeVS", solverMeshID );
    Vector2d data ( 0.0, 0.0 );
    interface.writeVectorData ( dataID, i0, data.data() );
    data << 2.0, 2.0;
    interface.writeVectorData ( dataID, i1, data.data() );
    data << 6.0, 6.0;
    interface.writeVectorData ( dataID, i2, data.data() );
    data << 2.0, 2.0;
    interface.writeVectorData ( dataID, i3, data.data() );
  }

  interface.advance(1.0);
  interface.finalize();
}

//void SolverInterfaceTestGeometry:: testConservativeDataMapping ()
//{
//# if defined(Dim2)
//  TRACE();
//
//  SolverInterface geoInterface ( "TestAccessor" );
//  geoInterface.configure ( _pathToTests + "SolverInterfaceTestGeometry-conservative-data-mapping-config.xml" );
//  geoInterface.initialize ();
//  int dataID = geoInterface.getDataID ( "Forces" );
//  int meshID = geoInterface.getMeshID ( "itest-cuboid" );
//
//  validateEquals ( geoInterface._participants.size(), 1,
//                   "testConservativeDataMapping()" );
//  impl::PtrParticipant participant = geoInterface._participants[0];
//  validateEquals ( participant->getName(), "TestAccessor",
//                   "testConservativeDataMapping" );
//  validateEquals ( participant->_meshContexts.size(), 1,
//                   "testConservativeDataMapping" );
//  validate ( participant->_meshContexts[0] != NULL,
//             "testConservativeDataMapping" );
//  validate ( participant->_meshContexts[0]->writeMappingContext.mapping.use_count() > 0,
//             "testConservativeDataMapping" );
//  validate ( participant->_meshContexts[0]->readMappingContext.mapping.use_count() > 0,
//             "testConservativeDataMapping" );
//
//  for ( int j=-10; j < 10; j+=2 ) {
//    for ( int i=-10; i < 10; i+=2 ) {
//      // From -10 to +10 in every dimension
//      Vector currentPoint ((double)i, (double)j);
//      Vector data (1.0);
//      geoInterface.setMeshVertex ( meshID, currentPoint );
//      geoInterface.writeData ( dataID, data );
//    }
//  }
//
//  geoInterface.advance ( 0.0 );
//  geoInterface.exportMesh ( "testConservativeDataMapping" );
//
//  Vector integral ( 0.0 );
//  Vector temp;
//  validateEquals ( geoInterface._participants.size(), 1,
//                   "testConservativeDataMapping" );
////  impl::PtrParticipant participant = geoInterface._server->_participants[0];
//  validate ( (int)participant->_dataContexts.size() > dataID,
//             "testConservativeDataMapping" );
//  mesh::Data::Values & values = participant->_dataContexts[dataID]->data->values();
//  for ( int index=0; index < values.size() / 2; index++ ) {
////    temp = tarch::la::dview(values, index);
//    la::readSlice ( temp, values, index*utils::Def::DIM );
//    integral += temp;
//  }
//
//  validate ( tarch::la::equals(integral, Vector(100.0)),
//             "testConservativeDataMapping" );
//  integral = 0.0;
//  geoInterface.integrateData ( dataID, integral );
//  validate ( tarch::la::equals(integral, Vector(100.0)),
//             "testConservativeDataMapping" );
//  geoInterface.advance ( 0.0 );
//  geoInterface.integrateData ( dataID, integral );
//  validate ( tarch::la::equals(integral, Vector(0.0)),
//             "testConservativeDataMapping" );
//# endif // Dim2
//}

void SolverInterfaceTestGeometry:: testPinelli()
{
  TRACE();
  SolverInterface interface("EOF", 0, 1);

  {
    std::string filename = "computeForce.py";
    std::ifstream  src((_pathToTests + filename).c_str(), std::ifstream::in);
    std::ofstream  dst(filename.c_str(), std::ifstream::out);
    dst << src.rdbuf();
  }
  configureSolverInterface(_pathToTests + "testPinelli.xml", interface);
  validateEquals(interface.getDimensions(), 2);
  int meshIdEOF = interface.getMeshID("EOFMesh");
  int dataIdEOFVeloc = interface.getDataID("Velocities",meshIdEOF);
  int dataIdEOFForces = interface.getDataID("Forces",meshIdEOF);

  std::vector<int> vertexIdsEOF;

  Eigen::Vector2d position ( 0.15, 0.15 );
  int vertexId = interface.setMeshVertex ( meshIdEOF, position.data() );
  vertexIdsEOF.push_back(vertexId);
  
  position << 0.15, 0.2;
  vertexId = interface.setMeshVertex ( meshIdEOF, position.data() );
  vertexIdsEOF.push_back(vertexId);

  position << 0.15, 0.25;
  vertexId = interface.setMeshVertex ( meshIdEOF, position.data() );
  vertexIdsEOF.push_back(vertexId);

  position << 0.2, 0.15;
  vertexId = interface.setMeshVertex ( meshIdEOF, position.data() );
  vertexIdsEOF.push_back(vertexId);

  position << 0.2, 0.25;
  vertexId = interface.setMeshVertex ( meshIdEOF, position.data() );
  vertexIdsEOF.push_back(vertexId);

  position << 0.25, 0.15;
  vertexId = interface.setMeshVertex ( meshIdEOF, position.data() );
  vertexIdsEOF.push_back(vertexId);

  position << 0.25, 0.2;
  vertexId = interface.setMeshVertex ( meshIdEOF, position.data() );
  vertexIdsEOF.push_back(vertexId);
  
  position << 0.25, 0.25;
  vertexId = interface.setMeshVertex ( meshIdEOF, position.data() );
  vertexIdsEOF.push_back(vertexId);

  interface.initialize();

  for (int vertexID : vertexIdsEOF){
    double data[2] = {2.0,2.0};
    interface.writeVectorData(dataIdEOFVeloc,vertexID,data);
  }

  interface.advance(0.1);

  Eigen::Vector2d totalForce ( 0.0, 0.0 );

  for (int vertexID : vertexIdsEOF){
    double data[2] = {0.0,0.0};
    interface.readVectorData(dataIdEOFForces,vertexID,data);
    totalForce[0] += data[0];
    totalForce[1] += data[1];
  }

  validateNumericalEquals(totalForce[0]+totalForce[1], 100.0);
}


void SolverInterfaceTestGeometry:: testCustomGeometryCreation()
{
  TRACE();
  using math::equals;
  { // 2D
    using Eigen::Vector2d;
    SolverInterface geo ( "TestAccessor", 0, 1 );
    configureSolverInterface (
        _pathToTests + "solvermesh-2D.xml", geo );
    validateEquals ( geo.getDimensions(), 2 );
    std::string meshName = "custom-geometry";
    int meshID = geo.getMeshID ( meshName );
    geo._impl->_accessor->meshContext(meshID).meshRequirement =
        mapping::Mapping::FULL;
    Vector2d coords0(0.0, 0.0);
    Vector2d coords1(1.0, 0.0);
    Vector2d coords2(1.0, 1.0);
    Vector2d coords3(0.0, 1.0);
    int v0 = geo.setMeshVertex ( meshID, coords0.data() );
    int v1 = geo.setMeshVertex ( meshID, coords1.data() );
    int v2 = geo.setMeshVertex ( meshID, coords2.data() );
    int v3 = geo.setMeshVertex ( meshID, coords3.data() );
    geo.setMeshEdge ( meshID, v0, v1 );
    geo.setMeshEdge ( meshID, v1, v2 );
    geo.setMeshEdge ( meshID, v2, v3 );
    geo.setMeshEdge ( meshID, v3, v0 );

    int size = geo.getMeshVertexSize(meshID);
    validateEquals(size, 4);

    MeshHandle handle = geo.getMeshHandle(meshName);
    VertexHandle vertices = handle.vertices();
    EdgeHandle edges = handle.edges();
    validateEquals (vertices.size(), 4);
    validateEquals (edges.size(), 4);

    VertexIterator vertexIter = vertices.begin();
    validateEquals (vertexIter.vertexID(), 0);
    Vector2d coords;
    coords = Eigen::Map<const Eigen::Vector2d>(vertexIter.vertexCoords());
    validateWithMessage (equals(coords, coords0), coords);
    vertexIter++;
    validateEquals (vertexIter.vertexID(), 1);
    coords = Eigen::Map<const Eigen::Vector2d>(vertexIter.vertexCoords());
    validateWithMessage (equals(coords, coords1), coords);
    vertexIter++;
    validateEquals (vertexIter.vertexID(), 2);
    coords = Eigen::Map<const Eigen::Vector2d>(vertexIter.vertexCoords());
    validateWithMessage (equals(coords, coords2), coords);
    vertexIter++;
    validateEquals (vertexIter.vertexID(), 3);
    coords = Eigen::Map<const Eigen::Vector2d>(vertexIter.vertexCoords());
    validateWithMessage (equals(coords, coords3), coords);
    vertexIter++;
    validate ( not (vertexIter != vertices.end()) );

    EdgeIterator edgeIter = edges.begin();
    validateEquals (edgeIter.vertexID(0), 0);
    validateEquals (edgeIter.vertexID(1), 1);
    edgeIter++;
    validateEquals (edgeIter.vertexID(0), 1);
    validateEquals (edgeIter.vertexID(1), 2);
    edgeIter++;
    validateEquals (edgeIter.vertexID(0), 2);
    validateEquals (edgeIter.vertexID(1), 3);
    edgeIter++;
    validateEquals (edgeIter.vertexID(0), 3);
    validateEquals (edgeIter.vertexID(1), 0);
    edgeIter++;
    validate ( not (edgeIter != edges.end()) );

    geo.exportMesh ("testCustomGeometryCreation-2D");
  }

  { // 3D, Test with manual creation of edges
    using Eigen::Vector3d;
    SolverInterface geo ( "TestAccessor", 0, 1 );
    configureSolverInterface (
        _pathToTests + "solvermesh-3D.xml", geo );
    validateEquals ( geo.getDimensions(), 3 );
    std::string meshName = "custom-geometry";
    int meshID = geo.getMeshID ( meshName );
    geo._impl->_accessor->meshContext(meshID).meshRequirement = mapping::Mapping::FULL;
    Vector3d coords0(0.0,  0.0, 0.0);
    Vector3d coords1(0.5, -0.5, 0.5);
    Vector3d coords2(1.0,  0.0, 1.0);
    Vector3d coords3(0.5,  0.5, 0.5);
    int v0 = geo.setMeshVertex ( meshID, coords0.data() );
    int v1 = geo.setMeshVertex ( meshID, coords1.data() );
    int v2 = geo.setMeshVertex ( meshID, coords2.data() );
    int v3 = geo.setMeshVertex ( meshID, coords3.data() );
    int e0 = geo.setMeshEdge ( meshID, v0, v1 );
    int e1 = geo.setMeshEdge ( meshID, v1, v2 );
    int e2 = geo.setMeshEdge ( meshID, v2, v3 );
    int e3 = geo.setMeshEdge ( meshID, v3, v0 );
    int e4 = geo.setMeshEdge ( meshID, v3, v1 );
    geo.setMeshTriangle ( meshID, e0, e4, e3 );
    geo.setMeshTriangle ( meshID, e1, e2, e4 );

    MeshHandle handle = geo.getMeshHandle(meshName);
    VertexHandle vertices = handle.vertices();
    EdgeHandle edges = handle.edges();
    TriangleHandle triangles = handle.triangles();
    validateEquals (vertices.size(), 4);
    validateEquals (edges.size(), 5);
    validateEquals (triangles.size(), 2);

    VertexIterator vertexIter = vertices.begin();
    validateEquals (vertexIter.vertexID(), 0);
    Vector3d coords;
    coords = Eigen::Map<const Eigen::Vector3d>(vertexIter.vertexCoords());
    validateWithMessage (equals(coords, coords0), coords);
    vertexIter++;
    validateEquals (vertexIter.vertexID(), 1);
    coords = Eigen::Map<const Eigen::Vector3d>(vertexIter.vertexCoords());
    validateWithMessage (equals(coords, coords1), coords);
    vertexIter++;
    validateEquals (vertexIter.vertexID(), 2);
    coords = Eigen::Map<const Eigen::Vector3d>(vertexIter.vertexCoords());
    validateWithMessage (equals(coords, coords2), coords);
    vertexIter++;
    validateEquals (vertexIter.vertexID(), 3);
    coords = Eigen::Map<const Eigen::Vector3d>(vertexIter.vertexCoords());
    validateWithMessage (equals(coords, coords3), coords);
    vertexIter++;
    validate (not (vertexIter != vertices.end()));

    EdgeIterator edgeIter = edges.begin();
    validateEquals (edgeIter.vertexID(0), 0);
    validateEquals (edgeIter.vertexID(1), 1);
    edgeIter++;
    validateEquals (edgeIter.vertexID(0), 1);
    validateEquals (edgeIter.vertexID(1), 2);
    edgeIter++;
    validateEquals (edgeIter.vertexID(0), 2);
    validateEquals (edgeIter.vertexID(1), 3);
    edgeIter++;
    validateEquals (edgeIter.vertexID(0), 3);
    validateEquals (edgeIter.vertexID(1), 0);
    edgeIter++;
    validateEquals (edgeIter.vertexID(0), 3);
    validateEquals (edgeIter.vertexID(1), 1);
    edgeIter++;
    validate (not (edgeIter != edges.end()));

    TriangleIterator triangleIter = triangles.begin();
    validateEquals (triangleIter.vertexID(0), 0);
    validateEquals (triangleIter.vertexID(1), 1);
    validateEquals (triangleIter.vertexID(2), 3);
    triangleIter++;
    validateEquals (triangleIter.vertexID(0), 1);
    validateEquals (triangleIter.vertexID(1), 2);
    validateEquals (triangleIter.vertexID(2), 3);
    triangleIter++;
    validate (not (triangleIter != triangles.end()));

    geo.exportMesh ( "testCustomGeometryCreation-3D-manual" );
  }
  { // 3D, Test with automatic creation of edges
    using Eigen::Vector3d;
    SolverInterface geo ( "TestAccessor", 0, 1 );
    configureSolverInterface ( _pathToTests + "solvermesh-3D.xml", geo );
    validateEquals ( geo.getDimensions(), 3 );
    std::string meshName = "custom-geometry";
    int meshID = geo.getMeshID ( meshName );
    geo._impl->_accessor->meshContext(meshID).meshRequirement = mapping::Mapping::FULL;
    Vector3d coords0(0.0,  0.0, 0.0);
    Vector3d coords1(0.5, -0.5, 0.5);
    Vector3d coords2(1.0,  0.0, 1.0);
    Vector3d coords3(0.5,  0.5, 0.5);
    int v0 = geo.setMeshVertex ( meshID, coords0.data() );
    int v1 = geo.setMeshVertex ( meshID, coords1.data() );
    int v2 = geo.setMeshVertex ( meshID, coords2.data() );
    int v3 = geo.setMeshVertex ( meshID, coords3.data() );
    geo.setMeshTriangleWithEdges ( meshID, v0, v1, v3 );
    geo.setMeshTriangleWithEdges ( meshID, v1, v2, v3 );

    MeshHandle handle = geo.getMeshHandle(meshName);
    VertexHandle vertices = handle.vertices();
    EdgeHandle edges = handle.edges();
    TriangleHandle triangles = handle.triangles();
    validateEquals (vertices.size(), 4);
    validateEquals (edges.size(), 5);
    validateEquals (triangles.size(), 2);

    VertexIterator vertexIter = vertices.begin();
    validateEquals (vertexIter.vertexID(), 0);
    Vector3d coords;
    coords = Eigen::Map<const Eigen::Vector3d>(vertexIter.vertexCoords());
    validateWithMessage (equals(coords, coords0), coords);
    vertexIter++;
    validateEquals (vertexIter.vertexID(), 1);
    coords = Eigen::Map<const Eigen::Vector3d>(vertexIter.vertexCoords());
    validateWithMessage (equals(coords, coords1), coords);
    vertexIter++;
    validateEquals (vertexIter.vertexID(), 2);
    coords = Eigen::Map<const Eigen::Vector3d>(vertexIter.vertexCoords());
    validateWithMessage (equals(coords, coords2), coords);
    vertexIter++;
    validateEquals (vertexIter.vertexID(), 3);
    coords = Eigen::Map<const Eigen::Vector3d>(vertexIter.vertexCoords());
    validateWithMessage (equals(coords, coords3), coords);
    vertexIter++;
    validate (not (vertexIter != vertices.end()));

    EdgeIterator edgeIter = edges.begin();
    validateEquals (edgeIter.vertexID(0), 0);
    validateEquals (edgeIter.vertexID(1), 1);
    edgeIter++;
    validateEquals (edgeIter.vertexID(0), 1);
    validateEquals (edgeIter.vertexID(1), 3);
    edgeIter++;
    validateEquals (edgeIter.vertexID(0), 3);
    validateEquals (edgeIter.vertexID(1), 0);
    edgeIter++;
    validateEquals (edgeIter.vertexID(0), 1);
    validateEquals (edgeIter.vertexID(1), 2);
    edgeIter++;
    validateEquals (edgeIter.vertexID(0), 2);
    validateEquals (edgeIter.vertexID(1), 3);
    edgeIter++;
    validate (not (edgeIter != edges.end()));

    TriangleIterator triangleIter = triangles.begin();
    validateEquals (triangleIter.vertexID(0), 0);
    validateEquals (triangleIter.vertexID(1), 1);
    validateEquals (triangleIter.vertexID(2), 3);
    triangleIter++;
    validateEquals (triangleIter.vertexID(0), 1);
    validateEquals (triangleIter.vertexID(1), 2);
    validateEquals (triangleIter.vertexID(2), 3);
    triangleIter++;
    validate (not (triangleIter != triangles.end()));

    geo.exportMesh ( "testCustomGeometryCreation-3D-auto" );
  }
}

void SolverInterfaceTestGeometry:: testBug()
{
  TRACE();

  {
    std::string filename = "testBug-geometry.wrl";
    std::ifstream  src((_pathToTests + filename).c_str(), std::ifstream::in);
    std::ofstream  dst(filename.c_str(), std::ifstream::out);
    dst << src.rdbuf();
  }

  SolverInterface interface ( "TestAccessor", 0, 1 );
  configureSolverInterface (
      _pathToTests + "testBug.xml", interface );
  validateEquals ( interface.getDimensions(), 3 );
  std::set<int> meshIDs = interface.getMeshIDs ();
  interface.initialize ();
  Eigen::Vector3d center ( 0.0, 0.0, 0.0 );
  Eigen::Vector3d h = Eigen::Vector3d::Constant ( 0.0166666666666666 );
  VoxelPosition pos;
  pos = interface.inquireVoxelPosition ( center.data(), h.data(), false, meshIDs );
  validateEquals ( pos.position(), constants::positionOnGeometry() );
  center.setConstant(-0.01666666666666666);
  int pointPos = interface.inquirePosition ( center.data(), meshIDs );
  validateEquals ( pointPos, constants::positionInsideOfGeometry() );
}

void SolverInterfaceTestGeometry:: testBug2()
{
  TRACE();
  SolverInterface interface("Participant-testBug2", 0, 1);
  configureSolverInterface(_pathToTests + "testBug2.xml", interface);
  validateEquals ( interface.getDimensions(), 3 );
  std::set<int> meshIDs = interface.getMeshIDs();
  interface.initialize();
  Eigen::Vector3d center(7.6875, 0.7625, 3.7125);
  Eigen::Vector3d h = Eigen::Vector3d::Constant(0.0125);
  VoxelPosition pos;
  pos = interface.inquireVoxelPosition ( center.data(), h.data(), false, meshIDs );
  validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
  center << 7.675, 0.75, 3.7;
  h.setConstant(0.025);
  pos = interface.inquireVoxelPosition ( center.data(), h.data(), false, meshIDs );
  validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
  int pointPos = interface.inquirePosition ( center.data(), meshIDs );
  validateEquals ( pointPos, constants::positionInsideOfGeometry() );
}

void SolverInterfaceTestGeometry:: testBug3()
{
  TRACE();
  SolverInterface interface("Peano", 0, 1);
  configureSolverInterface(_pathToTests + "testBug3.xml", interface);
  validateEquals(interface.getDimensions(), 3);
  std::set<int> meshIDs = interface.getMeshIDs();
  interface.initialize();
  
  // This query was buggy
  Eigen::Vector3d center(0.22549, 0.5, 0.205882);
  Eigen::Vector3d h(0.00980392, 0.00980392, 0.00980392);
  VoxelPosition pos;
  pos = interface.inquireVoxelPosition(center.data(), h.data(), false, meshIDs);
  validateEquals(pos.position(), constants::positionInsideOfGeometry());
//  interface.exportMesh("testBug3");
}

void SolverInterfaceTestGeometry:: testBug4()
{
  TRACE();
  SolverInterface interface("Peano", 0, 1);
  configureSolverInterface(_pathToTests + "testBug4.xml", interface);
  validateEquals(interface.getDimensions(), 3);
  std::set<int> meshIDs = interface.getMeshIDs();
  interface.initialize();
  interface.exportMesh("testBug4-init");
  
  // Completely inside in Peano fluid domain as answered by preCICE
  Eigen::Vector3d center(0.275, 0.525, 0.275);
  Eigen::Vector3d h(0.025, 0.025, 0.025);
  VoxelPosition pos;
  pos = interface.inquireVoxelPosition(center.data(), h.data(), false, meshIDs);
  validateEquals(pos.position(), constants::positionInsideOfGeometry());
  //interface.exportMesh("testBug4");
}

void SolverInterfaceTestGeometry:: testBug5()
{
  TRACE();

  {
    std::string filename = "testBug5-geometry.wrl";
    std::ifstream  src((_pathToTests + filename).c_str(), std::ifstream::in);
    std::ofstream  dst(filename.c_str(), std::ifstream::out);
    dst << src.rdbuf();
  }

  SolverInterface interface("Peano", 0, 1);
  configureSolverInterface(_pathToTests + "testBug5.xml", interface);
  validateEquals(interface.getDimensions(), 3);
  std::set<int> meshIDs = interface.getMeshIDs();
  interface.initialize();
  
  // On geometry in Peano fluid domain as answered wrongly by preCICE
  Eigen::Vector3d point(4.7777777777777777, 4.06666666666666666, 4.0666666666666666);
  int pos = interface.inquirePosition(point.data(), meshIDs);
  validateEquals(pos, constants::positionOutsideOfGeometry());

  Eigen::Vector3d h = Eigen::Vector3d::Constant(0.3555555555555555);
  VoxelPosition voxelPos = interface.inquireVoxelPosition(point.data(), h.data(), false, meshIDs);

  h.setConstant(3.950617283950000e-2);
  point << 4.343209876543000, 4.0666666666666666, 4.106172839509999;
//  INFO("----------------------------- START");
  voxelPos = interface.inquireVoxelPosition(point.data(), h.data(), false, meshIDs);
//  INFO("----------------------------- END, pos = " << voxelPos.position()
//               << ", ids.size = " << voxelPos.meshIDs().size());
  //mesh::Mesh found("Found", 3, false);
//  EdgeIterator it = voxelPos.contentHandle().edges().begin();
//  for (;it != voxelPos.contentHandle().edges().end(); it++){
//    INFO("Edge from " << wrap<3>(it.vertexCoords(0)) << " to " <<
//                 wrap<3>(it.vertexCoords(1)));
//    mesh::Vertex& v0 = found.createVertex(wrap<3>(it.vertexCoords(0)));
//    mesh::Vertex& v1 = found.createVertex(wrap<3>(it.vertexCoords(1)));
//    found.createEdge(v0, v1);
//  }
  //for (mesh::Edge& edge : voxelPos.contentHandle().triangles()){
  //  mesh::Vertex& v0 = found.createVertex(edge.vertex(0).getCoords());
  //  mesh::Vertex& v1 = found.createVertex(edge.vertex(1).getCoords());
  //  found.createEdge(v0, v1);
  //}
//  mesh::Mesh voxelMesh("VoxelMesh", 3, false);
//  geometry::Cuboid cuboid(point - h, 1.0, 2.0 * h);
//  cuboid.create(voxelMesh);
//  io::ExportVTK exportVTK(false);
//  exportVTK.doExport("testBug5-voxel", voxelMesh);
//  assertion(false);
  validateEquals(voxelPos.position(), constants::positionOutsideOfGeometry());
  validateEquals(voxelPos.meshIDs().size(), 0);
  interface.exportMesh("testBug5");
}

//void SolverInterfaceTestGeometry:: testBug6()
//{
//  TRACE();
//  SolverInterface interface("Peano", 0, 1);
//  configureSolverInterface(_pathToTests + "testBug6.xml", interface);
//  validateEquals(interface.getDimensions(), 3);
//  std::set<int> meshIDs;
//  interface.initialize();
//  using namespace tarch::la;
//  interface.exportMesh("testBug6");
//
//  // On geometry in Peano fluid domain as answered wrongly by preCICE
//  utils::Vector3D point(4.7777777777777777, 4.06666666666666666, 4.0666666666666666);
//  int pos = interface.inquirePosition(raw(point), meshIDs);
//  validateEquals(pos, constants::positionOutsideOfGeometry());
//
//  utils::Vector3D h(0.3555555555555555);
//  VoxelPosition voxelPos = interface.inquireVoxelPosition(raw(point), h.data(), false, meshIDs);
//  validateEquals(voxelPos.position(), constants::positionOnGeometry());
//  validateEquals(voxelPos.meshIDs().size(), 0);
//}

void SolverInterfaceTestGeometry:: testUpdateSpacetree()
{
  TRACE();
  for ( int dim=2; dim <= 3; dim++ ){
    SolverInterface interface("Accessor", 0, 1);
    std::string configName;
    if (dim == 2){
      configName = _pathToTests + "update-spacetree-2D.xml";
    }
    else {
      configName = _pathToTests + "update-spacetree-3D.xml";
    }
    configureSolverInterface(configName, interface);
    validateEquals(interface.getDimensions(), dim);
    interface.initialize();
    std::ostringstream filename;
    filename << "testUpdateSpacetree-" << dim << "D-";
    interface.exportMesh(filename.str() + "0");

    MeshHandle handle = interface.getMeshHandle("Mesh");
    int meshID = interface.getMeshID ("Mesh");
    int dataID = interface.getDataID(constants::dataDisplacements(), meshID);
    double value[dim];
    value[0] = 0.1;

    for (VertexIterator iter : handle.vertices()) {
      interface.writeVectorData(dataID, iter.vertexID(), value);
      interface.inquirePosition(iter.vertexCoords(), std::set<int>());
    }
    interface.exportMesh(filename.str() + "1");
    interface.advance(1.0);

    for (VertexIterator iter : handle.vertices()) {
      interface.writeVectorData(dataID, iter.vertexID(), value);
      interface.inquirePosition(iter.vertexCoords(), std::set<int>());
    }
    interface.exportMesh(filename.str() + "2");
    interface.advance(1.0);

    for (VertexIterator iter : handle.vertices()) {
      interface.writeVectorData(dataID, iter.vertexID(), value);
      interface.inquirePosition(iter.vertexCoords(), std::set<int>());
    }
    interface.exportMesh(filename.str() + "3");
    interface.advance(1.0);

    for (VertexIterator iter : handle.vertices()) {
      interface.writeVectorData(dataID, iter.vertexID(), value);
      interface.inquirePosition(iter.vertexCoords(), std::set<int>());
    }
    interface.exportMesh(filename.str() + "4");
    interface.advance(1.0);

    interface.finalize();
  }
}

void SolverInterfaceTestGeometry:: testMultipleMeshSpacetree()
{
  TRACE();
  using math::equals;
  
  { // Tests A: first mesh no spacetree, second, third spacetree
    SolverInterface interface("Accessor", 0, 1);
    impl::SolverInterfaceImpl* impl = interface._impl.get();
    std::string configName = _pathToTests + "multiple-mesh-spacetree-a.xml";
    configureSolverInterface(configName, interface);
    validateEquals(interface.getDimensions(), 2);
    interface.initialize();
    std::vector<int> meshIDs;
    meshIDs.push_back(interface.getMeshID("FirstMesh"));
    meshIDs.push_back(interface.getMeshID("SecondMesh"));
    meshIDs.push_back(interface.getMeshID("ThirdMesh"));
    std::set<int> allMeshIDs = interface.getMeshIDs();

    std::set<int> inputIDs;
    std::vector<int> outputIDs(3, -1);
    Eigen::Vector3i expected;
    expected[0] = impl->markedQueryDirectly();
    expected[1] = impl->markedQuerySpacetree();
    expected[2] = impl->markedSkip();
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(Eigen::Vector3i(outputIDs.data()), expected),
                        outputIDs);

    expected[0] = impl->markedQueryDirectly();
    expected[1] = impl->markedQuerySpacetree();
    expected[2] = impl->markedSkip();
    std::fill(outputIDs.begin(), outputIDs.end(), -1);
    impl->selectInquiryMeshIDs(allMeshIDs, outputIDs);
    validateWithMessage(equals(Eigen::Vector3i(outputIDs.data()), expected),
                        outputIDs);

    inputIDs.insert(0);
    expected[0] = impl->markedQueryDirectly();
    expected[1] = impl->markedSkip();
    expected[2] = impl->markedSkip();
    std::fill(outputIDs.begin(), outputIDs.end(), -1);
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(Eigen::Vector3i(outputIDs.data()), expected),
                        outputIDs);

    inputIDs.clear();
    inputIDs.insert(1);
    expected[0] = impl->markedSkip();
    expected[1] = impl->markedQueryDirectly();
    expected[2] = impl->markedSkip();
    std::fill(outputIDs.begin(), outputIDs.end(), -1);
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(Eigen::Vector3i(outputIDs.data()), expected),
                        outputIDs);

    inputIDs.clear();
    inputIDs.insert(2);
    expected[0] = impl->markedSkip();
    expected[1] = impl->markedSkip();
    expected[2] = impl->markedQueryDirectly();
    std::fill(outputIDs.begin(), outputIDs.end(), -1);
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(Eigen::Vector3i(outputIDs.data()), expected),
                        outputIDs);

    inputIDs.clear();
    inputIDs.insert(0);
    inputIDs.insert(2);
    expected[0] = impl->markedQueryDirectly();
    expected[1] = impl->markedSkip();
    expected[2] = impl->markedQueryDirectly();
    std::fill(outputIDs.begin(), outputIDs.end(), -1);
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(Eigen::Vector3i(outputIDs.data()), expected),
                        outputIDs);

    inputIDs.clear();
    inputIDs.insert(1);
    inputIDs.insert(2);
    expected[0] = impl->markedSkip();
    expected[1] = impl->markedQuerySpacetree();
    expected[2] = impl->markedSkip();
    std::fill(outputIDs.begin(), outputIDs.end(), -1);
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(Eigen::Vector3i(outputIDs.data()), expected),
                        outputIDs);

    inputIDs.clear();
    inputIDs.insert(0);
    inputIDs.insert(1);
    expected[0] = impl->markedQueryDirectly();
    expected[1] = impl->markedQueryDirectly();
    expected[2] = impl->markedSkip();
    std::fill(outputIDs.begin(), outputIDs.end(), -1);
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(Eigen::Vector3i(outputIDs.data()), expected),
                        outputIDs);

  }
  { // Tests B: first mesh has spacetree, second, third not
    SolverInterface interface("Accessor", 0, 1);
    impl::SolverInterfaceImpl* impl = interface._impl.get();
    std::string configName = _pathToTests + "multiple-mesh-spacetree-b.xml";
    configureSolverInterface(configName, interface);
    validateEquals(interface.getDimensions(), 2);
    interface.initialize();
    std::vector<int> meshIDs;
    meshIDs.push_back(interface.getMeshID("FirstMesh"));
    meshIDs.push_back(interface.getMeshID("SecondMesh"));
    meshIDs.push_back(interface.getMeshID("ThirdMesh"));
    std::set<int> allMeshIDs = interface.getMeshIDs();

    std::set<int> inputIDs;
    std::vector<int> outputIDs(3, -1);
    Eigen::Vector3i expected;
    expected[0] = impl->markedQuerySpacetree();
    expected[1] = impl->markedQueryDirectly();
    expected[2] = impl->markedQueryDirectly();
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(Eigen::Vector3i(outputIDs.data()), expected),
                        outputIDs);

    expected[0] = impl->markedQuerySpacetree();
    expected[1] = impl->markedQueryDirectly();
    expected[2] = impl->markedQueryDirectly();
    std::fill(outputIDs.begin(), outputIDs.end(), -1);
    impl->selectInquiryMeshIDs(allMeshIDs, outputIDs);
    validateWithMessage(equals(Eigen::Vector3i(outputIDs.data()), expected),
                        outputIDs);

    inputIDs.insert(0);
    expected[0] = impl->markedQuerySpacetree();
    expected[1] = impl->markedSkip();
    expected[2] = impl->markedSkip();
    std::fill(outputIDs.begin(), outputIDs.end(), -1);
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(Eigen::Vector3i(outputIDs.data()), expected),
                        outputIDs);

    inputIDs.clear();
    inputIDs.insert(1);
    expected[0] = impl->markedSkip();
    expected[1] = impl->markedQueryDirectly();
    expected[2] = impl->markedSkip();
    std::fill(outputIDs.begin(), outputIDs.end(), -1);
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(Eigen::Vector3i(outputIDs.data()), expected),
                        outputIDs);

    inputIDs.clear();
    inputIDs.insert(2);
    expected[0] = impl->markedSkip();
    expected[1] = impl->markedSkip();
    expected[2] = impl->markedQueryDirectly();
    std::fill(outputIDs.begin(), outputIDs.end(), -1);
    
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(Eigen::Vector3i(outputIDs.data()), expected),
                        outputIDs);

    inputIDs.clear();
    inputIDs.insert(0);
    inputIDs.insert(2);
    expected[0] = impl->markedQuerySpacetree();
    expected[1] = impl->markedSkip();
    expected[2] = impl->markedQueryDirectly();
    std::fill(outputIDs.begin(), outputIDs.end(), -1);
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(Eigen::Vector3i(outputIDs.data()), expected),
                        outputIDs);

    inputIDs.clear();
    inputIDs.insert(1);
    inputIDs.insert(2);
    expected[0] = impl->markedSkip();
    expected[1] = impl->markedQueryDirectly();
    expected[2] = impl->markedQueryDirectly();
    std::fill(outputIDs.begin(), outputIDs.end(), -1);
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(Eigen::Vector3i(outputIDs.data()), expected),
                        outputIDs);

    inputIDs.clear();
    inputIDs.insert(0);
    inputIDs.insert(1);
    expected[0] = impl->markedQuerySpacetree();
    expected[1] = impl->markedQueryDirectly();
    expected[2] = impl->markedSkip();
    std::fill(outputIDs.begin(), outputIDs.end(), -1);
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(Eigen::Vector3i(outputIDs.data()), expected),
                        outputIDs);
  }
}

}} // namespace precice, tests



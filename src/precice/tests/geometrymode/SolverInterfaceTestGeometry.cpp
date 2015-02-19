// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
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
#include "tarch/la/WrappedVector.h"
#include <vector>
#include <set>
#include <algorithm>
#include <boost/range/algorithm.hpp>

#include "tarch/tests/TestCaseFactory.h"
registerIntegrationTest(precice::tests::SolverInterfaceTestGeometry)

namespace precice {
namespace tests {

using namespace tarch::la;

tarch::logging::Log SolverInterfaceTestGeometry::
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
  preciceTrace("run ()");
  std::vector<int> ranks;
  ranks += 0;
  typedef utils::Parallel Par;
  Par::Communicator comm = Par::getRestrictedCommunicator(ranks);
  PRECICE_MASTER_ONLY {
    preciceTrace("run master ()");
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
    testMethod(testBug);
    testMethod(testBug2);
#   endif // not PRECICE_NO_SPIRIT2
    testMethod(testBug3);
    testMethod(testBug4);
    testMethod(testBug5);
//    testMethod(testBug6); // no bug actually
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
  preciceTrace1 ( "configureSolverInterface()", configFilename );
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
  preciceTrace ( "testConfiguration()" );
  mesh::Mesh::resetGeometryIDsGlobally ();
  preciceDebug ( "Test 2D configuration");
  { // 2D
    SolverInterface geoInterface ( "TestAccessor", 0, 1 );
    configureSolverInterface (
        _pathToTests + "configuration2D.xml", geoInterface );
    validateEquals ( geoInterface.getDimensions(), 2 );
    geoInterface.initialize ();
    geoInterface.exportMesh ( "testConfiguration2D" );
  }
  preciceDebug ( "Test 3D configuration");
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
  preciceTrace ( "testSearchQuery()" );
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
    utils::DynVector pos(dim,0.0);
    assign(pos) = 50.0;
    geoInterface.setMeshVertex(meshID,raw(pos));

    _geoID = geoInterface.getMeshID("itest-cuboid");

    std::set<int> ids;
    assign(pos) = 0.0;
    ClosestMesh closest = geoInterface.inquireClosestMesh (raw(pos), ids);
    validateEquals ( closest.meshIDs().size(), 1 );
    validateEquals ( closest.meshIDs()[0], _geoID );
    validateEquals ( closest.position(), constants::positionOutsideOfGeometry() );

    assign(pos) = 4.0;
    closest = geoInterface.inquireClosestMesh (raw(pos), ids);
    validateEquals ( closest.meshIDs().size(), 1 );
    validateEquals ( closest.meshIDs()[0], _geoID );
    validateEquals ( closest.position(), constants::positionOutsideOfGeometry() );

    assign(pos) = 5.0;
    closest = geoInterface.inquireClosestMesh (raw(pos), ids);
    validateEquals ( closest.meshIDs().size(), 1 );
    //validateEquals ( closest.meshIDs()[0], _geoID );  // why does this not work
    validateEquals ( closest.position(), constants::positionOnGeometry() );

    assign(pos) = 6.0;
    closest = geoInterface.inquireClosestMesh ( raw(pos), ids );
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
   preciceTrace("testVoxelQuery()");

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
     utils::DynVector posVertex(dim,50.0);
     geoInterface.setMeshVertex(meshID,raw(posVertex));


     _geoID = geoInterface.getMeshID ("itest-cuboid");

     // Voxel completely contained in center
     std::set<int> ids;
     bool include = false;
     utils::DynVector center(dim, 0.0);
     utils::DynVector h(dim, 0.1);
     VoxelPosition pos = geoInterface.inquireVoxelPosition (
                         raw(center), raw(h), include, ids );
     validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

     if (dim == 2){
       // Voxels in corners
       assignList(center) = -4.0, -4.0;
       assignList(h) = 1.0, 1.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids);
       validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
       assignList(center) = 4.0, -4.0;
       assignList(h) = 1.0, 1.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
       assignList(center) = -4.0, 4.0;
       assignList(h) = 1.0, 1.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
       assignList(center) = 4.0, 4.0;
       assignList(h) = 1.0, 1.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

       // Voxels on sides (outside of geometry)
       assignList(center) = -4.0, 0.0;
       assignList(h) = 1.0, 1.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
       assignList(center) = 4.0, 0.0;
       assignList(h) = 1.0, 1.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
       assignList(center) = -4.0, 4.0;
       assignList(h) = 1.0, 1.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
       assignList(center) = 4.0, 4.0;
       assignList(h) = 1.0, 1.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
       assignList(center) = -4.0, 0.0;
       assignList(h) = 1.0 + (tarch::la::NUMERICAL_ZERO_DIFFERENCE/2.0), 1.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

       // Voxels on sides (inside of geometry)
       assignList(center) = -6.0, 0.0;
       assignList(h) = 1.0, 1.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
       assignList(center) = 6.0, 0.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
       assignList(center) = -6.0, 6.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
       assignList(center) = 6.0, 6.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionInsideOfGeometry() );

       assignList(center) = 5.0, 5.0;
       assignList(h) = 5.0, 5.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.meshIDs().size(), 1 );
       validateEquals ( pos.meshIDs()[0], _geoID );
       geoInterface.exportMesh ("SolverInterfaceTestGeometry-testVoxelQuery-2D");
     }
     else { // 3D
       // Voxels in corners (outside of geometry)
       assignList(center) = -4.0, -4.0, -4.0;
       assignList(h) = 1.0, 1.0, 1.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
       assignList(center) = 4.0, -4.0, -4.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
       assignList(center) = -4.0, 4.0, -4.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
       assignList(center) = 4.0, 4.0, -4.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
       assignList(center) = -4.0, -4.0, 4.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
       assignList(center) = 4.0, -4.0, 4.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
       assignList(center) = -4.0, 4.0, 4.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
       assignList(center) = 4.0, 4.0, 4.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

       // Voxels in corners (inside of geometry)
       assignList(center) = -6.0, -6.0, -6.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
       assignList(center) = 6.0, -6.0, -6.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
       assignList(center) = -6.0, 6.0, -6.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
       assignList(center) = 6.0, 6.0, -6.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
       assignList(center) = -6.0, -6.0, 6.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
       assignList(center) = 6.0, -6.0, 6.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
       assignList(center) = -6.0, 6.0, 6.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
       assignList(center) = 6.0, 6.0, 6.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionInsideOfGeometry() );

       // Voxels on sides (outside of geometry)
       assignList(center) = -4.0, 0.0, 0.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionOutsideOfGeometry());
       assignList(center) = 4.0, 0.0, 0.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
       assignList(center) = 0.0, -4.0, 0.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
       assignList(center) = 0.0, 4.0, 0.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
       assignList(center) = 0.0, 0.0, -4.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );
       assignList(center) = 0.0, 0.0, 4.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

       // Voxels on sides (inside of geometry)
       assignList(center) = -6.0, 0.0, 0.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
       assignList(center) = 6.0, 0.0, 0.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
       assignList(center) = 0.0, -6.0, 0.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
       assignList(center) = 0.0, 6.0, 0.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
       assignList(center) = 0.0, 0.0, -6.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
       assignList(center) = 0.0, 0.0, 6.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.position(), constants::positionInsideOfGeometry() );

       // Voxels on geometry
       assignList(center) = 5.0, 5.0, 5.0;
       assignList(h) = 5.0, 5.0, 5.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.meshIDs().size(), 1 );
       validateEquals ( pos.meshIDs()[0], _geoID );

       assignList(center) = -5.0, -5.0, -5.0;
       pos = geoInterface.inquireVoxelPosition (raw(center), raw(h), include, ids );
       validateEquals ( pos.meshIDs().size(), 1 );
       validateEquals ( pos.meshIDs()[0], _geoID );
       geoInterface.exportMesh ("SolverInterfaceTestGeometry-testVoxelQuery-3D");
     }
   }
}

void SolverInterfaceTestGeometry:: testDataActions()
{
  preciceTrace("testDataActions()");
  using namespace tarch::la;
  SolverInterface geo("Accessor", 0, 1);
  configureSolverInterface(_pathToTests + "testDataActions.xml", geo);
  impl::SolverInterfaceImpl* impl = geo._impl;
  int meshID = geo.getMeshID("Box");
  int dataID = geo.getDataID("VectorData", meshID);
  geo.initialize();
  mesh::PtrMesh mesh = impl->_accessor->meshContext(meshID).mesh;
  std::vector<utils::Vector3D> coords ( mesh->vertices().size() );
  for ( size_t i=0; i < mesh->vertices().size(); i++ ){
    coords[i] = mesh->vertices()[i].getCoords();
  }
  mesh::PtrData data = mesh->data ( dataID );
  assign(data->values()) = 1.0;

  geo.advance ( 1.0 );

  for ( size_t i=0; i < mesh->vertices().size(); i++ ){
    validate ( equals(coords[i] + 1.0, mesh->vertices()[i].getCoords()) );
  }

  geo.advance ( 1.0 );

  for ( size_t i=0; i < mesh->vertices().size(); i++ ){
    validate ( equals(coords[i] + 2.0, mesh->vertices()[i].getCoords()) );
  }

  geo.finalize ();
}

void SolverInterfaceTestGeometry:: testVoxelQueryMultipleGeometryIDs()
{
  preciceTrace ( "testVoxelQueryMultipleGeometryIDsest()" );
  using utils::Vector2D;
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
    raw(Vector2D(0.0, 0.0)), raw(Vector2D(1.0, 1.0)), false, ids );
  validateEquals ( pos.position(), constants::positionOnGeometry() );
  validateEquals ( pos.meshIDs().size(), 2 );
  if ( pos.meshIDs()[0] == idQuad0Base ) {
    validateEquals ( pos.meshIDs()[1], idQuad0Side2 );
  }
  else {
    validateEquals ( pos.meshIDs()[0], idQuad0Side2 );
    validateEquals ( pos.meshIDs()[1], idQuad0Base );
  }

  pos = geoInterface.inquireVoxelPosition ( raw(Vector2D(5.0, 5.0)),
                                            raw(Vector2D(1.5, 1.5)),
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
  preciceTrace ( "testVoxelQueryDFGChannel()" );
  using utils::Vector2D;
  SolverInterface geoInterface ( "TestAccessor", 0, 1 );
  configureSolverInterface (
      _pathToTests + "dfgchannel.xml", geoInterface );
  validateEquals ( geoInterface.getDimensions(), 2 );
  geoInterface.initialize();
  _geoID = geoInterface.getMeshID ("itest-dfg-cuboid");
  Vector2D inquiryCenter (0.41/9.0, (46.0 * 0.41) / 9.0);
  Vector2D inquiryHalflengths (0.41/9.0, 0.41/9.0);
  std::set<int> ids;
  VoxelPosition pos = geoInterface.inquireVoxelPosition (
    raw(inquiryCenter), raw(inquiryHalflengths), false, ids );
  validateEquals ( pos.position(), constants::positionInsideOfGeometry() );

  assignList(inquiryCenter) = 0.0455556, 2.09556;
  assignList(inquiryHalflengths) = 0.0455556, 0.0455556;
  pos = geoInterface.inquireVoxelPosition ( raw(inquiryCenter), raw(inquiryHalflengths),
    false, ids );
  validateEquals ( pos.position(), constants::positionInsideOfGeometry() );

  assignList(inquiryCenter) = 0.41/9.0, (44.0 * 0.41) / 9.0;
  assignList(inquiryHalflengths) = 0.41/9.0, 0.41/9.0;
  pos = geoInterface.inquireVoxelPosition ( raw(inquiryCenter), raw(inquiryHalflengths),
    false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  geoInterface.exportMesh ("SolverInterfaceTestGeometry-testVoxelQueryDFGChannel");
}

void SolverInterfaceTestGeometry:: testVoxelQueryFSIChannel()
{
  preciceTrace ( "testVoxelQueryFSIChannel()" );
  using utils::Vector2D;
  SolverInterface geoInterface ( "TestAccessor", 0, 1 );
  configureSolverInterface (
      _pathToTests + "fsichannel.xml", geoInterface );
  validateEquals ( geoInterface.getDimensions(), 2 );
  geoInterface.initialize();
  geoInterface.exportMesh ("SolverInterfaceTestGeometry-testVoxelQueryFSIChannel");
  Vector2D inquiryCenter (0.230027434842249633995, 1.87790123456790114531);
  Vector2D inquiryHalflengths (0.000562414);
  std::set<int> ids;
  VoxelPosition pos = geoInterface.inquireVoxelPosition (
                      raw(inquiryCenter), raw(inquiryHalflengths), false, ids );
  validateEquals ( pos.position(), constants::positionOnGeometry() );
}

void SolverInterfaceTestGeometry:: testVoxelQueryChannelFour ()
{
   preciceTrace ( "testVoxelQueryChannelFour()");
   SolverInterface geoInterface ( "TestAccessor", 0, 1 );
   configureSolverInterface (
         _pathToTests + "four.xml", geoInterface );
   validateEquals ( geoInterface.getDimensions(), 2 );
   geoInterface.initialize();

   utils::Vector2D inquiryCenter (-4.44089e-16, -4.44089e-16);
   std::set<int> ids;
   using tarch::la::raw;
   ClosestMesh closest = geoInterface.inquireClosestMesh ( raw(inquiryCenter), ids );
   validateEquals ( closest.position(), constants::positionOnGeometry() );
}

void SolverInterfaceTestGeometry:: testVoxelQueryEpsBox()
{
  preciceTrace ( "testVoxelQueryEpsBox()" );
  using namespace tarch::la;
  SolverInterface geoInterface ( "TestAccessor", 0, 1 );
  configureSolverInterface (
      _pathToTests + "eps-box.xml", geoInterface );
  validateEquals ( geoInterface.getDimensions(), 2 );
  geoInterface.initialize();
  utils::Vector2D center(0.0);
  center[0] = 1.1;
  center[1] = 0.9;
  utils::Vector2D h(0.1);
  std::set<int> ids;
  VoxelPosition pos = geoInterface.inquireVoxelPosition (
                      raw(center), raw(h), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[0] = 1.1 - (NUMERICAL_ZERO_DIFFERENCE / 2.0);
  pos = geoInterface.inquireVoxelPosition ( raw(center), raw(h), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[1] = 1.0;
  pos = geoInterface.inquireVoxelPosition ( raw(center), raw(h), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[1] = 0.9;
  pos = geoInterface.inquireVoxelPosition ( raw(center), raw(h), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[1] = 0.5;
  pos = geoInterface.inquireVoxelPosition ( raw(center), raw(h), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[1] = 0.4;
  pos = geoInterface.inquireVoxelPosition ( raw(center), raw(h), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[1] = 0.0;
  pos = geoInterface.inquireVoxelPosition ( raw(center), raw(h), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[0] = 1.1 + (NUMERICAL_ZERO_DIFFERENCE / 2.0);
  pos = geoInterface.inquireVoxelPosition ( raw(center), raw(h), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[1] = 1.0;
  pos = geoInterface.inquireVoxelPosition ( raw(center), raw(h), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[1] = 0.9;
  pos = geoInterface.inquireVoxelPosition ( raw(center), raw(h), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[1] = 0.5;
  pos = geoInterface.inquireVoxelPosition ( raw(center), raw(h), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[1] = 0.4;
  pos = geoInterface.inquireVoxelPosition ( raw(center), raw(h), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[1] = 0.0;
  pos = geoInterface.inquireVoxelPosition ( raw(center), raw(h), false, ids );
  validateEquals ( pos.position(), constants::positionOutsideOfGeometry() );

  center[0] = 1.1 - (5.0 * NUMERICAL_ZERO_DIFFERENCE);
  pos = geoInterface.inquireVoxelPosition ( raw(center), raw(h), false, ids );
  validateEquals ( pos.position(), constants::positionOnGeometry() );

  center[1] = 1.0;
  pos = geoInterface.inquireVoxelPosition ( raw(center), raw(h), false, ids );
  validateEquals ( pos.position(), constants::positionOnGeometry() );

  center[1] = 0.9;
  pos = geoInterface.inquireVoxelPosition ( raw(center), raw(h), false, ids );
  validateEquals ( pos.position(), constants::positionOnGeometry() );

  center[1] = 0.5;
  pos = geoInterface.inquireVoxelPosition ( raw(center), raw(h), false, ids );
  validateEquals ( pos.position(), constants::positionOnGeometry() );

  center[1] = 0.4;
  pos = geoInterface.inquireVoxelPosition ( raw(center), raw(h), false, ids );
  validateEquals ( pos.position(), constants::positionOnGeometry() );

  center[1] = 0.0;
  pos = geoInterface.inquireVoxelPosition ( raw(center), raw(h), false, ids );
  validateEquals ( pos.position(), constants::positionOnGeometry() );
}

void SolverInterfaceTestGeometry:: testConservativeStationaryDataMapping()
{
  preciceTrace("testConservativeStationaryDataMapping()");
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
  preciceDebug ( "preCICE initialized");
  precice.mapWriteDataFrom(meshID);
  // Validate results
  impl::PtrParticipant p = precice._impl->_accessor;
  preciceDebug ( "Participant found");
  validate(p != NULL);
  preciceDebug ( "dataContexts: " << p->_dataContexts << " and dataID: " << dataID);
  validate(p->_dataContexts[dataID] != NULL);
  mesh::PtrData data = p->_dataContexts[dataID]->toData;
  preciceDebug ( "ToData found");
  validate(data.get() != NULL);
  utils::DynVector& writtenValues = data->values();

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
  preciceTrace("testMappingRBF()");
  SolverInterface interface("TestAccessor", 0, 1);
  configureSolverInterface(_pathToTests + "mapping-rbf.xml", interface);
  validateEquals(interface.getDimensions(), 2);
  interface.initialize();
  using utils::Vector2D;

  // Set write data for consistent thin plate splines
  {
    int solverMeshID = interface.getMeshID("SolverMesh-ConsistentTPS");
    // set positions
    Vector2D position ( 0.0, 0.0 );
    int i0 = interface.setMeshVertex ( solverMeshID, raw(position) );
    assignList(position) = 1.0, 0.0;
    int i1 = interface.setMeshVertex ( solverMeshID, raw(position) );
    assignList(position) = 1.0, 1.0;
    int i2 = interface.setMeshVertex ( solverMeshID, raw(position) );
    assignList(position) = 0.0, 1.0;
    int i3 = interface.setMeshVertex ( solverMeshID, raw(position) );
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
    Vector2D position ( 0.0, 0.0 );
    int i0 = interface.setMeshVertex ( solverMeshID, raw(position) );
    assignList(position) = 1.0, 0.0;
    int i1 = interface.setMeshVertex ( solverMeshID, raw(position) );
    assignList(position) = 1.0, 1.0;
    int i2 = interface.setMeshVertex ( solverMeshID, raw(position) );
    assignList(position) = 0.0, 1.0;
    int i3 = interface.setMeshVertex ( solverMeshID, raw(position) );
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
    Vector2D position ( 0.0, 0.0 );
    int i0 = interface.setMeshVertex ( solverMeshID, raw(position) );
    assignList(position) = 1.0, 0.0;
    int i1 = interface.setMeshVertex ( solverMeshID, raw(position) );
    assignList(position) = 1.0, 1.0;
    int i2 = interface.setMeshVertex ( solverMeshID, raw(position) );
    assignList(position) = 0.0, 1.0;
    int i3 = interface.setMeshVertex ( solverMeshID, raw(position) );
    int dataID = interface.getDataID ( "ConsistentVS", solverMeshID );
    Vector2D data ( 0.0, 0.0 );
    interface.writeVectorData ( dataID, i0, raw(data) );
    assignList(data) = 2.0, 2.0;
    interface.writeVectorData ( dataID, i1, raw(data) );
    assignList(data) = 6.0, 6.0;
    interface.writeVectorData ( dataID, i2, raw(data) );
    assignList(data) = 2.0, 2.0;
    interface.writeVectorData ( dataID, i3, raw(data) );
  }

  // Set write data for conservative volume splines
  {
    int solverMeshID = interface.getMeshID("SolverMesh-ConservativeVS");
    // set positions
    Vector2D position ( 0.0, 0.0 );
    int i0 = interface.setMeshVertex ( solverMeshID, raw(position) );
    assignList(position) = 1.0, 0.0;
    int i1 = interface.setMeshVertex ( solverMeshID, raw(position) );
    assignList(position) = 1.0, 1.0;
    int i2 = interface.setMeshVertex ( solverMeshID, raw(position) );
    assignList(position) = 0.0, 1.0;
    int i3 = interface.setMeshVertex ( solverMeshID, raw(position) );
    int dataID = interface.getDataID ( "ConservativeVS", solverMeshID );
    Vector2D data ( 0.0, 0.0 );
    interface.writeVectorData ( dataID, i0, raw(data) );
    assignList(data) = 2.0, 2.0;
    interface.writeVectorData ( dataID, i1, raw(data) );
    assignList(data) = 6.0, 6.0;
    interface.writeVectorData ( dataID, i2, raw(data) );
    assignList(data) = 2.0, 2.0;
    interface.writeVectorData ( dataID, i3, raw(data) );
  }

  interface.advance(1.0);
  interface.finalize();
}

//void SolverInterfaceTestGeometry:: testConservativeDataMapping ()
//{
//# if defined(Dim2)
//  preciceTrace ( "testConservativeDataMapping()" );
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
  preciceTrace("testPinelli()");
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
  int dataIsEOFVeloc = interface.getDataID("Velocities",meshIdEOF);
  int dataIsEOFForces = interface.getDataID("Forces",meshIdEOF);

  using utils::Vector2D;
  std::vector<int> vertexIdsEOF;

  Vector2D position ( 0.15, 0.15 );
  int vertexId = interface.setMeshVertex ( meshIdEOF, raw(position) );
  vertexIdsEOF.push_back(vertexId);

  assignList(position) = 0.15, 0.2;
  vertexId = interface.setMeshVertex ( meshIdEOF, raw(position) );
  vertexIdsEOF.push_back(vertexId);

  assignList(position) = 0.15, 0.25;
  vertexId = interface.setMeshVertex ( meshIdEOF, raw(position) );
  vertexIdsEOF.push_back(vertexId);

  assignList(position) = 0.2, 0.15;
  vertexId = interface.setMeshVertex ( meshIdEOF, raw(position) );
  vertexIdsEOF.push_back(vertexId);

  assignList(position) = 0.2, 0.25;
  vertexId = interface.setMeshVertex ( meshIdEOF, raw(position) );
  vertexIdsEOF.push_back(vertexId);

  assignList(position) = 0.25, 0.15;
  vertexId = interface.setMeshVertex ( meshIdEOF, raw(position) );
  vertexIdsEOF.push_back(vertexId);

  assignList(position) = 0.25, 0.2;
  vertexId = interface.setMeshVertex ( meshIdEOF, raw(position) );
  vertexIdsEOF.push_back(vertexId);

  assignList(position) = 0.25, 0.25;
  vertexId = interface.setMeshVertex ( meshIdEOF, raw(position) );
  vertexIdsEOF.push_back(vertexId);

  interface.initialize();

  for (int vertexID : vertexIdsEOF){
    double data[2] = {2.0,2.0};
    interface.writeVectorData(dataIsEOFVeloc,vertexID,data);
  }

  interface.advance(0.1);

  Vector2D totalForce ( 0.0, 0.0 );

  for (int vertexID : vertexIdsEOF){
    double data[2] = {0.0,0.0};
    interface.readVectorData(dataIsEOFForces,vertexID,data);
    totalForce[0] += data[0];
    totalForce[1] += data[1];
  }

  validateNumericalEquals(totalForce[0]+totalForce[1], 100.0);
}


void SolverInterfaceTestGeometry:: testCustomGeometryCreation()
{
  preciceTrace ( "testCustomGeometryCreation()" );
  using tarch::la::wrap;
  using tarch::la::raw;
  { // 2D
    using utils::Vector2D;
    SolverInterface geo ( "TestAccessor", 0, 1 );
    configureSolverInterface (
        _pathToTests + "solvermesh-2D.xml", geo );
    validateEquals ( geo.getDimensions(), 2 );
    std::string meshName = "custom-geometry";
    int meshID = geo.getMeshID ( meshName );
    geo._impl->_accessor->meshContext(meshID).meshRequirement =
        mapping::Mapping::FULL;
    Vector2D coords0(0.0, 0.0);
    Vector2D coords1(1.0, 0.0);
    Vector2D coords2(1.0, 1.0);
    Vector2D coords3(0.0, 1.0);
    int v0 = geo.setMeshVertex ( meshID, raw(coords0) );
    int v1 = geo.setMeshVertex ( meshID, raw(coords1) );
    int v2 = geo.setMeshVertex ( meshID, raw(coords2) );
    int v3 = geo.setMeshVertex ( meshID, raw(coords3) );
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
    Vector2D coords;
    coords = wrap<2,double>(vertexIter.vertexCoords());
    validateWithMessage (equals(coords, coords0), coords);
    vertexIter++;
    validateEquals (vertexIter.vertexID(), 1);
    coords = wrap<2,double>(vertexIter.vertexCoords());
    validateWithMessage (equals(coords, coords1), coords);
    vertexIter++;
    validateEquals (vertexIter.vertexID(), 2);
    coords = wrap<2,double>(vertexIter.vertexCoords());
    validateWithMessage (equals(coords, coords2), coords);
    vertexIter++;
    validateEquals (vertexIter.vertexID(), 3);
    coords = wrap<2,double>(vertexIter.vertexCoords());
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
    using utils::Vector3D;
    SolverInterface geo ( "TestAccessor", 0, 1 );
    configureSolverInterface (
        _pathToTests + "solvermesh-3D.xml", geo );
    validateEquals ( geo.getDimensions(), 3 );
    std::string meshName = "custom-geometry";
    int meshID = geo.getMeshID ( meshName );
    geo._impl->_accessor->meshContext(meshID).meshRequirement =
        mapping::Mapping::FULL;
    Vector3D coords0(0.0,  0.0, 0.0);
    Vector3D coords1(0.5, -0.5, 0.5);
    Vector3D coords2(1.0,  0.0, 1.0);
    Vector3D coords3(0.5,  0.5, 0.5);
    int v0 = geo.setMeshVertex ( meshID, raw(coords0) );
    int v1 = geo.setMeshVertex ( meshID, raw(coords1) );
    int v2 = geo.setMeshVertex ( meshID, raw(coords2) );
    int v3 = geo.setMeshVertex ( meshID, raw(coords3) );
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
    Vector3D coords;
    coords = wrap<3,double>(vertexIter.vertexCoords());
    validateWithMessage (equals(coords, coords0), coords);
    vertexIter++;
    validateEquals (vertexIter.vertexID(), 1);
    coords = wrap<3,double>(vertexIter.vertexCoords());
    validateWithMessage (equals(coords, coords1), coords);
    vertexIter++;
    validateEquals (vertexIter.vertexID(), 2);
    coords = wrap<3,double>(vertexIter.vertexCoords());
    validateWithMessage (equals(coords, coords2), coords);
    vertexIter++;
    validateEquals (vertexIter.vertexID(), 3);
    coords = wrap<3,double>(vertexIter.vertexCoords());
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
    using utils::Vector3D;
    SolverInterface geo ( "TestAccessor", 0, 1 );
    configureSolverInterface (
        _pathToTests + "solvermesh-3D.xml", geo );
    validateEquals ( geo.getDimensions(), 3 );
    std::string meshName = "custom-geometry";
    int meshID = geo.getMeshID ( meshName );
    geo._impl->_accessor->meshContext(meshID).meshRequirement =
        mapping::Mapping::FULL;
    Vector3D coords0(0.0,  0.0, 0.0);
    Vector3D coords1(0.5, -0.5, 0.5);
    Vector3D coords2(1.0,  0.0, 1.0);
    Vector3D coords3(0.5,  0.5, 0.5);
    int v0 = geo.setMeshVertex ( meshID, raw(coords0) );
    int v1 = geo.setMeshVertex ( meshID, raw(coords1) );
    int v2 = geo.setMeshVertex ( meshID, raw(coords2) );
    int v3 = geo.setMeshVertex ( meshID, raw(coords3) );
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
    Vector3D coords;
    coords = wrap<3,double>(vertexIter.vertexCoords());
    validateWithMessage (equals(coords, coords0), coords);
    vertexIter++;
    validateEquals (vertexIter.vertexID(), 1);
    coords = wrap<3,double>(vertexIter.vertexCoords());
    validateWithMessage (equals(coords, coords1), coords);
    vertexIter++;
    validateEquals (vertexIter.vertexID(), 2);
    coords = wrap<3,double>(vertexIter.vertexCoords());
    validateWithMessage (equals(coords, coords2), coords);
    vertexIter++;
    validateEquals (vertexIter.vertexID(), 3);
    coords = wrap<3,double>(vertexIter.vertexCoords());
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
  preciceTrace("testBug()");

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
  utils::Vector3D center ( 0.0, 0.0, 0.0 );
  utils::Vector3D h ( 0.0166666666666666 );
  using tarch::la::raw;
  VoxelPosition pos;
  pos = interface.inquireVoxelPosition ( raw(center), raw(h), false, meshIDs );
  validateEquals ( pos.position(), constants::positionOnGeometry() );
  assign(center) = -0.01666666666666666;
  int pointPos = interface.inquirePosition ( raw(center), meshIDs );
  validateEquals ( pointPos, constants::positionInsideOfGeometry() );
}

void SolverInterfaceTestGeometry:: testBug2()
{
  preciceTrace("testBug2()");
  SolverInterface interface("Participant-testBug2", 0, 1);
  configureSolverInterface(_pathToTests + "testBug2.xml", interface);
  validateEquals ( interface.getDimensions(), 3 );
  std::set<int> meshIDs = interface.getMeshIDs();
  interface.initialize();
  utils::Vector3D center(7.6875, 0.7625, 3.7125);
  utils::Vector3D h(0.0125);
  using tarch::la::raw;
  VoxelPosition pos;
  pos = interface.inquireVoxelPosition ( raw(center), raw(h), false, meshIDs );
  validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
  assignList(center) = 7.675, 0.75, 3.7;
  assign(h) = 0.025;
  pos = interface.inquireVoxelPosition ( raw(center), raw(h), false, meshIDs );
  validateEquals ( pos.position(), constants::positionInsideOfGeometry() );
  int pointPos = interface.inquirePosition ( raw(center), meshIDs );
  validateEquals ( pointPos, constants::positionInsideOfGeometry() );
}

void SolverInterfaceTestGeometry:: testBug3()
{
  preciceTrace("testBug3()");
  SolverInterface interface("Peano", 0, 1);
  configureSolverInterface(_pathToTests + "testBug3.xml", interface);
  validateEquals(interface.getDimensions(), 3);
  std::set<int> meshIDs = interface.getMeshIDs();
  interface.initialize();
  using namespace tarch::la;

  // This query was buggy
  utils::Vector3D center(0.22549, 0.5, 0.205882);
  utils::Vector3D h(0.00980392, 0.00980392, 0.00980392);
  VoxelPosition pos;
  pos = interface.inquireVoxelPosition(raw(center), raw(h), false, meshIDs);
  validateEquals(pos.position(), constants::positionInsideOfGeometry());
//  interface.exportMesh("testBug3");
}

void SolverInterfaceTestGeometry:: testBug4()
{
  preciceTrace("testBug4()");
  SolverInterface interface("Peano", 0, 1);
  configureSolverInterface(_pathToTests + "testBug4.xml", interface);
  validateEquals(interface.getDimensions(), 3);
  std::set<int> meshIDs = interface.getMeshIDs();
  interface.initialize();
  interface.exportMesh("testBug4-init");
  using namespace tarch::la;

  // Completely inside in Peano fluid domain as answered by preCICE
  utils::Vector3D center(0.275, 0.525, 0.275);
  utils::Vector3D h(0.025, 0.025, 0.025);
  VoxelPosition pos;
  pos = interface.inquireVoxelPosition(raw(center), raw(h), false, meshIDs);
  validateEquals(pos.position(), constants::positionInsideOfGeometry());
  //interface.exportMesh("testBug4");
}

void SolverInterfaceTestGeometry:: testBug5()
{
  preciceTrace("testBug5()");

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
  using namespace tarch::la;

  // On geometry in Peano fluid domain as answered wrongly by preCICE
  utils::Vector3D point(4.7777777777777777, 4.06666666666666666, 4.0666666666666666);
  int pos = interface.inquirePosition(raw(point), meshIDs);
  validateEquals(pos, constants::positionOutsideOfGeometry());

  utils::Vector3D h(0.3555555555555555);
  VoxelPosition voxelPos = interface.inquireVoxelPosition(raw(point), raw(h), false, meshIDs);

  assign(h) = 3.950617283950000e-2;
  assignList(point) = 4.343209876543000, 4.0666666666666666, 4.106172839509999;
//  precicePrint("----------------------------- START");
  voxelPos = interface.inquireVoxelPosition(raw(point), raw(h), false, meshIDs);
//  precicePrint("----------------------------- END, pos = " << voxelPos.position()
//               << ", ids.size = " << voxelPos.meshIDs().size());
  //mesh::Mesh found("Found", 3, false);
//  EdgeIterator it = voxelPos.contentHandle().edges().begin();
//  for (;it != voxelPos.contentHandle().edges().end(); it++){
//    precicePrint("Edge from " << wrap<3>(it.vertexCoords(0)) << " to " <<
//                 wrap<3>(it.vertexCoords(1)));
//    mesh::Vertex& v0 = found.createVertex(wrap<3>(it.vertexCoords(0)));
//    mesh::Vertex& v1 = found.createVertex(wrap<3>(it.vertexCoords(1)));
//    found.createEdge(v0, v1);
//  }
  //foreach (mesh::Edge& edge, voxelPos.contentHandle().triangles()){
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
//  preciceTrace("testBug6()");
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
//  VoxelPosition voxelPos = interface.inquireVoxelPosition(raw(point), raw(h), false, meshIDs);
//  validateEquals(voxelPos.position(), constants::positionOnGeometry());
//  validateEquals(voxelPos.meshIDs().size(), 0);
//}

void SolverInterfaceTestGeometry:: testUpdateSpacetree()
{
  preciceTrace ( "testUpdateSpacetree()" );
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

    foriter (VertexIterator, iter, handle.vertices()){
      interface.writeVectorData(dataID, iter.vertexID(), value);
      interface.inquirePosition(iter.vertexCoords(), std::set<int>());
    }
    interface.exportMesh(filename.str() + "1");
    interface.advance(1.0);

    foriter (VertexIterator, iter, handle.vertices()){
      interface.writeVectorData(dataID, iter.vertexID(), value);
      interface.inquirePosition(iter.vertexCoords(), std::set<int>());
    }
    interface.exportMesh(filename.str() + "2");
    interface.advance(1.0);

    foriter (VertexIterator, iter, handle.vertices()){
      interface.writeVectorData(dataID, iter.vertexID(), value);
      interface.inquirePosition(iter.vertexCoords(), std::set<int>());
    }
    interface.exportMesh(filename.str() + "3");
    interface.advance(1.0);

    foriter (VertexIterator, iter, handle.vertices()){
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
  preciceTrace("testMultipleMeshSpacetree()");
  using namespace tarch::la;
  { // Tests A: first mesh no spacetree, second, third spacetree
    SolverInterface interface("Accessor", 0, 1);
    impl::SolverInterfaceImpl* impl = interface._impl;
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
    std::vector<int> expected(3);
    expected[0] = impl->markedQueryDirectly();
    expected[1] = impl->markedQuerySpacetree();
    expected[2] = impl->markedSkip();
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(wrap<3>(&outputIDs[0]), expected), wrap<3>(&outputIDs[0]));

    expected[0] = impl->markedQueryDirectly();
    expected[1] = impl->markedQuerySpacetree();
    expected[2] = impl->markedSkip();
    std::fill(outputIDs.begin(), outputIDs.end(), -1);
    impl->selectInquiryMeshIDs(allMeshIDs, outputIDs);
    validateWithMessage(equals(wrap<3>(&outputIDs[0]), expected), wrap<3>(&outputIDs[0]));

    inputIDs.insert(0);
    expected[0] = impl->markedQueryDirectly();
    expected[1] = impl->markedSkip();
    expected[2] = impl->markedSkip();
    std::fill(outputIDs.begin(), outputIDs.end(), -1);
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(wrap<3>(&outputIDs[0]), expected), wrap<3>(&outputIDs[0]));

    inputIDs.clear();
    inputIDs.insert(1);
    expected[0] = impl->markedSkip();
    expected[1] = impl->markedQueryDirectly();
    expected[2] = impl->markedSkip();
    std::fill(outputIDs.begin(), outputIDs.end(), -1);
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(wrap<3>(&outputIDs[0]), expected), wrap<3>(&outputIDs[0]));

    inputIDs.clear();
    inputIDs.insert(2);
    expected[0] = impl->markedSkip();
    expected[1] = impl->markedSkip();
    expected[2] = impl->markedQueryDirectly();
    boost::fill(outputIDs, -1);
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(wrap<3>(&outputIDs[0]), expected), wrap<3>(&outputIDs[0]));

    inputIDs.clear();
    inputIDs.insert(0);
    inputIDs.insert(2);
    expected[0] = impl->markedQueryDirectly();
    expected[1] = impl->markedSkip();
    expected[2] = impl->markedQueryDirectly();
    boost::fill(outputIDs, -1);
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(wrap<3>(&outputIDs[0]), expected), wrap<3>(&outputIDs[0]));

    inputIDs.clear();
    inputIDs.insert(1);
    inputIDs.insert(2);
    expected[0] = impl->markedSkip();
    expected[1] = impl->markedQuerySpacetree();
    expected[2] = impl->markedSkip();
    boost::fill(outputIDs, -1);
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(wrap<3>(&outputIDs[0]), expected), wrap<3>(&outputIDs[0]));

    inputIDs.clear();
    inputIDs.insert(0);
    inputIDs.insert(1);
    expected[0] = impl->markedQueryDirectly();
    expected[1] = impl->markedQueryDirectly();
    expected[2] = impl->markedSkip();
    boost::fill(outputIDs, -1);
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(wrap<3>(&outputIDs[0]), expected), wrap<3>(&outputIDs[0]));
  }
  { // Tests B: first mesh has spacetree, second, third not
    SolverInterface interface("Accessor", 0, 1);
    impl::SolverInterfaceImpl* impl = interface._impl;
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
    std::vector<int> expected(3);
    expected[0] = impl->markedQuerySpacetree();
    expected[1] = impl->markedQueryDirectly();
    expected[2] = impl->markedQueryDirectly();
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(wrap<3>(&outputIDs[0]), expected), wrap<3>(&outputIDs[0]));

    expected[0] = impl->markedQuerySpacetree();
    expected[1] = impl->markedQueryDirectly();
    expected[2] = impl->markedQueryDirectly();
    std::fill(outputIDs.begin(), outputIDs.end(), -1);
    impl->selectInquiryMeshIDs(allMeshIDs, outputIDs);
    validateWithMessage(equals(wrap<3>(&outputIDs[0]), expected), wrap<3>(&outputIDs[0]));

    inputIDs.insert(0);
    expected[0] = impl->markedQuerySpacetree();
    expected[1] = impl->markedSkip();
    expected[2] = impl->markedSkip();
    std::fill(outputIDs.begin(), outputIDs.end(), -1);
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(wrap<3>(&outputIDs[0]), expected), wrap<3>(&outputIDs[0]));

    inputIDs.clear();
    inputIDs.insert(1);
    expected[0] = impl->markedSkip();
    expected[1] = impl->markedQueryDirectly();
    expected[2] = impl->markedSkip();
    std::fill(outputIDs.begin(), outputIDs.end(), -1);
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(wrap<3>(&outputIDs[0]), expected), wrap<3>(&outputIDs[0]));

    inputIDs.clear();
    inputIDs.insert(2);
    expected[0] = impl->markedSkip();
    expected[1] = impl->markedSkip();
    expected[2] = impl->markedQueryDirectly();
    boost::fill(outputIDs, -1);
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(wrap<3>(&outputIDs[0]), expected), wrap<3>(&outputIDs[0]));

    inputIDs.clear();
    inputIDs.insert(0);
    inputIDs.insert(2);
    expected[0] = impl->markedQuerySpacetree();
    expected[1] = impl->markedSkip();
    expected[2] = impl->markedQueryDirectly();
    boost::fill(outputIDs, -1);
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(wrap<3>(&outputIDs[0]), expected), wrap<3>(&outputIDs[0]));

    inputIDs.clear();
    inputIDs.insert(1);
    inputIDs.insert(2);
    expected[0] = impl->markedSkip();
    expected[1] = impl->markedQueryDirectly();
    expected[2] = impl->markedQueryDirectly();
    boost::fill(outputIDs, -1);
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(wrap<3>(&outputIDs[0]), expected), wrap<3>(&outputIDs[0]));

    inputIDs.clear();
    inputIDs.insert(0);
    inputIDs.insert(1);
    expected[0] = impl->markedQuerySpacetree();
    expected[1] = impl->markedQueryDirectly();
    expected[2] = impl->markedSkip();
    boost::fill(outputIDs, -1);
    impl->selectInquiryMeshIDs(inputIDs, outputIDs);
    validateWithMessage(equals(wrap<3>(&outputIDs[0]), expected), wrap<3>(&outputIDs[0]));
  }
}

}} // namespace precice, tests


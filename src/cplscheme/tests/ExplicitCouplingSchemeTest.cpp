// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ExplicitCouplingSchemeTest.hpp"
#include "cplscheme/ExplicitCouplingScheme.hpp"
#include "cplscheme/ImplicitCouplingScheme.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"
#include "cplscheme/Constants.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "geometry/config/GeometryConfiguration.hpp"
#include "com/Communication.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "com/config/CommunicationConfiguration.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"
#include "utils/xml/XMLTag.hpp"
#include "tarch/la/WrappedVector.h"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::cplscheme::tests::ExplicitCouplingSchemeTest)

namespace precice {
namespace cplscheme {
namespace tests {

using utils::Vector3D;

tarch::logging::Log ExplicitCouplingSchemeTest::
   _log ( "precice::cplscheme::tests::ExplicitCouplingSchemeTest" );

ExplicitCouplingSchemeTest:: ExplicitCouplingSchemeTest ()
:
   TestCase ( "cplscheme::ExplicitCouplingSchemeTest" ),
   _pathToTests()
{}

void ExplicitCouplingSchemeTest:: setUp ()
{
  _pathToTests = utils::Globals::getPathToSources() + "/cplscheme/tests/";
}

void ExplicitCouplingSchemeTest:: run ()
{
# ifndef PRECICE_NO_MPI
  typedef utils::Parallel Par;
  if ( Par::getCommunicatorSize() > 1 ) {
    std::vector<int> ranks;
    ranks += 0, 1;
    // Only use process 0 and 1 for the following tests
    Par::Communicator comm = Par::getRestrictedCommunicator(ranks);
    if ( Par::getProcessRank() <= 1 ){
      Par::setGlobalCommunicator(comm) ;
      validateEquals ( Par::getCommunicatorSize(), 2 );
      testMethod ( testSimpleExplicitCoupling );
      testMethod ( testConfiguredSimpleExplicitCoupling );
      testMethod ( testExplicitCouplingFirstParticipantSetsDt );
      //testMethod ( testExplicitCouplingSecondParticipantSetsDt );
      testMethod ( testExplicitCouplingWithSubcycling );
      testMethod ( testConfiguredExplicitCouplingWithSubcycling );
      Par::setGlobalCommunicator(Par::getCommunicatorWorld());
    }
  }
# endif // not PRECICE_NO_MPI
}

#ifndef PRECICE_NO_MPI

void ExplicitCouplingSchemeTest:: testSimpleExplicitCoupling()
{
  preciceTrace ( "testSimpleExplicitiCoupling()" );
  utils::Parallel::synchronizeProcesses();

  mesh::PropertyContainer::resetPropertyIDCounter();
  utils::XMLTag root = utils::getRootTag();
  mesh::PtrDataConfiguration dataConfig ( new mesh::DataConfiguration(root) );
  dataConfig->addData ( "data0", 1 );
  dataConfig->addData ( "data1", 3 );
  mesh::MeshConfiguration meshConfig ( root, dataConfig );
  mesh::PtrMesh mesh ( new mesh::Mesh("mesh", 3, false) );
  mesh->createData ( "data0", 1 );
  mesh->createData ( "data1", 3 );
  mesh->createVertex ( Vector3D(0.0) );
  mesh->allocateDataValues ();
  meshConfig.addMesh ( mesh );

  com::PtrCommunication communication ( new com::MPIDirectCommunication() );
  std::string nameParticipant0 ( "participant0" );
  std::string nameParticipant1 ( "participant1" );
  double maxTime = 1.0;
  int maxTimesteps = 10;
  double timestepLength = 0.1;
  std::string localParticipant ( "" );
  int sendDataIndex = -1;
  int receiveDataIndex = -1;
  if ( utils::Parallel::getProcessRank() == 0 ) {
    localParticipant = nameParticipant0;
    sendDataIndex = 0;
    receiveDataIndex = 1;
  }
  else if ( utils::Parallel::getProcessRank() == 1 ) {
    localParticipant = nameParticipant1;
    sendDataIndex = 1;
    receiveDataIndex = 0;
  }
  constants::TimesteppingMethod dtMethod = constants::FIXED_DT;
  cplscheme::ExplicitCouplingScheme cplScheme (
    maxTime, maxTimesteps, timestepLength, 12, nameParticipant0,
    nameParticipant1, localParticipant, communication, dtMethod );
  cplScheme.addDataToSend ( mesh->data()[sendDataIndex], false );
  cplScheme.addDataToReceive ( mesh->data()[receiveDataIndex], false );
  connect ( nameParticipant0, nameParticipant1, localParticipant, communication );
  runSimpleExplicitCoupling ( cplScheme, localParticipant, meshConfig );
}

void ExplicitCouplingSchemeTest:: testConfiguredSimpleExplicitCoupling ()
{
  preciceTrace ( "testConfiguredSimpleExplicitCoupling()" );
  using namespace mesh;
  utils::Parallel::synchronizeProcesses ();
  assertion ( utils::Parallel::getCommunicatorSize() > 1 );
  mesh::PropertyContainer::resetPropertyIDCounter ();

  std::string configurationPath ( _pathToTests + "explicit-coupling-scheme-1.xml" );

  std::string localParticipant ( "" );
  if ( utils::Parallel::getProcessRank() == 0 ) {
    localParticipant = "participant0";
  }
  else if ( utils::Parallel::getProcessRank() == 1 ) {
    localParticipant = "participant1";
  }
  utils::XMLTag root = utils::getRootTag();
  PtrDataConfiguration dataConfig(new DataConfiguration(root));
  dataConfig->setDimensions(3);
  PtrMeshConfiguration meshConfig(new MeshConfiguration(root, dataConfig));
  meshConfig->setDimensions(3);
  com::PtrCommunicationConfiguration comConfig(new com::CommunicationConfiguration(root));
  geometry::GeometryConfiguration geoConfig(root, meshConfig);
  geoConfig.setDimensions(3);
  CouplingSchemeConfiguration cplSchemeConfig(root, meshConfig, comConfig);

  utils::configure(root, configurationPath);
  //validate(success);
  //validate(dataConfig->isValid());
  //validate(meshConfig->isValid());
  //validate(comConfig->isValid());
  //validate(geoConfig.isValid());
  //validate(cplSchemeConfig.isValid());
  meshConfig->setMeshSubIDs();
  com::PtrCommunication com = comConfig->getCommunication("participant0", "participant1");

  geoConfig.geometries()[0]->create ( *meshConfig->meshes()[0] );
  connect ( "participant0", "participant1", localParticipant, com );
  runSimpleExplicitCoupling ( *cplSchemeConfig.getCouplingScheme(localParticipant),
                              localParticipant, *meshConfig );
}

void ExplicitCouplingSchemeTest:: testExplicitCouplingFirstParticipantSetsDt()
{
  preciceTrace ( "testExplicitCouplingFirstParticipantSetsDt()" );
  using namespace mesh;
  utils::Parallel::synchronizeProcesses();
  std::string configurationPath ( _pathToTests + "explicit-coupling-scheme-2.xml" );
  std::string localParticipant ( "" );
  if ( utils::Parallel::getProcessRank() == 0 ){
    localParticipant = "participant0";
  }
  else if ( utils::Parallel::getProcessRank() == 1 ){
    localParticipant = "participant1";
  }

  utils::XMLTag root = utils::getRootTag();
  PtrDataConfiguration dataConfig(new DataConfiguration(root));
  dataConfig->setDimensions(3);
  PtrMeshConfiguration meshConfig(new MeshConfiguration(root, dataConfig));
  meshConfig->setDimensions(3);
  com::PtrCommunicationConfiguration comConfig(new com::CommunicationConfiguration(root));
  geometry::GeometryConfiguration geoConfig(root, meshConfig);
  geoConfig.setDimensions(3);
  CouplingSchemeConfiguration cplSchemeConfig(root, meshConfig, comConfig);

  utils::configure(root, configurationPath);
  //validate(success);
  //validate(dataConfig->isValid());
  //validate(meshConfig->isValid());
  //validate(comConfig->isValid());
  //validate(geoConfig.isValid());
  //validate(cplSchemeConfig.isValid());
  meshConfig->setMeshSubIDs();
  com::PtrCommunication com = comConfig->getCommunication("participant0", "participant1");

  geoConfig.geometries()[0]->create( *meshConfig->meshes()[0] );
  connect ( "participant0", "participant1", localParticipant, com );
  CouplingScheme& cplScheme = *cplSchemeConfig.getCouplingScheme(localParticipant);

  double computedTime = 0.0;
  int computedTimesteps = 0;
  if ( localParticipant == std::string("participant0") ){
    double dt = 0.3;
    cplScheme.initialize ( 0.0, 0 );
    validateEquals ( cplScheme.isCouplingTimestepComplete(), false );
    validateEquals ( cplScheme.isCouplingOngoing(), true );
    while ( cplScheme.isCouplingOngoing() ){
      computedTime += dt;
      computedTimesteps++;
      cplScheme.addComputedTime(dt);
      cplScheme.advance();
      validate ( cplScheme.isCouplingTimestepComplete() );
      validateNumericalEquals ( computedTime, cplScheme.getTime() );
      validateEquals ( computedTimesteps, cplScheme.getTimesteps() );
      if ( cplScheme.isCouplingOngoing() ){
        validate ( cplScheme.hasDataBeenExchanged() );
      }
    }
    cplScheme.finalize();
    validateNumericalEquals ( computedTime, 1.2 );
    validateEquals ( computedTimesteps, 4 );
    validateEquals ( cplScheme.isCouplingTimestepComplete(), true );
    validateEquals ( cplScheme.isCouplingOngoing(), false );
  }
  else {
    assertion1 ( localParticipant == std::string("participant1"), localParticipant );
    cplScheme.initialize ( 0.0, 0 );
    validateEquals ( cplScheme.isCouplingTimestepComplete(), false );
    validateEquals ( cplScheme.isCouplingOngoing(), true );
    while ( cplScheme.isCouplingOngoing() ){
      computedTime += cplScheme.getTimestepLength();
      computedTimesteps++;
      cplScheme.addComputedTime(cplScheme.getTimestepLength());
      cplScheme.advance();
      validate ( cplScheme.isCouplingTimestepComplete() );
      validateNumericalEquals ( computedTime, cplScheme.getTime() );
      validateEquals ( computedTimesteps, cplScheme.getTimesteps() );
      if ( cplScheme.isCouplingOngoing() ){
        validate ( cplScheme.hasDataBeenExchanged() );
      }
    }
    cplScheme.finalize();
    validateNumericalEquals ( computedTime, 1.2 );
    validateEquals ( computedTimesteps, 4 );
    validateEquals ( cplScheme.isCouplingTimestepComplete(), true );
    validateEquals ( cplScheme.isCouplingOngoing(), false );
  }
}

//void ExplicitCouplingSchemeTest:: testExplicitCouplingSecondParticipantSetsDt()
//{
//  preciceTrace ( "testExplicitCouplingSecondParticipantSetsDt()" );
//  using namespace mesh;
//  utils::Parallel::synchronizeProcesses();
//  std::string configurationPath ( _pathToTests + "explicit-coupling-scheme-2.xml" );
//  std::string localParticipant ( "" );
//  if ( utils::Parallel::getProcessRank() == 0 ){
//    localParticipant = "participant0";
//  }
//  else if ( utils::Parallel::getProcessRank() == 1 ){
//    localParticipant = "participant1";
//  }
//
//  PtrDataConfiguration dataConfig ( new DataConfiguration(3) );
//  utils::configure ( dataConfig, configurationPath );
//  validate ( dataConfig->isValid() );
//
//  PtrMeshConfiguration meshConfig ( new MeshConfiguration(dataConfig) );
//  utils::configure ( meshConfig, configurationPath );
//  validate ( meshConfig->isValid() );
//  meshConfig->setMeshSubIDs();
//
//  com::PtrCommunicationConfiguration
//  comConfig ( new com::CommunicationConfiguration() );
//  utils::configure ( comConfig, configurationPath );
//  validate ( comConfig->isValid() );
//  com::PtrCommunication com =
//    comConfig->getCommunication ( "participant0", "participant1" );
//
//  geometry::GeometryConfiguration geoConfig ( meshConfig, 3 );
//  utils::configure ( geoConfig, configurationPath );
//  validate ( geoConfig.isValid() );
//
//  CouplingSchemeConfiguration cplSchemeConfig ( meshConfig, comConfig );
//  utils::configure ( cplSchemeConfig, configurationPath );
//  validate ( cplSchemeConfig.isValid() );
//
//  geoConfig.geometries()[0]->create( *meshConfig->meshes()[0] );
//  connect ( "participant0", "participant1", localParticipant, com );
//  CouplingScheme& cplScheme = *cplSchemeConfig.getCouplingScheme(localParticipant);
//
//  double computedTime = 0.0;
//  int computedTimesteps = 0;
//  if ( localParticipant == std::string("participant0") ){
//    double dt = 0.3;
//    cplScheme.initialize ( 0.0, 0 );
//    validateEquals ( cplScheme.isCouplingTimestepComplete(), false );
//    validateEquals ( cplScheme.isCouplingOngoing(), true );
//    while ( cplScheme.isCouplingOngoing() ){
//      computedTime += dt;
//      computedTimesteps++;
//      cplScheme.advance ( dt );
//      validate ( cplScheme.isCouplingTimestepComplete() );
//      validateNumericalEquals ( computedTime, cplScheme.getTime() );
//      validateEquals ( computedTimesteps, cplScheme.getTimesteps() );
//      if ( cplScheme.isCouplingOngoing() ){
//        validate ( cplScheme.hasDataBeenExchanged() );
//      }
//    }
//    cplScheme.finalize();
//    validateNumericalEquals ( computedTime, 1.2 );
//    validateEquals ( computedTimesteps, 4 );
//    validateEquals ( cplScheme.isCouplingTimestepComplete(), true );
//    validateEquals ( cplScheme.isCouplingOngoing(), false );
//  }
//  else {
//    assertion1 ( localParticipant == std::string("participant1"), localParticipant );
//    cplScheme.initialize ( 0.0, 0 );
//    validateEquals ( cplScheme.isCouplingTimestepComplete(), false );
//    validateEquals ( cplScheme.isCouplingOngoing(), true );
//    while ( cplScheme.isCouplingOngoing() ){
//      computedTime += cplScheme.getTimestepLength();
//      computedTimesteps++;
//      cplScheme.advance ( cplScheme.getTimestepLength() );
//      validate ( cplScheme.isCouplingTimestepComplete() );
//      validateNumericalEquals ( computedTime, cplScheme.getTime() );
//      validateEquals ( computedTimesteps, cplScheme.getTimesteps() );
//      if ( cplScheme.isCouplingOngoing() ){
//        validate ( cplScheme.hasDataBeenExchanged() );
//      }
//    }
//    cplScheme.finalize();
//    validateNumericalEquals ( computedTime, 1.2 );
//    validateEquals ( computedTimesteps, 4 );
//    validateEquals ( cplScheme.isCouplingTimestepComplete(), true );
//    validateEquals ( cplScheme.isCouplingOngoing(), false );
//  }
//}

void ExplicitCouplingSchemeTest:: runSimpleExplicitCoupling
(
  CouplingScheme&                cplScheme,
  const std::string&             participantName,
  const mesh::MeshConfiguration& meshConfig )
{
  preciceTrace1 ( "runSimpleExplicitCoupling()", participantName );
  using boost::get;

  validateEquals ( meshConfig.meshes().size(), 1 );
  mesh::PtrMesh mesh = meshConfig.meshes()[0];
  validateEquals ( mesh->data().size(), 2 );
  utils::DynVector& dataValues0 = mesh->data()[0]->values();
  utils::DynVector& dataValues1 = mesh->data()[1]->values();
  validate ( mesh->vertices().size() > 0 );
  mesh::Vertex& vertex = mesh->vertices()[0];
  double valueData0 = 1.0;
  Vector3D valueData1 ( 1.0 );

  double computedTime = 0.0;
  int computedTimesteps = 0;

  if ( participantName == std::string("participant0") ) {
    cplScheme.initialize ( 0.0, 0 );
    validate ( not cplScheme.hasDataBeenExchanged() );
    validateEquals ( cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()), false );
    validateEquals ( cplScheme.isActionRequired(constants::actionReadIterationCheckpoint()), false );
    validateEquals ( cplScheme.isCouplingTimestepComplete(), false );
    validateEquals ( cplScheme.isCouplingOngoing(), true );
    while ( cplScheme.isCouplingOngoing() ) {
      dataValues0[vertex.getID()] = valueData0;
      computedTime += cplScheme.getNextTimestepMaxLength();
      computedTimesteps ++;
      cplScheme.addComputedTime ( cplScheme.getNextTimestepMaxLength() );
      cplScheme.advance();
      validate ( cplScheme.isCouplingTimestepComplete() );
      validateNumericalEquals ( computedTime, cplScheme.getTime() );
      validateEquals ( computedTimesteps, cplScheme.getTimesteps() );
      validateEquals ( cplScheme.isActionRequired("WriteIterationCheckpoint"),
                       false );
      validateEquals ( cplScheme.isActionRequired("ReadIterationCheckpoint"),
                       false );
      validateEquals ( cplScheme.isCouplingTimestepComplete(), true );
      if ( cplScheme.isCouplingOngoing() ) {
        // No receive takes place for the participant that has started the
        // coupled simulation, in the last advance call
        Vector3D value;
        assign(value) = tarch::la::slice<3>(dataValues1, vertex.getID() * 3);
        validate ( tarch::la::equals(value, valueData1) );
      }
      validate ( cplScheme.hasDataBeenExchanged() );
      // Increment data values, to test if send/receive operations are also
      // correct in following timesteps.
      valueData0 += 1.0;
      valueData1 += Vector3D ( 1.0 );
    }
    cplScheme.finalize();
    // Validate results
    validateNumericalEquals ( computedTime, 1.0 );
    validateEquals ( computedTimesteps, 10 );
    validateEquals (
      cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"),
      false );
    validateEquals (
      cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"),
      false );
    validateEquals ( cplScheme.isCouplingTimestepComplete(), true );
    validateEquals ( cplScheme.isCouplingOngoing(), false );
    validate ( cplScheme.getNextTimestepMaxLength() > 0.0 );
  }
  else if ( participantName == std::string("participant1") ) {
    cplScheme.initialize ( 0.0, 0 );
    validate ( cplScheme.hasDataBeenExchanged() );
    double value = dataValues0[vertex.getID()];
    validateNumericalEquals ( value, valueData0 );
    valueData0 += 1.0;
    validateEquals (
      cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"),
      false );
    validateEquals (
      cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"),
      false );
    validateEquals ( cplScheme.isCouplingTimestepComplete(), false );
    validateEquals ( cplScheme.isCouplingOngoing(), true );
    while ( cplScheme.isCouplingOngoing() ) {
      tarch::la::slice<3>(dataValues1,vertex.getID()*3) = valueData1;
      computedTime += cplScheme.getNextTimestepMaxLength();
      computedTimesteps ++;
      cplScheme.addComputedTime ( cplScheme.getNextTimestepMaxLength() );
      cplScheme.advance();
      validateNumericalEquals ( computedTime, cplScheme.getTime() );
      validateEquals ( computedTimesteps, cplScheme.getTimesteps() );
      validateEquals (
        cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"),
        false );
      validateEquals (
        cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"),
        false );
      validateEquals ( cplScheme.isCouplingTimestepComplete(), true );
      if ( cplScheme.isCouplingOngoing() ) {
        // The participant not starting the coupled simulation does neither
        // receive nor send data in the last call to advance
        validate ( cplScheme.hasDataBeenExchanged() );
        double value = dataValues0[vertex.getID()];
        validateNumericalEquals ( value, valueData0 );
      }
      valueData0 += 1.0;
      valueData1 += Vector3D(1.0);
    }
    cplScheme.finalize ();
    // Validate results
    validateNumericalEquals ( computedTime, 1.0 );
    validateEquals ( computedTimesteps, 10 );
    validateEquals (
      cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"),
      false );
    validateEquals (
      cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"),
      false );
    validateEquals ( cplScheme.isCouplingTimestepComplete(), true );
    validateEquals ( cplScheme.isCouplingOngoing(), false );
    validate ( cplScheme.getNextTimestepMaxLength() > 0.0 );
  }
}

void ExplicitCouplingSchemeTest:: testExplicitCouplingWithSubcycling ()
{
  preciceTrace ( "testExplicitCouplingWithSubcycling()" );
  utils::Parallel::synchronizeProcesses ();
  assertion ( utils::Parallel::getCommunicatorSize() > 1 );

  mesh::PropertyContainer::resetPropertyIDCounter ();
  utils::XMLTag root = utils::getRootTag();
  mesh::PtrDataConfiguration dataConfig ( new mesh::DataConfiguration(root) );
  dataConfig->setDimensions(3);
  dataConfig->addData ( "data0", 1 );
  dataConfig->addData ( "data1", 3 );
  mesh::MeshConfiguration meshConfig ( root, dataConfig );
  meshConfig.setDimensions(3);
  mesh::PtrMesh mesh ( new mesh::Mesh("mesh", 3, false) );
  mesh->createData ( "data0", 1 );
  mesh->createData ( "data1", 3 );
  mesh->createVertex ( Vector3D(0.0) );
  mesh->allocateDataValues ();
  meshConfig.addMesh ( mesh );

  com::PtrCommunication communication ( new com::MPIDirectCommunication );
  std::string nameParticipant0 ( "participant0" );
  std::string nameParticipant1 ( "participant1" );
  double maxTime = 1.0;
  int maxTimesteps = 10;
  double timestepLength = 0.1;
  std::string localParticipant ( "" );
  int sendDataIndex = -1;
  int receiveDataIndex = -1;
  if ( utils::Parallel::getProcessRank() == 0 ) {
    localParticipant = nameParticipant0;
    sendDataIndex = 0;
    receiveDataIndex = 1;
  }
  else if ( utils::Parallel::getProcessRank() == 1 ) {
    localParticipant = nameParticipant1;
    sendDataIndex = 1;
    receiveDataIndex = 0;
  }
  constants::TimesteppingMethod dtMethod = constants::FIXED_DT;
  cplscheme::ExplicitCouplingScheme cplScheme (
    maxTime, maxTimesteps, timestepLength, 12, nameParticipant0,
    nameParticipant1, localParticipant, communication, dtMethod );
  cplScheme.addDataToSend ( mesh->data()[sendDataIndex], false );
  cplScheme.addDataToReceive ( mesh->data()[receiveDataIndex], false);
  connect ( nameParticipant0, nameParticipant1, localParticipant, communication );
  runExplicitCouplingWithSubcycling ( cplScheme, localParticipant, meshConfig );
}

void ExplicitCouplingSchemeTest:: testConfiguredExplicitCouplingWithSubcycling ()
{
  preciceTrace ( "testConfiguredExplicitCouplingWithSubcycling()" );
  using namespace mesh;
  utils::Parallel::synchronizeProcesses ();
  assertion ( utils::Parallel::getCommunicatorSize() > 1 );

  std::string configurationPath ( _pathToTests + "explicit-coupling-scheme-1.xml" );

  std::string localParticipant ( "" );
  if ( utils::Parallel::getProcessRank() == 0 ) {
    localParticipant = "participant0";
  }
  else if ( utils::Parallel::getProcessRank() == 1 ) {
    localParticipant = "participant1";
  }

  utils::XMLTag root = utils::getRootTag();
  PtrDataConfiguration dataConfig(new DataConfiguration(root));
  dataConfig->setDimensions(3);
  PtrMeshConfiguration meshConfig(new MeshConfiguration(root, dataConfig));
  meshConfig->setDimensions(3);
  com::PtrCommunicationConfiguration comConfig(new com::CommunicationConfiguration(root));
  geometry::GeometryConfiguration geoConfig(root, meshConfig);
  geoConfig.setDimensions(3);
  CouplingSchemeConfiguration cplSchemeConfig(root, meshConfig, comConfig);

  utils::configure(root, configurationPath);
  //validate(success);
  //validate(dataConfig->isValid());
  //validate(meshConfig->isValid());
  //validate(comConfig->isValid());
  //validate(geoConfig.isValid());
  //validate(cplSchemeConfig.isValid());
  meshConfig->setMeshSubIDs();
  com::PtrCommunication com = comConfig->getCommunication("participant0", "participant1");

  geoConfig.geometries()[0]->create ( *meshConfig->meshes()[0] );
  connect ( "participant0", "participant1", localParticipant, com );
  runExplicitCouplingWithSubcycling (
      *cplSchemeConfig.getCouplingScheme(localParticipant), localParticipant,
      *meshConfig );
}

void ExplicitCouplingSchemeTest:: runExplicitCouplingWithSubcycling
(
  CouplingScheme&                cplScheme,
  const std::string&             participantName,
  const mesh::MeshConfiguration& meshConfig )
{
  preciceTrace1 ( "runExplicitCouplingWithSubcycling", participantName );
  validateEquals ( meshConfig.meshes().size(), 1 );
  mesh::PtrMesh mesh = meshConfig.meshes()[0];
  validateEquals ( mesh->data().size(), 2 );
  validate ( mesh->vertices().size() > 0 );
  mesh::Vertex& vertex = mesh->vertices()[0];
  double valueData0 = 1.0;
  Vector3D valueData1(1.0);
  utils::DynVector& dataValues0 = mesh->data()[0]->values();
  utils::DynVector& dataValues1 = mesh->data()[1]->values();

  double computedTime = 0.0;
  int computedTimesteps = 0;
  std::string nameParticipant0 ( "participant0" );
  std::string nameParticipant1 ( "participant1" );
  assertion ( (participantName == nameParticipant0) ||
    (participantName == nameParticipant1) );
  if ( participantName == nameParticipant0 ) {
    cplScheme.initialize ( 0.0, 0 );
    double dtDesired = cplScheme.getNextTimestepMaxLength() / 2.0;
    double dtUsed = dtDesired;
    validate ( ! cplScheme.hasDataBeenExchanged() );
    validateEquals (
      cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"),
      false );
    validateEquals (
      cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"),
      false );
    validateEquals ( cplScheme.isCouplingTimestepComplete(), false );
    validateEquals ( cplScheme.isCouplingOngoing(), true );
    while ( cplScheme.isCouplingOngoing() ) {
      dataValues0[vertex.getID()] = valueData0;
      computedTime += dtUsed;
      computedTimesteps ++;
      cplScheme.addComputedTime(dtUsed);
      cplScheme.advance();
      // If the dt from preCICE is larger than the desired one, do subcycling,
      // else, use the dt from preCICE
      dtUsed = cplScheme.getNextTimestepMaxLength() > dtDesired
        ? dtDesired
        : cplScheme.getNextTimestepMaxLength();
      validateNumericalEquals ( computedTime, cplScheme.getTime() );
      validateEquals (
        cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"),
        false );
      validateEquals (
        cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"),
        false );
      if ( computedTimesteps % 2 == 0 ) {
        // Data exchange takes only place at every second local timestep,
        // since a subcycling of 2 is used.
        validateEquals ( cplScheme.isCouplingTimestepComplete(), true );
        if ( cplScheme.isCouplingOngoing() ) {
          // No receive takes place for the participant that has started the
          // coupled simulation, in the last advance call.
          Vector3D value;
          assign(value) = tarch::la::slice<3>(dataValues1, vertex.getID()*3);
          validate ( tarch::la::equals(value, valueData1) );
        }
        validate ( cplScheme.hasDataBeenExchanged() );
        // Increment data values, to test if send/receive operations are also
        // correct in following timesteps.
        valueData0 += 1.0;
        valueData1 += Vector3D(1.0);
      }
      else {
        validateEquals ( cplScheme.isCouplingTimestepComplete(), false );
      }
    }
    cplScheme.finalize ();
    validateNumericalEquals ( computedTime, 1.0 );
    validateEquals ( computedTimesteps, 20 );
    validateEquals (
      cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"),
      false );
    validateEquals (
      cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"),
      false );
    validateEquals ( cplScheme.isCouplingTimestepComplete(), true );
    validateEquals ( cplScheme.isCouplingOngoing(), false );
    validate ( cplScheme.getNextTimestepMaxLength() > 0.0 );
  }
  else if ( participantName == nameParticipant1 ) {
    // Start coupling
    cplScheme.initialize ( 0.0, 0 );
    // Validate current coupling status
    validate ( cplScheme.hasDataBeenExchanged() );
    validateNumericalEquals ( dataValues0[vertex.getID()], valueData0 );
    valueData0 += 1.0;
    validateEquals (
      cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"),
      false );
    validateEquals (
      cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"),
      false );
    validateEquals ( cplScheme.isCouplingTimestepComplete(), false );
    validateEquals ( cplScheme.isCouplingOngoing(), true );
    while ( cplScheme.isCouplingOngoing() ) {
      tarch::la::slice<3>(dataValues1,vertex.getID()*3) = valueData1;
      computedTime += cplScheme.getNextTimestepMaxLength ();
      computedTimesteps ++;
      cplScheme.addComputedTime ( cplScheme.getNextTimestepMaxLength() );
      cplScheme.advance();
      validateNumericalEquals ( computedTime, cplScheme.getTime() );
      validateEquals ( computedTimesteps, cplScheme.getTimesteps() );
      validateEquals (
        cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"),
        false );
      validateEquals (
        cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"),
        false );
      validateEquals ( cplScheme.isCouplingTimestepComplete(), true );
      if ( cplScheme.isCouplingOngoing() ) {
        // The participant not starting the coupled simulation does neither
        // receive nor send data in the last call to advance
        validate ( cplScheme.hasDataBeenExchanged() );
        validateNumericalEquals ( dataValues0[vertex.getID()], valueData0 );
        validate ( cplScheme.hasDataBeenExchanged() );
      }
      valueData0 += 1.0;
      valueData1 += Vector3D(1.0);
    }
    cplScheme.finalize ();
    validateNumericalEquals ( computedTime, 1.0 );
    validateEquals ( computedTimesteps, 10 );
    validateEquals (
      cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()), false );
    validateEquals (
      cplScheme.isActionRequired(constants::actionReadIterationCheckpoint()), false );
    validateEquals ( cplScheme.isCouplingTimestepComplete(), true );
    validateEquals ( cplScheme.isCouplingOngoing(), false );
    validate ( cplScheme.getNextTimestepMaxLength() > 0.0 );
  }
}

void ExplicitCouplingSchemeTest:: connect
(
  const std::string&     participant0,
  const std::string&     participant1,
  const std::string&     localParticipant,
  com::PtrCommunication& communication ) const
{
  preciceTrace3 ( "connect()", participant0, participant1, localParticipant );
  assertion ( communication.use_count() > 0 );
  assertion ( not communication->isConnected() );
  utils::Parallel::initialize ( NULL, NULL, localParticipant );
  if ( participant0 == localParticipant ) {
    communication->requestConnection ( participant1, participant0, 0, 1 );
  }
  else {
    assertion ( participant1 == localParticipant );
    communication->acceptConnection ( participant1, participant0, 0, 1 );
  }
}

#endif // not PRECICE_NO_MPI

}}} // namespace precice, cplscheme, tests

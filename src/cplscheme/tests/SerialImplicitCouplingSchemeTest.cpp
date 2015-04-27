// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "SerialImplicitCouplingSchemeTest.hpp"
#include "cplscheme/SerialCouplingScheme.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"
#include "cplscheme/config/PostProcessingConfiguration.hpp"
#include "cplscheme/impl/ConvergenceMeasure.hpp"
#include "cplscheme/impl/AbsoluteConvergenceMeasure.hpp"
#include "cplscheme/impl/MinIterationConvergenceMeasure.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "cplscheme/Constants.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "geometry/config/GeometryConfiguration.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "m2n/M2N.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"
#include "utils/xml/XMLTag.hpp"
#include "utils/Dimensions.hpp"
#include "tarch/la/Vector.h"
#include "tarch/la/WrappedVector.h"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::cplscheme::tests::SerialImplicitCouplingSchemeTest)

namespace precice {
namespace cplscheme {
namespace tests {

using utils::Vector3D;

tarch::logging::Log SerialImplicitCouplingSchemeTest::
  _log ( "precice::cplscheme::tests::SerialImplicitCouplingSchemeTest" );


SerialImplicitCouplingSchemeTest:: SerialImplicitCouplingSchemeTest ()
:
  TestCase ( "precice::cplscheme::tests::SerialImplicitCouplingSchemeTest" ),
  _pathToTests (),
  MY_WRITE_CHECKPOINT ( constants::actionWriteIterationCheckpoint() ),
  MY_READ_CHECKPOINT ( constants::actionReadIterationCheckpoint() )
{}

void SerialImplicitCouplingSchemeTest:: setUp ()
{
  _pathToTests = utils::Globals::getPathToSources() + "/cplscheme/tests/";
}

void SerialImplicitCouplingSchemeTest:: run ()
{
  preciceTrace("run()");
# ifndef PRECICE_NO_MPI
  PRECICE_MASTER_ONLY {
    testMethod(testParseConfigurationWithRelaxation);
    testMethod(testExtrapolateData);
  }
  typedef utils::Parallel Par;
  preciceDebug("CommunicatorSize: " << Par::getCommunicatorSize());
  if (Par::getCommunicatorSize() > 1){
    // Do only use process 0 and 1 for the following tests
    std::vector<int> ranks;
    ranks += 0, 1;
    MPI_Comm comm = Par::getRestrictedCommunicator(ranks);
    if (Par::getProcessRank() <= 1){
      Par::setGlobalCommunicator(comm);
      testMethod(testAbsConvergenceMeasureSynchronized);
      testMethod(testConfiguredAbsConvergenceMeasureSynchronized);
      testMethod(testMinIterConvergenceMeasureSynchronized);
      testMethod(testMinIterConvergenceMeasureSynchronizedWithSubcycling);
      testMethod(testInitializeData);
      Par::setGlobalCommunicator(Par::getCommunicatorWorld());
    }
  }
# endif // not PRECICE_NO_MPI
}

#ifndef PRECICE_NO_MPI

void SerialImplicitCouplingSchemeTest:: testParseConfigurationWithRelaxation()
{
  preciceTrace("testParseConfigurationWithRelaxation()");
  using namespace mesh;

  std::string path(_pathToTests + "serial-implicit-cplscheme-relax-const-config.xml");

  utils::XMLTag root = utils::getRootTag();
  PtrDataConfiguration dataConfig(new DataConfiguration(root));
  dataConfig->setDimensions(3);
  PtrMeshConfiguration meshConfig(new MeshConfiguration(root, dataConfig));
  meshConfig->setDimensions(3);
  m2n::M2NConfiguration::SharedPointer m2nConfig(
      new m2n::M2NConfiguration(root));
  CouplingSchemeConfiguration cplSchemeConfig(root, meshConfig, m2nConfig);

  utils::configure(root, path);
  validate(cplSchemeConfig._postProcConfig->getPostProcessing().get() != NULL);
  meshConfig->setMeshSubIDs();
}

void SerialImplicitCouplingSchemeTest:: testExtrapolateData()
{
  preciceTrace("testExtrapolateData()");
  using namespace mesh;

  PtrMesh mesh(new Mesh("MyMesh", 3, false));
  PtrData data = mesh->createData("MyData", 1);
  int dataID = data->getID();
  mesh->createVertex(Vector3D(0.0));
  mesh->allocateDataValues();
  validateEquals(data->values().size(), 1);

  double maxTime = CouplingScheme::UNDEFINED_TIME;
  int maxTimesteps = 1;
  double dt = 1.0;
  std::string first = "First";
  std::string second = "Second";
  std::string accessor = second;
  com::Communication::SharedPointer com(new com::MPIDirectCommunication());
  m2n::M2N::SharedPointer globalCom(new m2n::M2N(com, m2n::DistributedComFactory::SharedPointer()));
  int maxIterations = 1;

  // Test first order extrapolation
  SerialCouplingScheme scheme(maxTime, maxTimesteps, dt, 16, first, second,
                              accessor, globalCom, constants::FIXED_DT,
                              BaseCouplingScheme::Implicit, maxIterations);

  scheme.addDataToSend(data, mesh, true);
  scheme.setExtrapolationOrder(1);
  scheme.setupDataMatrices(scheme.getSendData());
  CouplingData* cplData = scheme.getSendData(dataID);
  validate(cplData != NULL);
  validateEquals(cplData->values->size(), 1);
  validateEquals(cplData->oldValues.cols(), 2);
  validateEquals(cplData->oldValues.rows(), 1);
  validateNumericalEquals((*cplData->values)[0], 0.0);
  validateNumericalEquals(cplData->oldValues(0,0), 0.0);
  validateNumericalEquals(cplData->oldValues(0,1), 0.0);

  (*cplData->values)[0] = 1.0;
  scheme.setTimesteps(scheme.getTimesteps() + 1);
  scheme.extrapolateData(scheme.getSendData());
  validateNumericalEquals((*cplData->values)[0], 2.0);
  validateNumericalEquals(cplData->oldValues(0,0), 2.0);
  validateNumericalEquals(cplData->oldValues(0,1), 1.0);

  (*cplData->values)[0] = 4.0;
  scheme.setTimesteps(scheme.getTimesteps() + 1);
  scheme.extrapolateData(scheme.getSendData());
  validateNumericalEquals((*cplData->values)[0], 7.0);
  validateNumericalEquals(cplData->oldValues(0,0), 7.0);
  validateNumericalEquals(cplData->oldValues(0,1), 4.0);

  // Test second order extrapolation
  assign(*cplData->values) = 0.0;
  assign(cplData->oldValues) = 0.0;
  SerialCouplingScheme scheme2 ( maxTime, maxTimesteps, dt, 16, first, second,
                                 accessor, globalCom, constants::FIXED_DT,
                                 BaseCouplingScheme::Implicit, maxIterations);

  scheme2.addDataToSend ( data, mesh, false );
  scheme2.setExtrapolationOrder ( 2 );
  scheme2.setupDataMatrices (scheme2.getSendData());
  cplData = scheme2.getSendData ( dataID );
  validate ( cplData != NULL );
  validateEquals ( cplData->values->size(), 1 );
  validateEquals ( cplData->oldValues.cols(), 3 );
  validateEquals ( cplData->oldValues.rows(), 1 );
  validateNumericalEquals ( (*cplData->values)[0], 0.0 );
  validateNumericalEquals ( cplData->oldValues(0,0), 0.0 );
  validateNumericalEquals ( cplData->oldValues(0,1), 0.0 );
  validateNumericalEquals ( cplData->oldValues(0,2), 0.0 );

  (*cplData->values)[0] = 1.0;
  scheme2.setTimesteps ( scheme2.getTimesteps() + 1 );
  scheme2.extrapolateData (scheme2.getSendData());
  validateNumericalEquals ( (*cplData->values)[0], 2.0 );
  validateNumericalEquals ( cplData->oldValues(0,0), 2.0 );
  validateNumericalEquals ( cplData->oldValues(0,1), 1.0 );
  validateNumericalEquals ( cplData->oldValues(0,2), 0.0 );

  (*cplData->values)[0] = 4.0;
  scheme2.setTimesteps ( scheme2.getTimesteps() + 1 );
  scheme2.extrapolateData (scheme2.getSendData());
  validateNumericalEquals ( (*cplData->values)[0], 8.0 );
  validateNumericalEquals ( cplData->oldValues(0,0), 8.0 );
  validateNumericalEquals ( cplData->oldValues(0,1), 4.0 );
  validateNumericalEquals ( cplData->oldValues(0,2), 1.0 );
}

void SerialImplicitCouplingSchemeTest:: testAbsConvergenceMeasureSynchronized ()
{
   preciceTrace ( "testAbsConvergenceMeasureSynchronized()" );
   using namespace mesh;
   utils::Parallel::synchronizeProcesses ();
   assertion ( utils::Parallel::getCommunicatorSize() > 1 );

   utils::XMLTag root = utils::getRootTag();
   // Create a data configuration, to simplify configuration of data
   PtrDataConfiguration dataConfig ( new DataConfiguration(root) );
   dataConfig->setDimensions(3);
   dataConfig->addData ( "data0", 1 );
   dataConfig->addData ( "data1", 3 );

   MeshConfiguration meshConfig ( root, dataConfig );
   meshConfig.setDimensions(3);
   mesh::PtrMesh mesh ( new Mesh("mesh", 3, false) );
   mesh->createData ( "data0", 1 );
   mesh->createData ( "data1", 3 );
   mesh->createVertex ( Vector3D(0.0) );
   mesh->allocateDataValues ();
   meshConfig.addMesh ( mesh );

   // Create all parameters necessary to create an ImplicitCouplingScheme object
   com::Communication::SharedPointer communication ( new com::MPIDirectCommunication() );
   m2n::M2N::SharedPointer globalCom(new m2n::M2N(communication, m2n::DistributedComFactory::SharedPointer()));
   double maxTime = 1.0;
   int maxTimesteps = 3;
   double timestepLength = 0.1;
   std::string nameParticipant0 ( "participant0" );
   std::string nameParticipant1 ( "participant1" );
   std::string nameLocalParticipant ( "" );
   int sendDataIndex = -1;
   int receiveDataIndex = -1;
   if ( utils::Parallel::getProcessRank() == 0 ) {
      nameLocalParticipant = nameParticipant0;
      sendDataIndex = 0;
      receiveDataIndex = 1;
   }
   else if ( utils::Parallel::getProcessRank() == 1 ) {
      nameLocalParticipant = nameParticipant1;
      sendDataIndex = 1;
      receiveDataIndex = 0;
   }

   // Create the coupling scheme object
   cplscheme::SerialCouplingScheme cplScheme (
       maxTime, maxTimesteps, timestepLength, 16, nameParticipant0,
       nameParticipant1, nameLocalParticipant, globalCom, constants::FIXED_DT,
       BaseCouplingScheme::Implicit, 100);
   cplScheme.addDataToSend ( mesh->data()[sendDataIndex], mesh, false );
   cplScheme.addDataToReceive ( mesh->data()[receiveDataIndex], mesh, false );

   double convergenceLimit1 = sqrt(3.0); // when diff_vector = (1.0, 1.0, 1.0)
   impl::PtrConvergenceMeasure absoluteConvMeasure1 (
         new impl::AbsoluteConvergenceMeasure(convergenceLimit1) );
   cplScheme.addConvergenceMeasure (
         mesh->data()[1]->getID(), false, absoluteConvMeasure1 );

   // Expected iterations per implicit timesptep
   std::vector<int> validIterations;
   validIterations += 5, 5, 5;
   connect ( "participant0", "participant1", nameLocalParticipant, globalCom );
   runCoupling ( cplScheme, nameLocalParticipant, meshConfig, validIterations );
   globalCom->closeConnection();
}

//void SerialImplicitCouplingSchemeTest:: testAbsConvergenceMeasureAsync ()
//{
//   preciceTrace ( "testAbsConvergenceMeasureAsync()" );
//   utils::Parallel::synchronizeProcesses ();
//   assertion ( utils::Parallel::getCommunicatorSize() > 1 );
//
//   // Create a data configuration, to simplify configuration of data
//   mesh::PtrDataConfiguration dataConfig ( new mesh::DataConfiguration() );
//   dataConfig->addData ( "data0", mesh::Data::TYPE_DOUBLE );
//   dataConfig->addData ( "data1", mesh::Data::TYPE_VECTOR );
//
//   mesh::MeshConfiguration meshConfig ( dataConfig );
//   mesh::PtrMesh mesh ( new mesh::Mesh("mesh", false) );
//   mesh->setData ( dataConfig->data()[0] );
//   mesh->setData ( dataConfig->data()[1] );
//   mesh->createVertex ( Vector(0.0) );
//   meshConfig.addMesh ( mesh );
//
//   // Create all parameters necessary to create an ImplicitCouplingScheme object
//   com::Communication::SharedPointer communication ( new com::MPIDirectCommunication );
//   double maxTime = 1.0;
//   int maxTimesteps = 3;
//   double timestepLength = 0.1;
//   std::string nameParticipant0 ( "participant0" );
//   std::string nameParticipant1 ( "participant1" );
//   std::string nameLocalParticipant ( "" );
//   int sendDataIndex = -1;
//   int receiveDataIndex = -1;
//   std::vector<int> validIterations;
//   if ( utils::Parallel::getProcessRank() == 0 ) {
//      nameLocalParticipant = nameParticipant0;
//      sendDataIndex = 0;
//      receiveDataIndex = 1;
//      validIterations += 6, 6, 6;
//   }
//   else if ( utils::Parallel::getProcessRank() == 1 ) {
//      nameLocalParticipant = nameParticipant1;
//      sendDataIndex = 1;
//      receiveDataIndex = 0;
//      validIterations += 5, 6, 6;
//   }
//
//   // Create the coupling scheme object
//   cplscheme::SerialImplicitCouplingScheme cplScheme (
//      maxTime, maxTimesteps, timestepLength, nameParticipant0, nameParticipant1,
//      nameLocalParticipant, communication, 100, true );
//   cplScheme.addDataToSend ( mesh->data()[sendDataIndex], mesh );
//   cplScheme.addDataToReceive ( mesh->data()[receiveDataIndex], mesh );
//
//   // Add absolute convergence measures to the coupling scheme
//   double convergenceLimit0 = 0.5;
//   PtrConvergenceMeasure absoluteConvMeasure0 (
//      new AbsoluteConvergenceMeasure(convergenceLimit0) );
//   cplScheme.addConvergenceMeasure ( dataConfig->data()[0], absoluteConvMeasure0 );
//   double convergenceLimit1 = sqrt(3.0); // when diff_vector = (1.0, 1.0, 1.0)
//   PtrConvergenceMeasure absoluteConvMeasure1 (
//      new AbsoluteConvergenceMeasure(convergenceLimit1) );
//   cplScheme.addConvergenceMeasure ( dataConfig->data()[1], absoluteConvMeasure1 );
//
//   // Expected iterations per implicit timesptep
//   runCoupling ( cplScheme, nameLocalParticipant, meshConfig, validIterations );
//
//   // Split up processes to individual parts and run coupling scheme
////   if ( utils::Parallel::getProcessRank() == 0 ) {
////      runCoupling ( cplScheme, nameParticipant0, dataConfig, validIterations0 );
////   }
////   else if ( utils::Parallel::getProcessRank() == 1 ) {
////      runCoupling ( cplScheme, nameParticipant1, dataConfig, validIterations1 );
////   }
//}

void SerialImplicitCouplingSchemeTest:: testConfiguredAbsConvergenceMeasureSynchronized ()
{
   preciceTrace ( "testConfiguredAbsoluteConvergenceMeasureSync()" );
   using namespace mesh;
   utils::Parallel::synchronizeProcesses ();
   assertion ( utils::Parallel::getCommunicatorSize() > 1 );

   std::string configurationPath (
      _pathToTests + "serial-implicit-cplscheme-absolute-config.xml" );

   std::string nameLocalParticipant ( "" );
   if ( utils::Parallel::getProcessRank() == 0 ) {
      nameLocalParticipant = "participant0";
   }
   else if ( utils::Parallel::getProcessRank() == 1 ) {
      nameLocalParticipant = "participant1";
   }

   utils::XMLTag root = utils::getRootTag();
   PtrDataConfiguration dataConfig(new DataConfiguration(root));
   dataConfig->setDimensions(3);
   PtrMeshConfiguration meshConfig(new MeshConfiguration(root, dataConfig));
   meshConfig->setDimensions(3);
   m2n::M2NConfiguration::SharedPointer m2nConfig(new m2n::M2NConfiguration(root));
   geometry::GeometryConfiguration geoConfig(root, meshConfig);
   geoConfig.setDimensions(3);
   CouplingSchemeConfiguration cplSchemeConfig(root, meshConfig, m2nConfig);

   utils::configure(root, configurationPath);
   //validate(success);
   //validate(dataConfig->isValid());
   //validate(meshConfig->isValid());
   //validate(comConfig->isValid());
   //validate(geoConfig.isValid());
   //validate(cplSchemeConfig.isValid());
   meshConfig->setMeshSubIDs();
   m2n::M2N::SharedPointer m2n = m2nConfig->getM2N("participant0", "participant1");

   geoConfig.geometries()[0]->create ( *meshConfig->meshes()[0] );

   std::vector<int> validIterations;
   validIterations += 5, 5, 5;
   connect ( "participant0", "participant1", nameLocalParticipant, m2n );
   runCoupling ( *cplSchemeConfig.getCouplingScheme(nameLocalParticipant),
                 nameLocalParticipant, *meshConfig, validIterations );
   m2n->closeConnection();
}

void SerialImplicitCouplingSchemeTest:: testMinIterConvergenceMeasureSynchronized ()
{
   preciceTrace ( "testMinIterConvergenceMeasureSynchronized()" );
   utils::Parallel::synchronizeProcesses ();

   utils::XMLTag root = utils::getRootTag();
   // Create a data configuration, to simplify configuration of data
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

   // Create all parameters necessary to create an ImplicitCouplingScheme object
   com::Communication::SharedPointer communication ( new com::MPIDirectCommunication );
   m2n::M2N::SharedPointer globalCom ( new m2n::M2N(communication, m2n::DistributedComFactory::SharedPointer()) );
   double maxTime = 1.0;
   int maxTimesteps = 3;
   double timestepLength = 0.1;
   std::string nameParticipant0 ( "participant0" );
   std::string nameParticipant1 ( "participant1" );
   std::string nameLocalParticipant ( "" );
   int sendDataIndex = -1;
   int receiveDataIndex = -1;
   if ( utils::Parallel::getProcessRank() == 0 ) {
      nameLocalParticipant = nameParticipant0;
      sendDataIndex = 0;
      receiveDataIndex = 1;
   }
   else if ( utils::Parallel::getProcessRank() == 1 ) {
      nameLocalParticipant = nameParticipant1;
      sendDataIndex = 1;
      receiveDataIndex = 0;
   }

   // Create the coupling scheme object
   cplscheme::SerialCouplingScheme cplScheme (
     maxTime, maxTimesteps, timestepLength, 16, nameParticipant0, nameParticipant1,
     nameLocalParticipant, globalCom, constants::FIXED_DT,
     BaseCouplingScheme::Implicit, 100);
   cplScheme.addDataToSend ( mesh->data()[sendDataIndex], mesh, false );
   cplScheme.addDataToReceive ( mesh->data()[receiveDataIndex], mesh, false );

   // Add convergence measures
   int minIterations = 3;
   impl::PtrConvergenceMeasure minIterationConvMeasure1 (
      new impl::MinIterationConvergenceMeasure(minIterations) );
   cplScheme.addConvergenceMeasure (
      mesh->data()[1]->getID(), false, minIterationConvMeasure1 );

   // Expected iterations per implicit timesptep
   std::vector<int> validIterations;
   validIterations += 3, 3, 3;
   connect ( "participant0", "participant1", nameLocalParticipant, globalCom );
   runCoupling ( cplScheme, nameLocalParticipant, meshConfig, validIterations );
   globalCom->closeConnection();
}

//void SerialImplicitCouplingSchemeTest:: testMinIterConvergenceMeasureAsync ()
//{
//   utils::Parallel::synchronizeProcesses ();
//   preciceDebug ( "testMinIterConvergenceMeasure()", "Entering" );
//
//   // Create a data configuration, to simplify configuration of data
//   mesh::PtrDataConfiguration dataConfig ( new mesh::DataConfiguration() );
//   dataConfig->addData ( "data0", mesh::Data::TYPE_DOUBLE );
//   dataConfig->addData ( "data1", mesh::Data::TYPE_VECTOR );
//
//   mesh::MeshConfiguration meshConfig ( dataConfig );
//   mesh::PtrMesh mesh ( new mesh::Mesh("mesh", false) );
//   mesh->setData ( dataConfig->data()[0] );
//   mesh->setData ( dataConfig->data()[1] );
//   mesh->createVertex ( Vector(0.0) );
//   meshConfig.addMesh ( mesh );
//
//   // Create all parameters necessary to create an ImplicitCouplingScheme object
//   com::Communication::SharedPointer communication ( new com::MPIDirectCommunication );
//   double maxTime = 1.0;
//   int maxTimesteps = 3;
//   double timestepLength = 0.1;
//   std::string nameParticipant0 ( "participant0" );
//   std::string nameParticipant1 ( "participant1" );
//   std::string nameLocalParticipant ( "" );
//   int sendDataIndex = -1;
//   int receiveDataIndex = -1;
//   std::vector<int> validIterations;
//   if ( utils::Parallel::getProcessRank() == 0 ) {
//      nameLocalParticipant = nameParticipant0;
//      sendDataIndex = 0;
//      receiveDataIndex = 1;
//      validIterations += 3, 3, 3;
//   }
//   else if ( utils::Parallel::getProcessRank() == 1 ) {
//      nameLocalParticipant = nameParticipant1;
//      sendDataIndex = 1;
//      receiveDataIndex = 0;
//      validIterations += 2, 3, 3;
//   }
//
//   // Create the coupling scheme object
//   cplscheme::SerialImplicitCouplingScheme cplScheme (
//      maxTime, maxTimesteps, timestepLength, nameParticipant0, nameParticipant1,
//      nameLocalParticipant, communication, 100, true );
//   cplScheme.addDataToSend ( mesh->data()[sendDataIndex], mesh );
//   cplScheme.addDataToReceive ( mesh->data()[receiveDataIndex], mesh );
//
//   // Add data to be exchanged
////   cplScheme.addDataToExchange ( nameParticipant0,
////                                 dataConfig->data()[0].getID() );
////   cplScheme.addDataToExchange ( nameParticipant1,
////                                 dataConfig->data()[1].getID() );
//
//   // Add convergence measures
//   int minIterations = 3;
//   PtrConvergenceMeasure minIterationConvMeasure0 (
//      new MinIterationConvergenceMeasure(minIterations) );
//   cplScheme.addConvergenceMeasure ( dataConfig->data()[0], minIterationConvMeasure0 );
//
//   runCoupling ( cplScheme, nameLocalParticipant, meshConfig, validIterations );
//
////   PtrConvergenceMeasure minIterationConvMeasure1 (
////      new AbsoluteConvergenceMeasure<Vector>(minIterations) );
////   cplScheme.addConvergenceMeasure ( minIterationConvMeasure1,
////                                     dataConfig->data()[1].getID() );
//
//   // Expected iterations per implicit timesptep
////   std::vector<int> validIterations;
//
//   // Split up processes to individual parts and run coupling scheme
////   if ( utils::Parallel::getProcessRank() == 0 ) {
////      validIterations += 3, 3, 3;
////      runCoupling ( cplScheme, nameParticipant0, dataConfig, validIterations );
////   }
////   else if ( utils::Parallel::getProcessRank() == 1 ) {
////      validIterations += 2, 3, 3;
////      runCoupling ( cplScheme, nameParticipant1, dataConfig, validIterations );
////   }
//
//   preciceDebug ( "testMinIterConvergenceMeasure()", "Leaving" );
//}

void SerialImplicitCouplingSchemeTest:: runCoupling
(
  CouplingScheme&                cplScheme,
  const std::string&             nameParticipant,
  const mesh::MeshConfiguration& meshConfig,
  const std::vector<int>&        validIterations )
{
  preciceTrace1 ( "runCoupling", nameParticipant );
  validateEquals ( meshConfig.meshes().size(), 1 );
  mesh::PtrMesh mesh = meshConfig.meshes()[0];
  validateEquals ( mesh->data().size(), 2 );
  validate ( mesh->vertices().size() > 0 );
  mesh::Vertex& vertex = mesh->vertices()[0];
  int index = vertex.getID();
  utils::DynVector& dataValues0 = mesh->data()[0]->values();
  utils::DynVector& dataValues1 = mesh->data()[1]->values();
  double initialStepsizeData0 = 5.0;
  double stepsizeData0 = 5.0;
  Vector3D initialStepsizeData1(5.0);
  Vector3D stepsizeData1(5.0);
  double computedTime = 0.0;
  int computedTimesteps = 0;
  std::string nameParticipant0 ( "participant0" );
  std::string nameParticipant1 ( "participant1" );
  assertion ( (nameParticipant == nameParticipant0)
              || (nameParticipant == nameParticipant1) );
  int iterationCount = 0;
  std::vector<int>::const_iterator iterValidIterations = validIterations.begin();

  if ( nameParticipant == nameParticipant0 ) {
    cplScheme.initialize ( 0.0, 1 );
    validate ( not cplScheme.isCouplingTimestepComplete() );
    validate ( cplScheme.isActionRequired(MY_WRITE_CHECKPOINT) );
    validate ( not cplScheme.isActionRequired(MY_READ_CHECKPOINT) );
    validate ( not cplScheme.hasDataBeenExchanged() );

    // Tells coupling scheme, that a checkpoint has been created.
    // All required actions have to be performed before calling advance().
    cplScheme.performedAction ( MY_WRITE_CHECKPOINT );
    validate ( not cplScheme.isActionRequired(MY_WRITE_CHECKPOINT) );

    while ( cplScheme.isCouplingOngoing() ) {
      dataValues0[index] += stepsizeData0;
      preciceDebug ( "Wrote data with stepsize " << stepsizeData0 );
      // The max timestep length is required to be obeyed.
      double maxLengthTimestep = cplScheme.getNextTimestepMaxLength();
      cplScheme.addComputedTime ( maxLengthTimestep );
      cplScheme.advance();
      iterationCount++;
      preciceDebug ( "increased iterations to " << iterationCount
                     << ", limit < " << *iterValidIterations );
      // A coupling timestep is complete, when the coupling iterations are
      // globally converged and if subcycling steps have filled one global
      // timestep.
      if ( cplScheme.isCouplingTimestepComplete() ) {
        preciceDebug ( "timestep complete" );
        // Advance participant time and timestep
        computedTime += maxLengthTimestep;
        computedTimesteps ++;
        validateNumericalEquals ( computedTime, cplScheme.getTime() );
        validateEquals ( computedTimesteps, cplScheme.getTimesteps()-1 );
        // The iteration number is enforced by the controlled decrease of the
        // change of data written
        validateEquals ( iterationCount, *iterValidIterations );
        if ( cplScheme.isCouplingOngoing() ) {
          validate ( cplScheme.isActionRequired(MY_WRITE_CHECKPOINT) );
          cplScheme.performedAction(MY_WRITE_CHECKPOINT);
          validate ( not cplScheme.isActionRequired(MY_WRITE_CHECKPOINT) );
        }
        else {
          validate ( not cplScheme.isActionRequired(MY_WRITE_CHECKPOINT) );
          validate ( not cplScheme.isActionRequired(MY_READ_CHECKPOINT ) );
        }
        iterationCount = 0;
        iterValidIterations++;
        if ( iterValidIterations == validIterations.end() ) {
          validate ( not cplScheme.isCouplingOngoing() );
        }
        // Reset data values, to simulate same convergence behavior of
        // interface values in next timestep.
        stepsizeData0 = initialStepsizeData0;
      }
      else { // coupling timestep is not yet complete
        preciceDebug ( "timestep not complete" );
        validate ( cplScheme.isCouplingOngoing() );
        validate ( iterationCount < *iterValidIterations );
        validate ( cplScheme.isActionRequired(MY_READ_CHECKPOINT) );
        cplScheme.performedAction ( MY_READ_CHECKPOINT );
        validate ( not cplScheme.isActionRequired(MY_READ_CHECKPOINT) );
        // The written data value is decreased in a regular manner, in order
        // to achieve a predictable convergence.
        stepsizeData0 -= 1.0;
      }
      // In every coupling cycle, data is sent
      validate ( cplScheme.hasDataBeenExchanged() );
    }
    cplScheme.finalize (); // Ends the coupling scheme
    validateNumericalEquals ( computedTime, 0.3 );
    validateEquals ( computedTimesteps, 3 );
  }

  else if ( nameParticipant == nameParticipant1 ) {
    cplScheme.initialize ( 0.0, 1 );
    validate ( not cplScheme.isCouplingTimestepComplete() );
    validate ( cplScheme.isActionRequired(MY_WRITE_CHECKPOINT) );
    validate ( not cplScheme.isActionRequired(MY_READ_CHECKPOINT) );
    validate ( cplScheme.hasDataBeenExchanged() );

    // Tells coupling scheme, that a checkpoint has been created.
    // All required actions have to be performed before calling advance().
    cplScheme.performedAction ( MY_WRITE_CHECKPOINT );
    validate ( not cplScheme.isActionRequired(MY_WRITE_CHECKPOINT) );

    while ( cplScheme.isCouplingOngoing() ) {
      Vector3D currentData;
      assign(currentData) = tarch::la::slice<3>(dataValues1, index *3);
      currentData += stepsizeData1;
      tarch::la::slice<3>(dataValues1,index*3) = currentData;
      preciceDebug ( "Wrote data with stepsize " << stepsizeData1 );
      // The max timestep length is required to be obeyed.
      double maxLengthTimestep = cplScheme.getNextTimestepMaxLength();
      cplScheme.addComputedTime ( maxLengthTimestep );
      cplScheme.advance();
      iterationCount++;
      preciceDebug ( "increased iterations to " << iterationCount
        << ", limit < " << *iterValidIterations );
      // A coupling timestep is complete, when the coupling iterations are
      // globally converged and if subcycling steps have filled one global
      // timestep.
      if ( cplScheme.isCouplingTimestepComplete() ) {
        preciceDebug ( "timestep complete" );
        // Advance participant time and timestep
        computedTime += maxLengthTimestep;
        computedTimesteps ++;
        validateNumericalEquals ( computedTime, cplScheme.getTime() );
        validateEquals ( computedTimesteps, cplScheme.getTimesteps()-1 );
        // The iterations are enforced by the controlled decrease of the
        // change of data written
        validateEquals ( iterationCount, *iterValidIterations );
        if ( cplScheme.isCouplingOngoing() ) {
          validate ( cplScheme.isActionRequired(MY_WRITE_CHECKPOINT) );
          cplScheme.performedAction ( MY_WRITE_CHECKPOINT );
          validate ( not cplScheme.isActionRequired(MY_WRITE_CHECKPOINT) );
        }
        else {
          validate ( not cplScheme.isActionRequired(MY_WRITE_CHECKPOINT) );
          validate ( not cplScheme.isActionRequired(MY_READ_CHECKPOINT ) );
        }
        iterationCount = 0;
        iterValidIterations++;
        if ( iterValidIterations == validIterations.end() ) {
          validate ( ! cplScheme.isCouplingOngoing() );
        }
        // Reset data values, to simulate same convergence behavior of
        // interface values in next timestep.
        stepsizeData1 = initialStepsizeData1;
      }
      else { // coupling timestep is not yet complete
        preciceDebug ( "timestep not complete" );
        validate ( cplScheme.isCouplingOngoing() );
        validate ( iterationCount < *iterValidIterations );
        validate ( cplScheme.isActionRequired(MY_READ_CHECKPOINT) );
        // The load checkpoint action requires to fallback to the cplScheme of the
        // first implicit iteration of the current timestep/time.
        cplScheme.performedAction ( MY_READ_CHECKPOINT );
        validate ( not cplScheme.isActionRequired(MY_READ_CHECKPOINT) );
        // The written data value is decreased in a regular manner, in order
        // to achieve a predictable convergence.
        stepsizeData1 -= 1.0;
      }
      // In every coupling cycle, data is sent
      validate ( cplScheme.hasDataBeenExchanged() );
    }
    cplScheme.finalize (); // Ends the coupling scheme
    validateNumericalEquals ( computedTime, 0.3 );
    validateEquals ( computedTimesteps, 3 );
  }
}

void SerialImplicitCouplingSchemeTest::
     testMinIterConvergenceMeasureSynchronizedWithSubcycling ()
{
   preciceTrace ( "testMinIterConvergenceMeasureSynchronizedWithSubcycling()" );
   utils::Parallel::synchronizeProcesses ();

   utils::XMLTag root = utils::getRootTag();
   // Create a data configuration, to simplify configuration of data
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

   // Create all parameters necessary to create an ImplicitCouplingScheme object
   com::Communication::SharedPointer communication ( new com::MPIDirectCommunication );
   m2n::M2N::SharedPointer globalCom ( new m2n::M2N(communication, m2n::DistributedComFactory::SharedPointer()));
   double maxTime = 1.0;
   int maxTimesteps = 3;
   double timestepLength = 0.1;
   std::string nameParticipant0 ( "participant0" );
   std::string nameParticipant1 ( "participant1" );
   std::string nameLocalParticipant ( "" );
   int sendDataIndex = -1;
   int receiveDataIndex = -1;
   std::vector<int> validIterations;
   if ( utils::Parallel::getProcessRank() == 0 ) {
      nameLocalParticipant = nameParticipant0;
      sendDataIndex = 0;
      receiveDataIndex = 1;
      validIterations += 3, 3, 3;
   }
   else if ( utils::Parallel::getProcessRank() == 1 ) {
      nameLocalParticipant = nameParticipant1;
      sendDataIndex = 1;
      receiveDataIndex = 0;
      validIterations += 3, 3, 3;
   }

   // Create the coupling scheme object
   cplscheme::SerialCouplingScheme cplScheme (
      maxTime, maxTimesteps, timestepLength, 16, nameParticipant0, nameParticipant1,
      nameLocalParticipant, globalCom, constants::FIXED_DT,
      BaseCouplingScheme::Implicit, 100);
   cplScheme.addDataToSend ( mesh->data()[sendDataIndex], mesh, false );
   cplScheme.addDataToReceive ( mesh->data()[receiveDataIndex], mesh, false );

   // Add convergence measures
   int minIterations = 3;
   impl::PtrConvergenceMeasure minIterationConvMeasure1 (
         new impl::MinIterationConvergenceMeasure(minIterations) );
   cplScheme.addConvergenceMeasure (
         mesh->data()[1]->getID(), false, minIterationConvMeasure1 );
   connect ( "participant0", "participant1", nameLocalParticipant, globalCom );
   runCouplingWithSubcycling (
      cplScheme, nameLocalParticipant, meshConfig, validIterations );
   globalCom->closeConnection();
}

void SerialImplicitCouplingSchemeTest:: testInitializeData()
{
  preciceTrace("testInitializeData()");
  utils::Parallel::synchronizeProcesses();

  utils::XMLTag root = utils::getRootTag();

  // Create a data configuration, to simplify configuration of data

  mesh::PtrDataConfiguration dataConfig(new mesh::DataConfiguration(root));
  dataConfig->setDimensions(3);
  dataConfig->addData("Data0", 1);
  dataConfig->addData("Data1", 3);

  mesh::MeshConfiguration meshConfig(root, dataConfig);
  meshConfig.setDimensions(3);
  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 3, false));
  mesh->createData("Data0", 1);
  mesh->createData("Data1", 3);
  mesh->createVertex(Vector3D(0.0));
  mesh->allocateDataValues();
  meshConfig.addMesh(mesh);

  // Create all parameters necessary to create an ImplicitCouplingScheme object
  com::Communication::SharedPointer communication(new com::MPIDirectCommunication);
  m2n::M2N::SharedPointer globalCom ( new m2n::M2N(communication, m2n::DistributedComFactory::SharedPointer())  );
  double maxTime = 1.0;
  int maxTimesteps = 3;
  double timestepLength = 0.1;
  std::string nameParticipant0("participant0");
  std::string nameParticipant1("participant1");
  std::string nameLocalParticipant("");
  int sendDataIndex = -1;
  int receiveDataIndex = -1;
  bool initData = false;
  if (utils::Parallel::getProcessRank() == 0){
     nameLocalParticipant = nameParticipant0;
     sendDataIndex = 0;
     receiveDataIndex = 1;
  }
  else if (utils::Parallel::getProcessRank() == 1){
     nameLocalParticipant = nameParticipant1;
     sendDataIndex = 1;
     receiveDataIndex = 0;
     initData = true;
  }

  // Create the coupling scheme object
  cplscheme::SerialCouplingScheme cplScheme(
     maxTime, maxTimesteps, timestepLength, 16, nameParticipant0, nameParticipant1,
     nameLocalParticipant, globalCom, constants::FIXED_DT,
     BaseCouplingScheme::Implicit, 100);
  cplScheme.addDataToSend(mesh->data()[sendDataIndex], mesh, initData);
  cplScheme.addDataToReceive(mesh->data()[receiveDataIndex], mesh, not initData);

  // Add convergence measures
  int minIterations = 3;
  impl::PtrConvergenceMeasure minIterationConvMeasure1 (
        new impl::MinIterationConvergenceMeasure(minIterations) );
  cplScheme.addConvergenceMeasure (
        mesh->data()[1]->getID(), false, minIterationConvMeasure1 );
  connect(nameParticipant0, nameParticipant1, nameLocalParticipant, globalCom);

  std::string writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterationCheckpoint(constants::actionReadIterationCheckpoint());

  cplScheme.initialize(0.0, 1);

  if (nameLocalParticipant == nameParticipant0){
    cplScheme.initializeData();
    validate(cplScheme.hasDataBeenExchanged());
    utils::DynVector& values = mesh->data(1)->values();
    validateWithParams1(tarch::la::equals(values, Vector3D(1.0, 2.0, 3.0)), values);
    mesh->data(0)->values() = 4.0;
    while (cplScheme.isCouplingOngoing()){
      if (cplScheme.isActionRequired(writeIterationCheckpoint)){
        cplScheme.performedAction(writeIterationCheckpoint);
      }
      if (cplScheme.isActionRequired(readIterationCheckpoint)){
        cplScheme.performedAction(readIterationCheckpoint);
      }
      cplScheme.addComputedTime(timestepLength);
      cplScheme.advance();
    }
  }
  else {
    assertion(nameLocalParticipant == nameParticipant1);
    validate(cplScheme.isActionRequired(constants::actionWriteInitialData()));
    cplScheme.performedAction(constants::actionWriteInitialData());
    utils::DynVector& values = mesh->data(0)->values();
    validateWithParams1(tarch::la::equals(values(0), 0.0), values);
    mesh->data(1)->values() = Vector3D(1.0, 2.0, 3.0);
    cplScheme.initializeData();
    validate(cplScheme.hasDataBeenExchanged());
    validateWithParams1(tarch::la::equals(values(0), 4.0), values);
    while (cplScheme.isCouplingOngoing()){
      if (cplScheme.isActionRequired(writeIterationCheckpoint)){
        cplScheme.performedAction(writeIterationCheckpoint);
      }
      cplScheme.addComputedTime(timestepLength);
      cplScheme.advance();
      if (cplScheme.isActionRequired(readIterationCheckpoint)){
        cplScheme.performedAction(readIterationCheckpoint);
      }
    }
  }
  cplScheme.finalize();
}

void SerialImplicitCouplingSchemeTest:: runCouplingWithSubcycling
(
  CouplingScheme&                cplScheme,
  const std::string&             nameParticipant,
  const mesh::MeshConfiguration& meshConfig,
  const std::vector<int>&        validIterations )
{
  preciceTrace1 ( "runCouplingWithSubcycling()", nameParticipant );

  validateEquals ( meshConfig.meshes().size(), 1 );
  mesh::PtrMesh mesh = meshConfig.meshes()[0];
  validateEquals ( mesh->data().size(), 2 );
  validate ( mesh->vertices().size() > 0 );
  double initialStepsizeData0 = 5.0;
  double stepsizeData0 = 5.0;
  Vector3D initialStepsizeData1 ( 5.0 );
  Vector3D stepsizeData1 ( 5.0 );
  double computedTime = 0.0;
  int computedTimesteps = 0;
  std::string nameParticipant0 ( "participant0" );
  std::string nameParticipant1 ( "participant1" );
  assertion ( (nameParticipant == nameParticipant0)
    || (nameParticipant == nameParticipant1) );
  int iterationCount = 0;
  std::vector<int>::const_iterator iterValidIterations =
    validIterations.begin();

  if ( nameParticipant == nameParticipant0 ) {
    iterationCount++; // different handling due to subcycling
    cplScheme.initialize ( 0.0, 1 );
    validate ( not cplScheme.isCouplingTimestepComplete() );
    validate ( cplScheme.isActionRequired(MY_WRITE_CHECKPOINT) );
    validate ( not cplScheme.isActionRequired(MY_READ_CHECKPOINT) );
    validate ( not cplScheme.hasDataBeenExchanged() );

    // Tells coupling scheme, that a checkpoint has been created.
    // All required actions have to be performed before calling advance().
    cplScheme.performedAction ( MY_WRITE_CHECKPOINT );
    validate ( not cplScheme.isActionRequired(MY_WRITE_CHECKPOINT) );

    double maxTimestepLength = cplScheme.getNextTimestepMaxLength();
    double computedTimestepLength = maxTimestepLength / 2.0;
    int subcyclingStep = 0;

    // Main coupling loop
    while ( cplScheme.isCouplingOngoing() ){
      cplScheme.addComputedTime ( computedTimestepLength );
      cplScheme.advance();
      // A coupling timestep is complete, when the coupling iterations are
      // globally converged and if subcycling steps have filled one global
      // timestep.
      if ( cplScheme.isCouplingTimestepComplete() ){
        // Advance participant time and timestep
        computedTime += maxTimestepLength;
        computedTimesteps ++;
        validateNumericalEquals ( computedTime, cplScheme.getTime() );
        validateEquals ( computedTimesteps, cplScheme.getTimesteps()-1 );
        // The iteration number is enforced by the controlled decrease of the
        // change of data written
        validateEquals ( iterationCount, *iterValidIterations );
        if ( cplScheme.isCouplingOngoing() ) {
          validate ( cplScheme.isActionRequired(MY_WRITE_CHECKPOINT) );
          cplScheme.performedAction(MY_WRITE_CHECKPOINT);
          validate ( ! cplScheme.isActionRequired(MY_WRITE_CHECKPOINT) );
        }
        else {
          validate ( ! cplScheme.isActionRequired(MY_WRITE_CHECKPOINT) );
          validate ( ! cplScheme.isActionRequired(MY_READ_CHECKPOINT ) );
        }
        iterationCount = 1;
        iterValidIterations++;
        if ( iterValidIterations == validIterations.end() ) {
          validate ( ! cplScheme.isCouplingOngoing() );
        }
        // Reset data values, to simulate same convergence behavior of
        // interface values in next timestep.
        stepsizeData0 = initialStepsizeData0;
        validateEquals ( subcyclingStep, 1 );
        subcyclingStep = 0;
      }
      else { // coupling timestep is not yet complete
        validate ( cplScheme.isCouplingOngoing() );
        // If length of global timestep is reached
        if ( cplScheme.hasDataBeenExchanged() ) {
          validate ( iterationCount <= *iterValidIterations );
          validate ( cplScheme.isActionRequired(MY_READ_CHECKPOINT) );
          validate ( ! cplScheme.isActionRequired(MY_WRITE_CHECKPOINT) );
          cplScheme.performedAction ( MY_READ_CHECKPOINT );
          validate ( ! cplScheme.isActionRequired(MY_READ_CHECKPOINT) );
          // The written data value is decreased in a regular manner, in order
          // to achieve a predictable convergence.
          stepsizeData0 -= 1.0;
          subcyclingStep = 0; // Subcycling steps
          iterationCount++; // Implicit coupling iterations
          //precicePrint ( "increased iterations to " << iterationCount );
        }
        else { // If subcycling
          validate ( iterationCount <= *iterValidIterations );
          validate ( ! cplScheme.isActionRequired(MY_READ_CHECKPOINT) );
          validate ( ! cplScheme.isActionRequired(MY_WRITE_CHECKPOINT) );
          validate ( subcyclingStep < 2 );
          subcyclingStep++;
        }
      }
    }
    cplScheme.finalize (); // Ends the coupling scheme
    validateNumericalEquals ( computedTime, 0.3 );
    validateEquals ( computedTimesteps, 3 );
  }

  else if ( nameParticipant == nameParticipant1 ) {
    iterationCount++; // different handling due to subcycling
    cplScheme.initialize ( 0.0, 1 );
    validate ( not cplScheme.isCouplingTimestepComplete() );
    validate ( cplScheme.isActionRequired(MY_WRITE_CHECKPOINT) );
    validate ( not cplScheme.isActionRequired(MY_READ_CHECKPOINT) );
    validate ( cplScheme.hasDataBeenExchanged() );

    // Tells coupling scheme, that a checkpoint has been created.
    // All required actions have to be performed before calling advance().
    cplScheme.performedAction ( MY_WRITE_CHECKPOINT );
    validate ( ! cplScheme.isActionRequired(MY_WRITE_CHECKPOINT) );

    double maxTimestepLength = cplScheme.getNextTimestepMaxLength ();
    double preferredTimestepLength = maxTimestepLength / 2.5;
    double computedTimestepLength = preferredTimestepLength;
    int subcyclingStep = 0;

    // Main coupling loop
    while ( cplScheme.isCouplingOngoing() ){
      cplScheme.addComputedTime ( computedTimestepLength );
      cplScheme.advance();
      computedTimestepLength =
          cplScheme.getNextTimestepMaxLength() < preferredTimestepLength
          ? cplScheme.getNextTimestepMaxLength()
          : preferredTimestepLength;
      // A coupling timestep is complete, when the coupling iterations are
      // globally converged and if subcycling steps have filled one global
      // timestep.
      if ( cplScheme.isCouplingTimestepComplete() ){
        //            precicePrint ( "timestep complete" );
        // Advance participant time and timestep
        computedTime += maxTimestepLength;
        computedTimesteps ++;
        validateNumericalEquals ( computedTime, cplScheme.getTime() );
        validateEquals ( computedTimesteps, cplScheme.getTimesteps()-1 );
        // The iteration number is enforced by the controlled decrease of the
        // change of data written
        validateEquals ( iterationCount, *iterValidIterations );
        if ( cplScheme.isCouplingOngoing() ) {
          validate ( cplScheme.isActionRequired(MY_WRITE_CHECKPOINT) );
          cplScheme.performedAction(MY_WRITE_CHECKPOINT);
          validate ( ! cplScheme.isActionRequired(MY_WRITE_CHECKPOINT) );
        }
        else {
          validate ( ! cplScheme.isActionRequired(MY_WRITE_CHECKPOINT) );
          validate ( ! cplScheme.isActionRequired(MY_READ_CHECKPOINT ) );
        }
        iterationCount = 1;
        iterValidIterations++;
        if ( iterValidIterations == validIterations.end() ) {
          validate ( ! cplScheme.isCouplingOngoing() );
        }
        // Reset data values, to simulate same convergence behavior of
        // interface values in next timestep.
        stepsizeData1 = initialStepsizeData1;
        validateEquals ( subcyclingStep, 2 );
        subcyclingStep = 0;
      }
      else { // coupling timestep is not yet complete
        validate ( cplScheme.isCouplingOngoing() );
        // If length of global timestep is reached
        if ( cplScheme.hasDataBeenExchanged() ) {
          validate ( iterationCount <= *iterValidIterations );
          validate ( cplScheme.isActionRequired(MY_READ_CHECKPOINT) );
          validate ( ! cplScheme.isActionRequired(MY_WRITE_CHECKPOINT) );
          cplScheme.performedAction ( MY_READ_CHECKPOINT );
          validate ( ! cplScheme.isActionRequired(MY_READ_CHECKPOINT) );
          // The written data value is decreased in a regular manner, in order
          // to achieve a predictable convergence.
          stepsizeData1 -= 1.0;
          subcyclingStep = 0; // Subcycling steps
          iterationCount++; // Implicit coupling iterations
          //precicePrint ( "increased iterations to " << iterationCount );
        }
        else { // If subcycling
          validate ( iterationCount <= *iterValidIterations );
          validate ( ! cplScheme.isActionRequired(MY_READ_CHECKPOINT) );
          validate ( ! cplScheme.isActionRequired(MY_WRITE_CHECKPOINT) );
          validate ( subcyclingStep < 3 );
          subcyclingStep++;
        }
      }
    }
    cplScheme.finalize (); // Ends the coupling scheme
    validateNumericalEquals ( computedTime, 0.3 );
    validateEquals ( computedTimesteps, 3 );
  }
}

void SerialImplicitCouplingSchemeTest:: connect
(
  const std::string&      participant0,
  const std::string&      participant1,
  const std::string&      localParticipant,
  m2n::M2N::SharedPointer& communication ) const
{
  assertion ( communication.use_count() > 0 );
  assertion ( not communication->isConnected() );
  utils::Parallel::initialize ( NULL, NULL, localParticipant );
  if ( participant0 == localParticipant ) {
    communication->requestMasterConnection ( participant1, participant0);
  }
  else {
    assertion ( participant1 == localParticipant );
    communication->acceptMasterConnection ( participant1, participant0);
  }
}

#endif // not PRECICE_NO_MPI

}}} // namespace precice, cplscheme, tests

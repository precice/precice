// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ParallelImplicitCouplingSchemeTest.hpp"
#include "cplscheme/ParallelImplicitCouplingScheme.hpp"
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
#include "com/config/CommunicationConfiguration.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"
#include "utils/xml/XMLTag.hpp"
#include "utils/Dimensions.hpp"
#include "tarch/la/Vector.h"
#include "tarch/la/WrappedVector.h"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::cplscheme::tests::ParallelImplicitCouplingSchemeTest)

namespace precice {
namespace cplscheme {
namespace tests {

using utils::Vector3D;

tarch::logging::Log ParallelImplicitCouplingSchemeTest::
  _log ( "precice::cplscheme::tests::ParallelImplicitCouplingSchemeTest" );


ParallelImplicitCouplingSchemeTest:: ParallelImplicitCouplingSchemeTest ()
:
  TestCase ( "precice::cplscheme::tests::ParallelImplicitCouplingSchemeTest" ),
  _pathToTests (),
  MY_WRITE_CHECKPOINT ( constants::actionWriteIterationCheckpoint() ),
  MY_READ_CHECKPOINT ( constants::actionReadIterationCheckpoint() )
{}

void ParallelImplicitCouplingSchemeTest:: setUp ()
{
  _pathToTests = utils::Globals::getPathToSources() + "/cplscheme/tests/";
}

void ParallelImplicitCouplingSchemeTest:: run ()
{
# ifndef PRECICE_NO_MPI
  PRECICE_MASTER_ONLY {
    testMethod(testParseConfigurationWithRelaxation);
  }
  typedef utils::Parallel Par;
    if (Par::getCommunicatorSize() > 1){
      // Do only use process 0 and 1 for the following tests
      std::vector<int> ranks;
      ranks += 0, 1;
      MPI_Comm comm = Par::getRestrictedCommunicator(ranks);
      if (Par::getProcessRank() <= 1){
        Par::setGlobalCommunicator(comm);
        testMethod(testInitializeData);
        Par::setGlobalCommunicator(Par::getCommunicatorWorld());
      }
    }
# endif // not PRECICE_NO_MPI
}

#ifndef PRECICE_NO_MPI

void ParallelImplicitCouplingSchemeTest:: testParseConfigurationWithRelaxation()
{
  preciceTrace("testParseConfigurationWithRelaxation()");
  using namespace mesh;

  std::string path(_pathToTests + "parallel-implicit-cplscheme-relax-const-config.xml");

  utils::XMLTag root = utils::getRootTag();
  PtrDataConfiguration dataConfig(new DataConfiguration(root));
  dataConfig->setDimensions(3);
  PtrMeshConfiguration meshConfig(new MeshConfiguration(root, dataConfig));
  meshConfig->setDimensions(3);
  com::PtrCommunicationConfiguration comConfig(
      new com::CommunicationConfiguration(root));
  CouplingSchemeConfiguration cplSchemeConfig(root, meshConfig, comConfig);

  utils::configure(root, path);
  validate(cplSchemeConfig._postProcConfig->getPostProcessing().get() != NULL);
  meshConfig->setMeshSubIDs();
}

void ParallelImplicitCouplingSchemeTest:: testInitializeData()
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

  // Create all parameters necessary to create a ParallelImplicitCouplingScheme object
  com::PtrCommunication communication(new com::MPIDirectCommunication);
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
     initData = true;
  }
  else if (utils::Parallel::getProcessRank() == 1){
     nameLocalParticipant = nameParticipant1;
     sendDataIndex = 1;
     receiveDataIndex = 0;
     initData = true;
  }

  // Create the coupling scheme object
  cplscheme::ParallelImplicitCouplingScheme cplScheme(
     maxTime, maxTimesteps, timestepLength, 16, nameParticipant0, nameParticipant1,
     nameLocalParticipant, communication, 100, constants::FIXED_DT);
  cplScheme.addDataToSend(mesh->data()[sendDataIndex], initData);
  cplScheme.addDataToReceive(mesh->data()[receiveDataIndex], initData);



  // Add convergence measures
  int minIterations = 3;
  impl::PtrConvergenceMeasure minIterationConvMeasure1 (
        new impl::MinIterationConvergenceMeasure(minIterations) );
  impl::PtrConvergenceMeasure minIterationConvMeasure2 (
          new impl::MinIterationConvergenceMeasure(minIterations) );
  cplScheme.addConvergenceMeasure (
        mesh->data()[1]->getID(), false, minIterationConvMeasure1 );
  cplScheme.addConvergenceMeasure (
          mesh->data()[0]->getID(), false, minIterationConvMeasure2 );
  connect(nameParticipant0, nameParticipant1, nameLocalParticipant, communication);

  std::string writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterationCheckpoint(constants::actionReadIterationCheckpoint());

  cplScheme.initialize(0.0, 0);




  if (nameLocalParticipant == nameParticipant0){
    validate(cplScheme.isActionRequired(constants::actionWriteInitialData()));
    mesh->data(0)->values() = 4.0;
    cplScheme.performedAction(constants::actionWriteInitialData());
    cplScheme.initializeData();
    validate(cplScheme.hasDataBeenExchanged());
    utils::DynVector& values = mesh->data(1)->values();
    validateWithParams1(tarch::la::equals(values, Vector3D(1.0, 2.0, 3.0)), values);

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
    utils::DynVector& values = mesh->data(0)->values();
    validate(cplScheme.isActionRequired(constants::actionWriteInitialData()));
    mesh->data(1)->values() = Vector3D(1.0, 2.0, 3.0);
    cplScheme.performedAction(constants::actionWriteInitialData());
    validateWithParams1(tarch::la::equals(values(0), 0.0), values);
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

void ParallelImplicitCouplingSchemeTest:: connect
(
  const std::string&      participant0,
  const std::string&      participant1,
  const std::string&      localParticipant,
  com::PtrCommunication& communication ) const
{
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

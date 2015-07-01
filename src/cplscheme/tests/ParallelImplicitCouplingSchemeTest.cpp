// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ParallelImplicitCouplingSchemeTest.hpp"
#include "cplscheme/ParallelCouplingScheme.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"
#include "cplscheme/config/PostProcessingConfiguration.hpp"
#include "cplscheme/impl/ConvergenceMeasure.hpp"
#include "cplscheme/impl/AbsoluteConvergenceMeasure.hpp"
#include "cplscheme/impl/MinIterationConvergenceMeasure.hpp"
#include "cplscheme/impl/IQNILSPostProcessing.hpp"
#include "cplscheme/impl/MVQNPostProcessing.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "cplscheme/Constants.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "geometry/config/GeometryConfiguration.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "m2n/GatherScatterCommunication.hpp"
#include "m2n/config/M2NConfiguration.hpp"
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
    testMethod(testVIQNPP);
    testMethod(testMVQNPP);
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
  m2n::M2NConfiguration::SharedPointer m2nConfig(
      new m2n::M2NConfiguration(root));
  CouplingSchemeConfiguration cplSchemeConfig(root, meshConfig, m2nConfig);
  
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
  com::Communication::SharedPointer communication(new com::MPIDirectCommunication);
  m2n::M2N::SharedPointer globalCom(new m2n::M2N(communication, m2n::DistributedComFactory::SharedPointer()));
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
  cplscheme::ParallelCouplingScheme cplScheme(
    maxTime, maxTimesteps, timestepLength, 16, nameParticipant0, nameParticipant1,
    nameLocalParticipant, globalCom, constants::FIXED_DT, BaseCouplingScheme::Implicit, 100);
  cplScheme.addDataToSend(mesh->data()[sendDataIndex], mesh, initData);
  cplScheme.addDataToReceive(mesh->data()[receiveDataIndex], mesh, initData);

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
  connect(nameParticipant0, nameParticipant1, nameLocalParticipant, globalCom);

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

void ParallelImplicitCouplingSchemeTest:: testVIQNPP()
{
  preciceTrace("testVIQNPP()");

  //use two vectors and see if underrelaxation works

  double initialRelaxation = 0.01;
  int    maxIterationsUsed = 50;
  int    timestepsReused = 6;
  double singularityLimit = 1e-10;
  std::vector<int> dataIDs;
  dataIDs.push_back(0);
  dataIDs.push_back(1);
  std::map<int, double> scalings;
  scalings.insert(std::make_pair(0,1.0));
  scalings.insert(std::make_pair(1,1.0));
  mesh::PtrMesh dummyMesh ( new mesh::Mesh("dummyMesh", 3, false) );

  cplscheme::impl::IQNILSPostProcessing pp(initialRelaxation,maxIterationsUsed,
                                           timestepsReused, singularityLimit, dataIDs, scalings);

  //init displacements
  utils::DynVector dvalues;
  dvalues.append(1.0);
  dvalues.append(2.0);
  dvalues.append(3.0);
  dvalues.append(4.0);

  utils::DynVector dcol1;
  dcol1.append(1.0);
  dcol1.append(1.0);
  dcol1.append(1.0);
  dcol1.append(1.0);

  PtrCouplingData dpcd(new CouplingData(&dvalues,dummyMesh,false,1));

  //init forces
  utils::DynVector fvalues;
  fvalues.append(0.1);
  fvalues.append(0.1);
  fvalues.append(0.1);

  utils::DynVector fcol1;
  fcol1.append(0.2);
  fcol1.append(0.2);
  fcol1.append(0.2);

  PtrCouplingData fpcd(new CouplingData(&fvalues,dummyMesh,false,1));

  DataMap data;
  data.insert(std::pair<int,PtrCouplingData>(0,dpcd));
  data.insert(std::pair<int,PtrCouplingData>(1,fpcd));

//  foreach (DataMap::value_type& pair, data){
//    std::cout << *pair.second->values << "\n";
//    std::cout << pair.second->oldValues << "\n";
//  }

  pp.initialize(data);

  dpcd->oldValues.column(0) = dcol1;
  fpcd->oldValues.column(0) = fcol1;

  pp.performPostProcessing(data);

  validateWithParams1(tarch::la::equals((*data.at(0)->values)(0), 1.00), (*data.at(0)->values)(0));
  validateWithParams1(tarch::la::equals((*data.at(0)->values)(1), 1.01), (*data.at(0)->values)(1));
  validateWithParams1(tarch::la::equals((*data.at(0)->values)(2), 1.02), (*data.at(0)->values)(2));
  validateWithParams1(tarch::la::equals((*data.at(0)->values)(3), 1.03), (*data.at(0)->values)(3));
  validateWithParams1(tarch::la::equals((*data.at(1)->values)(0), 0.199), (*data.at(1)->values)(0));
  validateWithParams1(tarch::la::equals((*data.at(1)->values)(1), 0.199), (*data.at(1)->values)(1));
  validateWithParams1(tarch::la::equals((*data.at(1)->values)(2), 0.199), (*data.at(1)->values)(2));

  utils::DynVector newdvalues;
  newdvalues.append(10.0);
  newdvalues.append(10.0);
  newdvalues.append(10.0);
  newdvalues.append(10.0);
  data.begin()->second->values = &newdvalues;

  pp.performPostProcessing(data);

  validateWithParams1(tarch::la::equals((*data.at(0)->values)(0), -5.63855295490201413600e-01), (*data.at(0)->values)(0));
  validateWithParams1(tarch::la::equals((*data.at(0)->values)(1), 6.09906404008709657205e-01), (*data.at(0)->values)(1));
  validateWithParams1(tarch::la::equals((*data.at(0)->values)(2), 1.78366810350762072801e+0), (*data.at(0)->values)(2));
  validateWithParams1(tarch::la::equals((*data.at(0)->values)(3), 2.95742980300653179881e+00), (*data.at(0)->values)(3));
  validateWithParams1(tarch::la::equals((*data.at(1)->values)(0), 8.27975917496077823410e-02), (*data.at(1)->values)(0));
  validateWithParams1(tarch::la::equals((*data.at(1)->values)(1), 8.27975917496077823410e-02), (*data.at(1)->values)(1));
  validateWithParams1(tarch::la::equals((*data.at(1)->values)(2), 8.27975917496077823410e-02), (*data.at(1)->values)(2));
}


void ParallelImplicitCouplingSchemeTest:: testMVQNPP()
{
  preciceTrace("testMVQNPP()");
  
  //use two vectors and see if underrelaxation works
  
  double initialRelaxation = 0.01;
  int    maxIterationsUsed = 50;
  int    timestepsReused = 6;
  double singularityLimit = 1e-10;
  std::vector<int> dataIDs;
  dataIDs.push_back(0);
  dataIDs.push_back(1);
  std::map<int, double> scalings;
  scalings.insert(std::make_pair(0,1.0));
  scalings.insert(std::make_pair(1,1.0));
  mesh::PtrMesh dummyMesh ( new mesh::Mesh("dummyMesh", 3, false) );

  
  cplscheme::impl::MVQNPostProcessing pp(initialRelaxation,maxIterationsUsed,
                                         timestepsReused, singularityLimit, dataIDs, scalings);
  
  //init displacements
  utils::DynVector dvalues;
  dvalues.append(1.0);
  dvalues.append(2.0);
  dvalues.append(3.0);
  dvalues.append(4.0);
  
  utils::DynVector dcol1;
  dcol1.append(1.0);
  dcol1.append(1.0);
  dcol1.append(1.0);
  dcol1.append(1.0);
  
  PtrCouplingData dpcd(new CouplingData(&dvalues,dummyMesh,false,1));
  
  //init forces
  utils::DynVector fvalues;
  fvalues.append(0.1);
  fvalues.append(0.1);
  fvalues.append(0.1);
  
  utils::DynVector fcol1;
  fcol1.append(0.2);
  fcol1.append(0.2);
  fcol1.append(0.2);
  
  PtrCouplingData fpcd(new CouplingData(&fvalues,dummyMesh,false,1));
  
  DataMap data;
  data.insert(std::pair<int,PtrCouplingData>(0,dpcd));
  data.insert(std::pair<int,PtrCouplingData>(1,fpcd));
  
//  foreach (DataMap::value_type& pair, data){
//    std::cout << *pair.second->values << "\n";
//    std::cout << pair.second->oldValues << "\n";
//  }
  
  pp.initialize(data);
  
  dpcd->oldValues.column(0) = dcol1;
  fpcd->oldValues.column(0) = fcol1;
  
  pp.performPostProcessing(data);
  
  validateWithParams1(tarch::la::equals((*data.at(0)->values)(0), 1.00000000000000000000), (*data.at(0)->values)(0));
  validateWithParams1(tarch::la::equals((*data.at(0)->values)(1), 1.01000000000000000888), (*data.at(0)->values)(1));
  validateWithParams1(tarch::la::equals((*data.at(0)->values)(2), 1.02000000000000001776), (*data.at(0)->values)(2));
  validateWithParams1(tarch::la::equals((*data.at(0)->values)(3), 1.03000000000000002665), (*data.at(0)->values)(3));
  validateWithParams1(tarch::la::equals((*data.at(1)->values)(0), 0.199000000000000010214), (*data.at(1)->values)(0));
  validateWithParams1(tarch::la::equals((*data.at(1)->values)(1), 0.199000000000000010214), (*data.at(1)->values)(1));
  validateWithParams1(tarch::la::equals((*data.at(1)->values)(2), 0.199000000000000010214), (*data.at(1)->values)(2));
  
  
  utils::DynVector newdvalues;
  newdvalues.append(10.0);
  newdvalues.append(10.0);
  newdvalues.append(10.0);
  newdvalues.append(10.0);
  data.begin()->second->values = &newdvalues;
  
  pp.performPostProcessing(data);
  
  
  validateWithParams1(tarch::la::equals((*data.at(0)->values)(0), -5.63855295490201413600e-01), (*data.at(0)->values)(0));
  validateWithParams1(tarch::la::equals((*data.at(0)->values)(1), 6.09906404008707880848e-01), (*data.at(0)->values)(1));
  validateWithParams1(tarch::la::equals((*data.at(0)->values)(2), 1.78366810350762250437e+00), (*data.at(0)->values)(2));
  validateWithParams1(tarch::la::equals((*data.at(0)->values)(3), 2.95742980300653179881e+00), (*data.at(0)->values)(3));
  validateWithParams1(tarch::la::equals((*data.at(1)->values)(0), 8.27975917496077962188e-02), (*data.at(1)->values)(0));
  validateWithParams1(tarch::la::equals((*data.at(1)->values)(1), 8.27975917496077962188e-02), (*data.at(1)->values)(1));
  validateWithParams1(tarch::la::equals((*data.at(1)->values)(2), 8.27975917496077962188e-02), (*data.at(1)->values)(2));
  
}

#endif // not PRECICE_NO_MPI

}}}// namespace precice, cplscheme, tests

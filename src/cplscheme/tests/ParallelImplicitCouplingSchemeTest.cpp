#include "ParallelImplicitCouplingSchemeTest.hpp"
#include "cplscheme/ParallelCouplingScheme.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"
#include "cplscheme/config/PostProcessingConfiguration.hpp"
#include "cplscheme/impl/ConvergenceMeasure.hpp"
#include "cplscheme/impl/AbsoluteConvergenceMeasure.hpp"
#include "cplscheme/impl/MinIterationConvergenceMeasure.hpp"
#include "cplscheme/impl/IQNILSPostProcessing.hpp"
#include "cplscheme/impl/MVQNPostProcessing.hpp"
#include "cplscheme/impl/BaseQNPostProcessing.hpp"
#include "cplscheme/impl/ConstantPreconditioner.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "cplscheme/impl/SharedPointer.hpp"
#include "cplscheme/Constants.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "m2n/GatherScatterCommunication.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "m2n/M2N.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"
#include "xml/XMLTag.hpp"
#include <Eigen/Core>
#include "utils/EigenHelperFunctions.hpp"
#include "math/math.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::cplscheme::tests::ParallelImplicitCouplingSchemeTest)

namespace precice {
namespace cplscheme {
namespace tests {

logging::Logger ParallelImplicitCouplingSchemeTest::
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
  _pathToTests = utils::getPathToSources() + "/cplscheme/tests/";
}

void ParallelImplicitCouplingSchemeTest:: run ()
{
# ifndef PRECICE_NO_MPI
  PRECICE_MASTER_ONLY {
    testMethod(testParseConfigurationWithRelaxation);
    testMethod(testMVQNPP);
    testMethod(testVIQNPP);
  }
  typedef utils::Parallel Par;
  if (Par::getCommunicatorSize() > 1){
    // Do only use process 0 and 1 for the following tests
    MPI_Comm comm = Par::getRestrictedCommunicator({0, 1});
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
  TRACE();
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
  validate(cplSchemeConfig._postProcConfig->getPostProcessing().get() != nullptr);
  meshConfig->setMeshSubIDs();
}

void ParallelImplicitCouplingSchemeTest:: testInitializeData()
{
  TRACE();
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
  mesh->createVertex(Eigen::Vector3d::Zero());
  mesh->allocateDataValues();
  meshConfig.addMesh(mesh);

  // Create all parameters necessary to create a ParallelImplicitCouplingScheme object
  com::PtrCommunication communication(new com::MPIDirectCommunication);
  m2n::PtrM2N globalCom(new m2n::M2N(communication, m2n::DistributedComFactory::SharedPointer()));
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
    mesh->data()[1]->getID(), false, false, minIterationConvMeasure1 );
  cplScheme.addConvergenceMeasure (
    mesh->data()[0]->getID(), false, false, minIterationConvMeasure2 );
  connect(nameParticipant0, nameParticipant1, nameLocalParticipant, globalCom);

  std::string writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterationCheckpoint(constants::actionReadIterationCheckpoint());

  cplScheme.initialize(0.0, 0);

  if (nameLocalParticipant == nameParticipant0){
    validate(cplScheme.isActionRequired(constants::actionWriteInitialData()));
    mesh->data(0)->values() = Eigen::VectorXd::Constant(1, 4.0);
    cplScheme.performedAction(constants::actionWriteInitialData());
    cplScheme.initializeData();
    validate(cplScheme.hasDataBeenExchanged());
    auto& values = mesh->data(1)->values();
    validateWithParams1(math::equals(values, Eigen::Vector3d(1.0, 2.0, 3.0)), values);

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
    auto& values = mesh->data(0)->values();
    validate(cplScheme.isActionRequired(constants::actionWriteInitialData()));
    Eigen::VectorXd v(3); v << 1.0, 2.0, 3.0;
    mesh->data(1)->values() = v;
    cplScheme.performedAction(constants::actionWriteInitialData());
    validateWithParams1(math::equals(values(0), 0.0), values);
    cplScheme.initializeData();
    validate(cplScheme.hasDataBeenExchanged());
    validateWithParams1(math::equals(values(0), 4.0), values);

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
  utils::Parallel::clearGroups();
}

void ParallelImplicitCouplingSchemeTest:: connect
(
  const std::string&      participant0,
  const std::string&      participant1,
  const std::string&      localParticipant,
  m2n::PtrM2N& communication ) const
{
  assertion ( communication.use_count() > 0 );
  assertion ( not communication->isConnected() );
  utils::Parallel::splitCommunicator( localParticipant );
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
  TRACE();

  //use two vectors and see if underrelaxation works

  double initialRelaxation = 0.01;
  int    maxIterationsUsed = 50;
  int    timestepsReused = 6;
  int filter = impl::BaseQNPostProcessing::QR1FILTER;
  double singularityLimit = 1e-10;
  bool enforceInitialRelaxation = false;
  std::vector<int> dataIDs;
  dataIDs.push_back(0);
  dataIDs.push_back(1);
  std::vector<double> factors;
  factors.resize(2,1.0);
  impl::PtrPreconditioner prec(new impl::ConstantPreconditioner(factors));

  std::map<int, double> scalings;
  scalings.insert(std::make_pair(0,1.0));
  scalings.insert(std::make_pair(1,1.0));
  mesh::PtrMesh dummyMesh ( new mesh::Mesh("dummyMesh", 3, false) );

  cplscheme::impl::IQNILSPostProcessing pp(initialRelaxation, enforceInitialRelaxation, maxIterationsUsed,
                                           timestepsReused, filter, singularityLimit, dataIDs, prec);


  Eigen::VectorXd dvalues;
  Eigen::VectorXd dcol1;
  Eigen::VectorXd fvalues;
  Eigen::VectorXd fcol1;

  //init displacements
  utils::append(dvalues, 1.0);
  utils::append(dvalues, 2.0);
  utils::append(dvalues, 3.0);
  utils::append(dvalues, 4.0);

  utils::append(dcol1, 1.0);
  utils::append(dcol1, 1.0);
  utils::append(dcol1, 1.0);
  utils::append(dcol1, 1.0);

  PtrCouplingData dpcd(new CouplingData(&dvalues,dummyMesh,false,1));

  //init forces
  utils::append(fvalues, 0.1);
  utils::append(fvalues, 0.1);
  utils::append(fvalues, 0.1);
  utils::append(fvalues, 0.1);

  utils::append(fcol1, 0.2);
  utils::append(fcol1, 0.2);
  utils::append(fcol1, 0.2);
  utils::append(fcol1, 0.2);

  PtrCouplingData fpcd(new CouplingData(&fvalues,dummyMesh,false,1));

  DataMap data;
  data.insert(std::pair<int,PtrCouplingData>(0,dpcd));
  data.insert(std::pair<int,PtrCouplingData>(1,fpcd));

//  for (DataMap::value_type& pair : data){
//    std::cout << *pair.second->values << "\n";
//    std::cout << pair.second->oldValues << "\n";
//  }

  pp.initialize(data);

  dpcd->oldValues.col(0) = dcol1;
  fpcd->oldValues.col(0) = fcol1;

  pp.performPostProcessing(data);

  validateWithParams1(math::equals((*data.at(0)->values)(0), 1.00), (*data.at(0)->values)(0));
  validateWithParams1(math::equals((*data.at(0)->values)(1), 1.01), (*data.at(0)->values)(1));
  validateWithParams1(math::equals((*data.at(0)->values)(2), 1.02), (*data.at(0)->values)(2));
  validateWithParams1(math::equals((*data.at(0)->values)(3), 1.03), (*data.at(0)->values)(3));
  validateWithParams1(math::equals((*data.at(1)->values)(0), 0.199), (*data.at(1)->values)(0));
  validateWithParams1(math::equals((*data.at(1)->values)(1), 0.199), (*data.at(1)->values)(1));
  validateWithParams1(math::equals((*data.at(1)->values)(2), 0.199), (*data.at(1)->values)(2));
  validateWithParams1(math::equals((*data.at(1)->values)(3), 0.199), (*data.at(1)->values)(3));

  Eigen::VectorXd newdvalues;
  utils::append(newdvalues, 10.0);
  utils::append(newdvalues, 10.0);
  utils::append(newdvalues, 10.0);
  utils::append(newdvalues, 10.0);
  data.begin()->second->values = &newdvalues;

  pp.performPostProcessing(data);

  validateWithParams1(math::equals((*data.at(0)->values)(0), -5.63401340929692295845e-01), (*data.at(0)->values)(0));
  validateWithParams1(math::equals((*data.at(0)->values)(1), 6.10309919173607440257e-01), (*data.at(0)->values)(1));
  validateWithParams1(math::equals((*data.at(0)->values)(2), 1.78402117927690717636e+00), (*data.at(0)->values)(2));
  validateWithParams1(math::equals((*data.at(0)->values)(3), 2.95773243938020513610e+00), (*data.at(0)->values)(3));
  validateWithParams1(math::equals((*data.at(1)->values)(0), 8.28025852497733944046e-02), (*data.at(1)->values)(0));
  validateWithParams1(math::equals((*data.at(1)->values)(1), 8.28025852497733944046e-02), (*data.at(1)->values)(1));
  validateWithParams1(math::equals((*data.at(1)->values)(2), 8.28025852497733944046e-02), (*data.at(1)->values)(2));
  validateWithParams1(math::equals((*data.at(1)->values)(3), 8.28025852497733944046e-02), (*data.at(1)->values)(3));

}


void ParallelImplicitCouplingSchemeTest:: testMVQNPP()
{
  TRACE();
  
  //use two vectors and see if underrelaxation works
  
  double initialRelaxation = 0.01;
  int    maxIterationsUsed = 50;
  int    timestepsReused = 6;
  int    reusedTimestepsAtRestart = 0;
  int    chunkSize = 0;
  int filter = impl::BaseQNPostProcessing::QR1FILTER;
  int restartType = impl::MVQNPostProcessing::NO_RESTART;
  double singularityLimit = 1e-10;
  double svdTruncationEps = 0.0;
  bool enforceInitialRelaxation = false;
  bool alwaysBuildJacobian = false;
  std::vector<int> dataIDs;
  dataIDs.push_back(0);
  dataIDs.push_back(1);
  std::vector<double> factors;
  factors.resize(2,1.0);
  impl::PtrPreconditioner prec(new impl::ConstantPreconditioner(factors));
  mesh::PtrMesh dummyMesh ( new mesh::Mesh("dummyMesh", 3, false) );

  
  cplscheme::impl::MVQNPostProcessing pp(initialRelaxation, enforceInitialRelaxation, maxIterationsUsed,
                                         timestepsReused, filter, singularityLimit, dataIDs, prec, alwaysBuildJacobian,
                                         restartType, chunkSize, reusedTimestepsAtRestart, svdTruncationEps);
  
  Eigen::VectorXd dvalues;
  Eigen::VectorXd dcol1;
  Eigen::VectorXd fvalues;
  Eigen::VectorXd fcol1;

  //init displacements
  utils::append(dvalues, 1.0);
  utils::append(dvalues, 2.0);
  utils::append(dvalues, 3.0);
  utils::append(dvalues, 4.0);
  
  utils::append(dcol1, 1.0);
  utils::append(dcol1, 1.0);
  utils::append(dcol1, 1.0);
  utils::append(dcol1, 1.0);

  PtrCouplingData dpcd(new CouplingData(&dvalues,dummyMesh,false,1));

  //init forces
  utils::append(fvalues, 0.1);
  utils::append(fvalues, 0.1);
  utils::append(fvalues, 0.1);
  utils::append(fvalues, 0.1);

  utils::append(fcol1, 0.2);
  utils::append(fcol1, 0.2);
  utils::append(fcol1, 0.2);
  utils::append(fcol1, 0.2);

  PtrCouplingData fpcd(new CouplingData(&fvalues,dummyMesh,false,1));
  
  DataMap data;
  data.insert(std::pair<int,PtrCouplingData>(0,dpcd));
  data.insert(std::pair<int,PtrCouplingData>(1,fpcd));
  
//  for (DataMap::value_type& pair : data){
//    std::cout << *pair.second->values << "\n";
//    std::cout << pair.second->oldValues << "\n";
//  }
  
  pp.initialize(data);
  
  dpcd->oldValues.col(0) = dcol1;
  fpcd->oldValues.col(0) = fcol1;
  
  pp.performPostProcessing(data);
  
  validateWithParams1(math::equals((*data.at(0)->values)(0), 1.00000000000000000000), (*data.at(0)->values)(0));
  validateWithParams1(math::equals((*data.at(0)->values)(1), 1.01000000000000000888), (*data.at(0)->values)(1));
  validateWithParams1(math::equals((*data.at(0)->values)(2), 1.02000000000000001776), (*data.at(0)->values)(2));
  validateWithParams1(math::equals((*data.at(0)->values)(3), 1.03000000000000002665), (*data.at(0)->values)(3));
  validateWithParams1(math::equals((*data.at(1)->values)(0), 0.199000000000000010214), (*data.at(1)->values)(0));
  validateWithParams1(math::equals((*data.at(1)->values)(1), 0.199000000000000010214), (*data.at(1)->values)(1));
  validateWithParams1(math::equals((*data.at(1)->values)(2), 0.199000000000000010214), (*data.at(1)->values)(2));
  validateWithParams1(math::equals((*data.at(1)->values)(3), 0.199000000000000010214), (*data.at(1)->values)(3));
    
  Eigen::VectorXd newdvalues;
  utils::append(newdvalues, 10.0);
  utils::append(newdvalues, 10.0);
  utils::append(newdvalues, 10.0);
  utils::append(newdvalues, 10.0);

  data.begin()->second->values = &newdvalues;
  
  pp.performPostProcessing(data);
  
  validateWithParams1(math::equals((*data.at(0)->values)(0), -5.63401340929695848558e-01), (*data.at(0)->values)(0));
  validateWithParams1(math::equals((*data.at(0)->values)(1), 6.10309919173602111186e-01), (*data.at(0)->values)(1));
  validateWithParams1(math::equals((*data.at(0)->values)(2), 1.78402117927690184729e+00), (*data.at(0)->values)(2));
  validateWithParams1(math::equals((*data.at(0)->values)(3), 2.95773243938020247157e+00), (*data.at(0)->values)(3));
  validateWithParams1(math::equals((*data.at(1)->values)(0), 8.28025852497733250157e-02), (*data.at(1)->values)(0));
  validateWithParams1(math::equals((*data.at(1)->values)(1), 8.28025852497733250157e-02), (*data.at(1)->values)(1));
  validateWithParams1(math::equals((*data.at(1)->values)(2), 8.28025852497733250157e-02), (*data.at(1)->values)(2));
  validateWithParams1(math::equals((*data.at(1)->values)(3), 8.28025852497733250157e-02), (*data.at(1)->values)(3));

}

#endif // not PRECICE_NO_MPI

}}}// namespace precice, cplscheme, tests


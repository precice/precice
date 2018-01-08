#include "cplscheme/SerialCouplingScheme.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"
#include "cplscheme/Constants.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "com/Communication.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "m2n/M2N.hpp"
#include "m2n/GatherScatterCommunication.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"
#include "xml/XMLTag.hpp"
#include <Eigen/Core>

#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::cplscheme;

BOOST_AUTO_TEST_SUITE(CplSchemeTests)

struct ExplicitCouplingSchemeFixture  // TODO fixtures in cplscheme/boosttests are a candidate for refactoring, lots of copy paste code.
{
  std::string _pathToTests;

  ExplicitCouplingSchemeFixture(){
    _pathToTests = utils::getPathToSources() + "/cplscheme/boosttests/";
  }

  void runSimpleExplicitCoupling
  (
      CouplingScheme&                cplScheme,
      const std::string&             participantName,
      const mesh::MeshConfiguration& meshConfig )
  {
    BOOST_TEST(meshConfig.meshes().size() == 1);
    mesh::PtrMesh mesh = meshConfig.meshes()[0];
    BOOST_TEST(mesh->data().size() == 2);
    auto& dataValues0 = mesh->data()[0]->values();
    auto& dataValues1 = mesh->data()[1]->values();
    BOOST_TEST(mesh->vertices().size() > 0);
    mesh::Vertex& vertex = mesh->vertices()[0];
    double valueData0 = 1.0;
    Eigen::VectorXd valueData1 = Eigen::VectorXd::Constant(3, 1.0 );

    double computedTime = 0.0;
    int computedTimesteps = 0;

    if ( participantName == std::string("participant0") ) {
      cplScheme.initialize ( 0.0, 1 );
      BOOST_TEST(not cplScheme.hasDataBeenExchanged());
      BOOST_TEST(not cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));
      BOOST_TEST(not cplScheme.isActionRequired(constants::actionReadIterationCheckpoint()));
      BOOST_TEST(not cplScheme.isCouplingTimestepComplete());
      BOOST_TEST(cplScheme.isCouplingOngoing());
      while ( cplScheme.isCouplingOngoing() ) {
        dataValues0(vertex.getID()) = valueData0;
        computedTime += cplScheme.getNextTimestepMaxLength();
        computedTimesteps ++;
        cplScheme.addComputedTime ( cplScheme.getNextTimestepMaxLength() );
        cplScheme.advance();
        BOOST_TEST(cplScheme.isCouplingTimestepComplete());
        BOOST_TEST(testing::equals(computedTime, cplScheme.getTime()));
        BOOST_TEST(computedTimesteps == cplScheme.getTimesteps()-1);
        BOOST_TEST(not cplScheme.isActionRequired("WriteIterationCheckpoint"));
        BOOST_TEST(not cplScheme.isActionRequired("ReadIterationCheckpoint"));
        BOOST_TEST(cplScheme.isCouplingTimestepComplete());
        if ( cplScheme.isCouplingOngoing() ) {
          // No receive takes place for the participant that has started the
          // coupled simulation, in the last advance call
          Eigen::VectorXd value = dataValues1.segment(vertex.getID() * 3, 3);
          BOOST_TEST(testing::equals(value, valueData1));
        }
        BOOST_TEST(cplScheme.hasDataBeenExchanged());
        // Increment data values, to test if send/receive operations are also
        // correct in following timesteps.
        valueData0 += 1.0;
        valueData1 += Eigen::VectorXd::Constant(3, 1.0 );
      }
      cplScheme.finalize();
      // Validate results
      BOOST_TEST(testing::equals(computedTime, 1.0));
      BOOST_TEST(computedTimesteps == 10);
      BOOST_TEST(not cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"));
      BOOST_TEST(not cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"));
      BOOST_TEST(cplScheme.isCouplingTimestepComplete());
      BOOST_TEST(not cplScheme.isCouplingOngoing());
      BOOST_TEST(cplScheme.getNextTimestepMaxLength() > 0.0);
    }
    else if ( participantName == std::string("participant1") ) {
      cplScheme.initialize ( 0.0, 1 );
      BOOST_TEST(cplScheme.hasDataBeenExchanged());
      double value = dataValues0(vertex.getID());
      BOOST_TEST(testing::equals(value, valueData0));
      valueData0 += 1.0;
      BOOST_TEST(not cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"));
      BOOST_TEST(not cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"));
      BOOST_TEST(not cplScheme.isCouplingTimestepComplete());
      BOOST_TEST(cplScheme.isCouplingOngoing());
      while ( cplScheme.isCouplingOngoing() ) {
        dataValues1.segment(vertex.getID()*3, 3) = valueData1;
        computedTime += cplScheme.getNextTimestepMaxLength();
        computedTimesteps ++;
        cplScheme.addComputedTime ( cplScheme.getNextTimestepMaxLength() );
        cplScheme.advance();
        BOOST_TEST(testing::equals(computedTime, cplScheme.getTime()));
        BOOST_TEST(computedTimesteps == cplScheme.getTimesteps()-1);
        BOOST_TEST(not cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"));
        BOOST_TEST(not cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"));
        BOOST_TEST(cplScheme.isCouplingTimestepComplete());
        if ( cplScheme.isCouplingOngoing() ) {
          // The participant not starting the coupled simulation does neither
          // receive nor send data in the last call to advance
          BOOST_TEST(cplScheme.hasDataBeenExchanged());
          double value = dataValues0[vertex.getID()];
          BOOST_TEST(testing::equals(value, valueData0));
        }
        valueData0 += 1.0;
        valueData1 += Eigen::VectorXd::Constant(3, 1.0 );
      }
      cplScheme.finalize ();
      // Validate results
      BOOST_TEST(testing::equals(computedTime, 1.0));
      BOOST_TEST(computedTimesteps == 10);
      BOOST_TEST(not cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"));
      BOOST_TEST(not cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"));
      BOOST_TEST(cplScheme.isCouplingTimestepComplete());
      BOOST_TEST(not cplScheme.isCouplingOngoing());
      BOOST_TEST(cplScheme.getNextTimestepMaxLength() > 0.0);
    }
  }

  void runExplicitCouplingWithSubcycling
  (
      CouplingScheme&                cplScheme,
      const std::string&             participantName,
      const mesh::MeshConfiguration& meshConfig )
  {
    BOOST_TEST(meshConfig.meshes().size() == 1);
    mesh::PtrMesh mesh = meshConfig.meshes()[0];
    BOOST_TEST(mesh->data().size() == 2);
    BOOST_TEST(mesh->vertices().size() > 0);
    mesh::Vertex& vertex = mesh->vertices()[0];
    double valueData0 = 1.0;
    Eigen::VectorXd valueData1 = Eigen::VectorXd::Constant(3, 1.0);
    auto& dataValues0 = mesh->data()[0]->values();
    auto& dataValues1 = mesh->data()[1]->values();

    double computedTime = 0.0;
    int computedTimesteps = 0;
    std::string nameParticipant0 ( "participant0" );
    std::string nameParticipant1 ( "participant1" );
    assertion ( (participantName == nameParticipant0) ||
        (participantName == nameParticipant1) );
    if ( participantName == nameParticipant0 ) {
      cplScheme.initialize ( 0.0, 1 );
      double dtDesired = cplScheme.getNextTimestepMaxLength() / 2.0;
      double dtUsed = dtDesired;
      BOOST_TEST(not cplScheme.hasDataBeenExchanged());
      BOOST_TEST(not cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"));
      BOOST_TEST(not cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"));
      BOOST_TEST(not cplScheme.isCouplingTimestepComplete());
      BOOST_TEST(cplScheme.isCouplingOngoing());
      while ( cplScheme.isCouplingOngoing() ) {
        dataValues0(vertex.getID()) = valueData0;
        computedTime += dtUsed;
        computedTimesteps ++;
        cplScheme.addComputedTime(dtUsed);
        cplScheme.advance();
        // If the dt from preCICE is larger than the desired one, do subcycling,
        // else, use the dt from preCICE
        dtUsed = cplScheme.getNextTimestepMaxLength() > dtDesired
            ? dtDesired
                : cplScheme.getNextTimestepMaxLength();
        BOOST_TEST(testing::equals(computedTime, cplScheme.getTime()));
        BOOST_TEST(not cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"));
        BOOST_TEST(not cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"));
        if ( computedTimesteps % 2 == 0 ) {
          // Data exchange takes only place at every second local timestep,
          // since a subcycling of 2 is used.
          BOOST_TEST(cplScheme.isCouplingTimestepComplete());
          if ( cplScheme.isCouplingOngoing() ) {
            // No receive takes place for the participant that has started the
            // coupled simulation, in the last advance call.
            Eigen::VectorXd value = dataValues1.segment(vertex.getID()*3, 3);
            BOOST_TEST(testing::equals(value, valueData1));
          }
          BOOST_TEST(cplScheme.hasDataBeenExchanged());
          // Increment data values, to test if send/receive operations are also
          // correct in following timesteps.
          valueData0 += 1.0;
          valueData1 += Eigen::VectorXd::Constant(3, 1.0);
        }
        else {
          BOOST_TEST(not cplScheme.isCouplingTimestepComplete());
        }
      }
      cplScheme.finalize ();
      BOOST_TEST(testing::equals(computedTime, 1.0));
      BOOST_TEST(computedTimesteps == 20);
      BOOST_TEST(not cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"));
      BOOST_TEST(not cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"));
      BOOST_TEST(cplScheme.isCouplingTimestepComplete());
      BOOST_TEST(not cplScheme.isCouplingOngoing());
      BOOST_TEST(cplScheme.getNextTimestepMaxLength() > 0.0);
    }
    else if(participantName == nameParticipant1) {
      // Start coupling
      cplScheme.initialize ( 0.0, 1 );
      // Validate current coupling status
      BOOST_TEST(cplScheme.hasDataBeenExchanged());
      BOOST_TEST(testing::equals(dataValues0(vertex.getID()), valueData0));
      valueData0 += 1.0;
      BOOST_TEST(not cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"));
      BOOST_TEST(not cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"));
      BOOST_TEST(not cplScheme.isCouplingTimestepComplete());
      BOOST_TEST(cplScheme.isCouplingOngoing());
      while ( cplScheme.isCouplingOngoing() ) {
        dataValues1.segment(vertex.getID()*3, 3) = valueData1;
        computedTime += cplScheme.getNextTimestepMaxLength ();
        computedTimesteps ++;
        cplScheme.addComputedTime ( cplScheme.getNextTimestepMaxLength() );
        cplScheme.advance();
        BOOST_TEST(testing::equals(computedTime, cplScheme.getTime()));
        BOOST_TEST(computedTimesteps == cplScheme.getTimesteps()-1);
        BOOST_TEST(not cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"));
        BOOST_TEST(not cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"));
        BOOST_TEST(cplScheme.isCouplingTimestepComplete());
        if ( cplScheme.isCouplingOngoing() ) {
          // The participant not starting the coupled simulation does neither
          // receive nor send data in the last call to advance
          BOOST_TEST(cplScheme.hasDataBeenExchanged());
          BOOST_TEST(testing::equals(dataValues0(vertex.getID()), valueData0));
          BOOST_TEST(cplScheme.hasDataBeenExchanged());
        }
        valueData0 += 1.0;
        valueData1 += Eigen::VectorXd::Constant(3, 1.0);
      }
      cplScheme.finalize ();
      BOOST_TEST(testing::equals(computedTime, 1.0));
      BOOST_TEST(computedTimesteps == 10);
      BOOST_TEST(not cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));
      BOOST_TEST(not cplScheme.isActionRequired(constants::actionReadIterationCheckpoint()));
      BOOST_TEST(cplScheme.isCouplingTimestepComplete());
      BOOST_TEST(not cplScheme.isCouplingOngoing());
      BOOST_TEST(cplScheme.getNextTimestepMaxLength() > 0.0);
    }
  }

  void connect( // TODO this function occurs in multiple tests. Move this to a common fixture? see https://github.com/precice/precice/issues/90
      const std::string&      participant0,
      const std::string&      participant1,
      const std::string&      localParticipant,
      m2n::PtrM2N& communication )
  {
    assertion(communication.use_count() > 0);
    assertion(not communication->isConnected());
    utils::Parallel::splitCommunicator( localParticipant );
    if ( participant0 == localParticipant ) {
      communication->requestMasterConnection ( participant1, participant0);
    }
    else {
      assertion ( participant1 == localParticipant );
      communication->acceptMasterConnection ( participant1, participant0);
    }
  }
};

BOOST_FIXTURE_TEST_SUITE(ExplicitCouplingSchemeTests, ExplicitCouplingSchemeFixture)

#ifndef PRECICE_NO_MPI

/// Test that runs on 2 processors.
BOOST_AUTO_TEST_CASE(testSimpleExplicitCoupling, * testing::MinRanks(2) * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))
{
  if (utils::Parallel::getCommunicatorSize() != 2) // only run test on ranks {0,1}, for other ranks return
    return;

  mesh::PropertyContainer::resetPropertyIDCounter();
  xml::XMLTag root = xml::getRootTag();
  mesh::PtrDataConfiguration dataConfig ( new mesh::DataConfiguration(root) );
  dataConfig->addData ( "data0", 1 );
  dataConfig->addData ( "data1", 3 );
  mesh::MeshConfiguration meshConfig ( root, dataConfig );
  mesh::PtrMesh mesh ( new mesh::Mesh("mesh", 3, false) );
  mesh->createData ( "data0", 1 );
  mesh->createData ( "data1", 3 );
  mesh->createVertex ( Eigen::Vector3d::Zero() );
  mesh->allocateDataValues ();
  meshConfig.addMesh ( mesh );

  com::PtrCommunication communication ( new com::MPIDirectCommunication() );
  m2n::PtrM2N globalCom( new m2n::M2N(communication,m2n::DistributedComFactory::SharedPointer()) );
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
  cplscheme::SerialCouplingScheme cplScheme(
      maxTime, maxTimesteps, timestepLength, 12, nameParticipant0,
      nameParticipant1, localParticipant, globalCom, dtMethod, BaseCouplingScheme::Explicit );
  cplScheme.addDataToSend ( mesh->data()[sendDataIndex], mesh , false );
  cplScheme.addDataToReceive( mesh->data()[receiveDataIndex], mesh , false );
  connect(nameParticipant0, nameParticipant1, localParticipant, globalCom);
  runSimpleExplicitCoupling( cplScheme, localParticipant, meshConfig );
}

/// Test that runs on 2 processors.
BOOST_AUTO_TEST_CASE(testConfiguredSimpleExplicitCoupling, * testing::MinRanks(2) * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))
{
  if (utils::Parallel::getCommunicatorSize() != 2) // only run test on ranks {0,1}, for other ranks return
    return;

  using namespace mesh;
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
  xml::XMLTag root = xml::getRootTag();
  PtrDataConfiguration dataConfig(new DataConfiguration(root));
  dataConfig->setDimensions(3);
  PtrMeshConfiguration meshConfig(new MeshConfiguration(root, dataConfig));
  meshConfig->setDimensions(3);
  m2n::M2NConfiguration::SharedPointer m2nConfig(new m2n::M2NConfiguration(root));
  CouplingSchemeConfiguration cplSchemeConfig(root, meshConfig, m2nConfig);

  xml::configure(root, configurationPath);
  meshConfig->setMeshSubIDs();
  m2n::PtrM2N m2n = m2nConfig->getM2N("participant0", "participant1");

  // some dummy mesh
  meshConfig->meshes()[0]->createVertex(Eigen::Vector3d(1.0, 1.0, 1.0));
  meshConfig->meshes()[0]->createVertex(Eigen::Vector3d(2.0, 1.0, -1.0));
  meshConfig->meshes()[0]->createVertex(Eigen::Vector3d(3.0, 1.0, 1.0));
  meshConfig->meshes()[0]->createVertex(Eigen::Vector3d(4.0, 1.0, -1.0));
  meshConfig->meshes()[0]->allocateDataValues();

  connect ( "participant0", "participant1", localParticipant, m2n );
  runSimpleExplicitCoupling ( *cplSchemeConfig.getCouplingScheme(localParticipant),
      localParticipant, *meshConfig );
}

/// Test that runs on 2 processors.
BOOST_AUTO_TEST_CASE(testExplicitCouplingFirstParticipantSetsDt, * testing::MinRanks(2) * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))
{
  if (utils::Parallel::getCommunicatorSize() != 2) // only run test on ranks {0,1}, for other ranks return
    return;

  using namespace mesh;
  std::string configurationPath ( _pathToTests + "explicit-coupling-scheme-2.xml" );
  std::string localParticipant ( "" );
  if ( utils::Parallel::getProcessRank() == 0 ){
    localParticipant = "participant0";
  }
  else if ( utils::Parallel::getProcessRank() == 1 ){
    localParticipant = "participant1";
  }

  xml::XMLTag root = xml::getRootTag();
  PtrDataConfiguration dataConfig(new DataConfiguration(root));
  dataConfig->setDimensions(3);
  PtrMeshConfiguration meshConfig(new MeshConfiguration(root, dataConfig));
  meshConfig->setDimensions(3);
  m2n::M2NConfiguration::SharedPointer m2nConfig(new m2n::M2NConfiguration(root));
  CouplingSchemeConfiguration cplSchemeConfig(root, meshConfig, m2nConfig);

  xml::configure(root, configurationPath);
  meshConfig->setMeshSubIDs();
  m2n::PtrM2N m2n = m2nConfig->getM2N("participant0", "participant1");

  // some dummy mesh
  meshConfig->meshes()[0]->createVertex(Eigen::Vector3d(1.0, 1.0, 1.0));
  meshConfig->meshes()[0]->createVertex(Eigen::Vector3d(2.0, 1.0, -1.0));
  meshConfig->meshes()[0]->createVertex(Eigen::Vector3d(3.0, 1.0, 1.0));
  meshConfig->meshes()[0]->createVertex(Eigen::Vector3d(4.0, 1.0, -1.0));
  meshConfig->meshes()[0]->allocateDataValues();

  connect("participant0", "participant1", localParticipant, m2n);
  CouplingScheme& cplScheme = *cplSchemeConfig.getCouplingScheme(localParticipant);

  double computedTime = 0.0;
  int computedTimesteps = 0;
  if ( localParticipant == std::string("participant0") ){
    double dt = 0.3;
    cplScheme.initialize ( 0.0, 1 );
    BOOST_TEST(not cplScheme.isCouplingTimestepComplete());
    BOOST_TEST(cplScheme.isCouplingOngoing());
    while ( cplScheme.isCouplingOngoing() ){
      computedTime += dt;
      computedTimesteps++;
      cplScheme.addComputedTime(dt);
      cplScheme.advance();
      BOOST_TEST(cplScheme.isCouplingTimestepComplete());
      BOOST_TEST(testing::equals(computedTime, cplScheme.getTime()));
      BOOST_TEST(computedTimesteps == cplScheme.getTimesteps()-1);
      if ( cplScheme.isCouplingOngoing() ){
        BOOST_TEST(cplScheme.hasDataBeenExchanged());
      }
    }
    cplScheme.finalize();
    BOOST_TEST(testing::equals(computedTime, 1.2));
    BOOST_TEST(computedTimesteps == 4);
    BOOST_TEST(cplScheme.isCouplingTimestepComplete());
    BOOST_TEST(not cplScheme.isCouplingOngoing());
  }
  else {
    assertion ( localParticipant == std::string("participant1"), localParticipant );
    cplScheme.initialize ( 0.0, 1 );
    BOOST_TEST(not cplScheme.isCouplingTimestepComplete());
    BOOST_TEST(cplScheme.isCouplingOngoing());
    while ( cplScheme.isCouplingOngoing() ){
      computedTime += cplScheme.getTimestepLength();
      computedTimesteps++;
      cplScheme.addComputedTime(cplScheme.getTimestepLength());
      cplScheme.advance();
      BOOST_TEST(cplScheme.isCouplingTimestepComplete());
      BOOST_TEST(testing::equals(computedTime, cplScheme.getTime()));
      BOOST_TEST(computedTimesteps == cplScheme.getTimesteps()-1);
      if ( cplScheme.isCouplingOngoing() ){
        BOOST_TEST(cplScheme.hasDataBeenExchanged());
      }
    }
    cplScheme.finalize();
    BOOST_TEST(testing::equals(computedTime, 1.2));
    BOOST_TEST(computedTimesteps == 4);
    BOOST_TEST(cplScheme.isCouplingTimestepComplete());
    BOOST_TEST(not cplScheme.isCouplingOngoing());
  }
}

/// Test that runs on 2 processors.
BOOST_AUTO_TEST_CASE(testSerialDataInitialization, * testing::MinRanks(2) * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))
{
  if (utils::Parallel::getCommunicatorSize() != 2) // only run test on ranks {0,1}, for other ranks return
    return;

  using namespace mesh;
  mesh::PropertyContainer::resetPropertyIDCounter();

  std::string configurationPath(_pathToTests + "serial-explicit-coupling-datainit.xml");

  std::string localParticipant("");
  if (utils::Parallel::getProcessRank() == 0){
    localParticipant = "participant0";
  }
  else if (utils::Parallel::getProcessRank() == 1){
    localParticipant = "participant1";
  }
  xml::XMLTag root = xml::getRootTag();
  PtrDataConfiguration dataConfig(new DataConfiguration(root));
  dataConfig->setDimensions(2);
  PtrMeshConfiguration meshConfig(new MeshConfiguration(root, dataConfig));
  meshConfig->setDimensions(2);
  m2n::M2NConfiguration::SharedPointer m2nConfig(new m2n::M2NConfiguration(root));
  CouplingSchemeConfiguration cplSchemeConfig(root, meshConfig, m2nConfig);

  xml::configure(root, configurationPath);
  meshConfig->setMeshSubIDs();
  m2n::PtrM2N m2n = m2nConfig->getM2N("participant0", "participant1");

  // some dummy mesh
  meshConfig->meshes()[0]->createVertex(Eigen::Vector2d(1.0, 1.0));
  meshConfig->meshes()[0]->createVertex(Eigen::Vector2d(2.0,-1.0));
  meshConfig->meshes()[0]->createVertex(Eigen::Vector2d(3.0, 1.0));
  meshConfig->meshes()[0]->createVertex(Eigen::Vector2d(4.0,-1.0));
  meshConfig->meshes()[0]->allocateDataValues();

  connect("participant0", "participant1", localParticipant, m2n);
  CouplingScheme& cplScheme = *cplSchemeConfig.getCouplingScheme(localParticipant);

  BOOST_TEST(meshConfig->meshes().size() == 1);
  mesh::PtrMesh mesh = meshConfig->meshes()[0];
  BOOST_TEST(mesh->data().size() == 3);
  auto& dataValues0 = mesh->data()[0]->values();
  auto& dataValues1 = mesh->data()[1]->values();
  auto& dataValues2 = mesh->data()[2]->values();

  if (localParticipant == std::string("participant0")){
    cplScheme.initialize(0.0, 1);
    BOOST_TEST(not cplScheme.isActionRequired(constants::actionWriteInitialData()));
    cplScheme.initializeData();
    BOOST_TEST(cplScheme.hasDataBeenExchanged());
    BOOST_TEST(testing::equals(dataValues0(0), 0.0));
    BOOST_TEST(testing::equals(dataValues1(0), 1.0));
    dataValues2(0) = 2.0;
    cplScheme.addComputedTime(cplScheme.getNextTimestepMaxLength());
    cplScheme.advance();
    BOOST_TEST(not cplScheme.isCouplingOngoing());
    cplScheme.finalize();
  }
  else if (localParticipant == std::string("participant1")){
    cplScheme.initialize(0.0, 1);
    BOOST_TEST(not cplScheme.hasDataBeenExchanged());
    BOOST_TEST(cplScheme.isActionRequired(constants::actionWriteInitialData()));
    dataValues1(0) = 1.0;
    cplScheme.performedAction(constants::actionWriteInitialData());
    cplScheme.initializeData();
    BOOST_TEST(testing::equals(dataValues2(0), 2.0));
    cplScheme.addComputedTime(cplScheme.getNextTimestepMaxLength());
    cplScheme.advance();
    BOOST_TEST(not cplScheme.isCouplingOngoing());
    cplScheme.finalize();
  }
}

/// Test that runs on 2 processors.
BOOST_AUTO_TEST_CASE(testParallelDataInitialization, * testing::MinRanks(2) * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))
{
  if (utils::Parallel::getCommunicatorSize() != 2) // only run test on ranks {0,1}, for other ranks return
    return;

  using namespace mesh;
  mesh::PropertyContainer::resetPropertyIDCounter();

  std::string configurationPath(_pathToTests + "parallel-explicit-coupling-datainit.xml");

  std::string localParticipant("");
  if (utils::Parallel::getProcessRank() == 0){
    localParticipant = "participant0";
  }
  else if (utils::Parallel::getProcessRank() == 1){
    localParticipant = "participant1";
  }
  xml::XMLTag root = xml::getRootTag();
  PtrDataConfiguration dataConfig(new DataConfiguration(root));
  dataConfig->setDimensions(2);
  PtrMeshConfiguration meshConfig(new MeshConfiguration(root, dataConfig));
  meshConfig->setDimensions(2);
  m2n::M2NConfiguration::SharedPointer m2nConfig(new m2n::M2NConfiguration(root));
  CouplingSchemeConfiguration cplSchemeConfig(root, meshConfig, m2nConfig);

  xml::configure(root, configurationPath);
  meshConfig->setMeshSubIDs();
  m2n::PtrM2N m2n = m2nConfig->getM2N("participant0", "participant1");

  // some dummy mesh
  meshConfig->meshes()[0]->createVertex(Eigen::Vector2d(1.0, 1.0));
  meshConfig->meshes()[0]->createVertex(Eigen::Vector2d(2.0,-1.0));
  meshConfig->meshes()[0]->createVertex(Eigen::Vector2d(3.0, 1.0));
  meshConfig->meshes()[0]->createVertex(Eigen::Vector2d(4.0,-1.0));
  meshConfig->meshes()[0]->allocateDataValues();

  connect("participant0", "participant1", localParticipant, m2n);
  CouplingScheme& cplScheme = *cplSchemeConfig.getCouplingScheme(localParticipant);

  BOOST_TEST(meshConfig->meshes().size() == 1);
  mesh::PtrMesh mesh = meshConfig->meshes()[0];
  BOOST_TEST(mesh->data().size() == 3);
  auto& dataValues0 = mesh->data()[0]->values();
  auto& dataValues1 = mesh->data()[1]->values();
  auto& dataValues2 = mesh->data()[2]->values();

  if (localParticipant == std::string("participant0")){
    cplScheme.initialize(0.0, 1);
    BOOST_TEST(cplScheme.isActionRequired(constants::actionWriteInitialData()));
    dataValues2(0) = 3.0;
    cplScheme.performedAction(constants::actionWriteInitialData());
    cplScheme.initializeData();
    BOOST_TEST(cplScheme.hasDataBeenExchanged());
    BOOST_TEST(testing::equals(dataValues0(0), 0.0));
    BOOST_TEST(testing::equals(dataValues1(0), 1.0));
    dataValues2(0) = 2.0;
    cplScheme.addComputedTime(cplScheme.getNextTimestepMaxLength());
    cplScheme.advance();
    BOOST_TEST(testing::equals(dataValues0(0), 4.0));
    BOOST_TEST(not cplScheme.isCouplingOngoing());
    cplScheme.finalize();
  }
  else if (localParticipant == std::string("participant1")){
    cplScheme.initialize(0.0, 1);
    BOOST_TEST(not cplScheme.hasDataBeenExchanged());
    BOOST_TEST(cplScheme.isActionRequired(constants::actionWriteInitialData()));
    dataValues1(0) = 1.0;
    cplScheme.performedAction(constants::actionWriteInitialData());
    cplScheme.initializeData();
    BOOST_TEST(cplScheme.hasDataBeenExchanged());
    BOOST_TEST(testing::equals(dataValues2(0), 3.0));
    dataValues0(0) = 4.0;
    cplScheme.addComputedTime(cplScheme.getNextTimestepMaxLength());
    cplScheme.advance();
    BOOST_TEST(testing::equals(dataValues2(0), 2.0));
    BOOST_TEST(not cplScheme.isCouplingOngoing());
    cplScheme.finalize();
  }
}

/// Test that runs on 2 processors.
BOOST_AUTO_TEST_CASE(testExplicitCouplingWithSubcycling, * testing::MinRanks(2) * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))
{
  if (utils::Parallel::getCommunicatorSize() != 2) // only run test on ranks {0,1}, for other ranks return
    return;

  mesh::PropertyContainer::resetPropertyIDCounter ();
  xml::XMLTag root = xml::getRootTag();
  mesh::PtrDataConfiguration dataConfig ( new mesh::DataConfiguration(root) );
  dataConfig->setDimensions(3);
  dataConfig->addData ( "data0", 1 );
  dataConfig->addData ( "data1", 3 );
  mesh::MeshConfiguration meshConfig ( root, dataConfig );
  meshConfig.setDimensions(3);
  mesh::PtrMesh mesh ( new mesh::Mesh("mesh", 3, false) );
  mesh->createData ( "data0", 1 );
  mesh->createData ( "data1", 3 );
  mesh->createVertex ( Eigen::Vector3d::Zero() );
  mesh->allocateDataValues ();
  meshConfig.addMesh ( mesh );

  com::PtrCommunication communication ( new com::MPIDirectCommunication );
  m2n::PtrM2N globalCom (new m2n::M2N(communication,m2n::DistributedComFactory::SharedPointer()));
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
  cplscheme::SerialCouplingScheme cplScheme (
      maxTime, maxTimesteps, timestepLength, 12, nameParticipant0,
      nameParticipant1, localParticipant, globalCom, dtMethod, BaseCouplingScheme::Explicit );
  cplScheme.addDataToSend ( mesh->data()[sendDataIndex], mesh , false);
  cplScheme.addDataToReceive ( mesh->data()[receiveDataIndex], mesh , false);
  connect(nameParticipant0, nameParticipant1, localParticipant, globalCom);
  runExplicitCouplingWithSubcycling ( cplScheme, localParticipant, meshConfig );
}

/// Test that runs on 2 processors.
BOOST_AUTO_TEST_CASE(testConfiguredExplicitCouplingWithSubcycling, * testing::MinRanks(2) * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))
{
  if (utils::Parallel::getCommunicatorSize() != 2) // only run test on ranks {0,1}, for other ranks return
    return;

  using namespace mesh;

  std::string configurationPath ( _pathToTests + "explicit-coupling-scheme-1.xml" );

  std::string localParticipant ( "" );
  if ( utils::Parallel::getProcessRank() == 0 ) {
    localParticipant = "participant0";
  }
  else if ( utils::Parallel::getProcessRank() == 1 ) {
    localParticipant = "participant1";
  }

  xml::XMLTag root = xml::getRootTag();
  PtrDataConfiguration dataConfig(new DataConfiguration(root));
  dataConfig->setDimensions(3);
  PtrMeshConfiguration meshConfig(new MeshConfiguration(root, dataConfig));
  meshConfig->setDimensions(3);
  m2n::M2NConfiguration::SharedPointer m2nConfig(new m2n::M2NConfiguration(root));
  CouplingSchemeConfiguration cplSchemeConfig(root, meshConfig, m2nConfig);

  xml::configure(root, configurationPath);
  meshConfig->setMeshSubIDs();
  m2n::PtrM2N m2n = m2nConfig->getM2N("participant0", "participant1");

  // some dummy mesh
  meshConfig->meshes()[0]->createVertex(Eigen::Vector3d(1.0, 1.0, 1.0));
  meshConfig->meshes()[0]->createVertex(Eigen::Vector3d(2.0,-1.0, 1.0));
  meshConfig->meshes()[0]->createVertex(Eigen::Vector3d(3.0, 1.0, 1.0));
  meshConfig->meshes()[0]->createVertex(Eigen::Vector3d(4.0,-1.0, 1.0));
  meshConfig->meshes()[0]->allocateDataValues();

  connect ( "participant0", "participant1", localParticipant, m2n );
  runExplicitCouplingWithSubcycling (
      *cplSchemeConfig.getCouplingScheme(localParticipant), localParticipant,
      *meshConfig );
}

# endif // not PRECICE_NO_MPI

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

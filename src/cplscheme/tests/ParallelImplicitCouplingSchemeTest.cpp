#include <Eigen/Core>
#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include "acceleration/ConstantRelaxationAcceleration.hpp"
#include "acceleration/config/AccelerationConfiguration.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "cplscheme/BaseCouplingScheme.hpp"
#include "cplscheme/Constants.hpp"
#include "cplscheme/CouplingData.hpp"
#include "cplscheme/ParallelCouplingScheme.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"
#include "cplscheme/impl/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "math/differences.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "precice/config/ParticipantConfiguration.hpp"
#include "testing/ParallelCouplingSchemeFixture.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "xml/XMLTag.hpp"

using namespace precice;
using namespace precice::cplscheme;

BOOST_AUTO_TEST_SUITE(CplSchemeTests)

struct ParallelImplicitCouplingSchemeFixture {
  std::string _pathToTests;

  ParallelImplicitCouplingSchemeFixture()
  {
    _pathToTests = testing::getPathToSources() + "/cplscheme/tests/";
  }
};

BOOST_FIXTURE_TEST_SUITE(ParallelImplicitCouplingSchemeTests, ParallelImplicitCouplingSchemeFixture)

#ifndef PRECICE_NO_MPI

BOOST_AUTO_TEST_CASE(testParseConfigurationWithRelaxation)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;

  std::string path(_pathToTests + "parallel-implicit-cplscheme-relax-const-config.xml");

  xml::XMLTag                          root = xml::getRootTag();
  PtrDataConfiguration                 dataConfig(new DataConfiguration(root));
  PtrMeshConfiguration                 meshConfig(new MeshConfiguration(root, dataConfig));
  m2n::M2NConfiguration::SharedPointer m2nConfig(
      new m2n::M2NConfiguration(root));
  precice::config::PtrParticipantConfiguration participantConfig(new precice::config::ParticipantConfiguration(root, meshConfig));
  CouplingSchemeConfiguration                  cplSchemeConfig(root, meshConfig, m2nConfig, participantConfig);

  xml::configure(root, xml::ConfigurationContext{}, path);
  BOOST_CHECK(cplSchemeConfig._accelerationConfig->getAcceleration().get());
}

BOOST_AUTO_TEST_CASE(testInitializeData)
{
  PRECICE_TEST("Participant0"_on(1_rank), "Participant1"_on(1_rank), Require::Events);
  testing::ConnectionOptions options;
  options.useOnlyPrimaryCom = true;
  auto m2n                  = context.connectPrimaryRanks("Participant0", "Participant1", options);

  xml::XMLTag root = xml::getRootTag();

  int dimensions = 3;

  // Create a data configuration, to simplify configuration of data
  mesh::PtrDataConfiguration dataConfig(new mesh::DataConfiguration(root));
  dataConfig->addData("Data0", mesh::Data::typeName::SCALAR);
  dataConfig->addData("Data1", mesh::Data::typeName::VECTOR);

  mesh::MeshConfiguration meshConfig(root, dataConfig);
  mesh::PtrMesh           mesh(new mesh::Mesh("Mesh", 3, testing::nextMeshID()));
  const auto              dataID0 = mesh->createData("Data0", 1, 0_dataID)->getID();
  const auto              dataID1 = mesh->createData("Data1", 3, 1_dataID)->getID();
  mesh->createVertex(Eigen::Vector3d::Zero());
  mesh->allocateDataValues();
  meshConfig.insertMeshToMeshDimensionsMap(mesh->getName(), mesh->getDimensions());
  meshConfig.addMesh(mesh);

  // Create all parameters necessary to create a ParallelImplicitCouplingScheme object
  double       maxTime        = 1.0;
  int          maxTimeWindows = 3;
  const double timeWindowSize = 0.1;
  const double timeStepSize   = timeWindowSize; // solver is not subcycling
  std::string  nameParticipant0("Participant0");
  std::string  nameParticipant1("Participant1");
  int          sendDataIndex              = -1;
  int          receiveDataIndex           = -1;
  bool         dataRequiresInitialization = false;
  if (context.isNamed(nameParticipant0)) {
    sendDataIndex              = dataID0;
    receiveDataIndex           = dataID1;
    dataRequiresInitialization = true;
  } else {
    sendDataIndex              = dataID1;
    receiveDataIndex           = dataID0;
    dataRequiresInitialization = true;
  }

  // Create the coupling scheme object
  const int              minIterations = 1;
  const int              maxIterations = 3;
  ParallelCouplingScheme cplScheme(maxTime, maxTimeWindows, timeWindowSize, nameParticipant0, nameParticipant1, context.name, m2n, constants::FIXED_TIME_WINDOW_SIZE, BaseCouplingScheme::Implicit, minIterations, maxIterations);

  using Fixture = testing::ParallelCouplingSchemeFixture;
  cplScheme.addDataToSend(mesh->data(sendDataIndex), mesh, dataRequiresInitialization, true);
  CouplingData *sendCouplingData = Fixture::getSendData(cplScheme, sendDataIndex);
  cplScheme.addDataToReceive(mesh->data(receiveDataIndex), mesh, dataRequiresInitialization, true);
  CouplingData *receiveCouplingData = Fixture::getReceiveData(cplScheme, receiveDataIndex);
  cplScheme.determineInitialDataExchange();

  if (context.isNamed(nameParticipant0)) {
    BOOST_TEST(testing::equals(receiveCouplingData->values(), Eigen::Vector3d(0.0, 0.0, 0.0)));
    BOOST_TEST(receiveCouplingData->values().size() == 3);
    BOOST_TEST(testing::equals(sendCouplingData->values()(0), 0.0));
    BOOST_TEST(sendCouplingData->values().size() == 1);
    BOOST_TEST(Fixture::isImplicitCouplingScheme(cplScheme));
    BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::InitializeData));
    sendCouplingData->setSampleAtTime(0, time::Sample{1, Eigen::VectorXd::Constant(1, 4.0)});
    cplScheme.markActionFulfilled(CouplingScheme::Action::InitializeData);
    cplScheme.initialize();
    BOOST_TEST(cplScheme.hasDataBeenReceived());
    BOOST_TEST(testing::equals(receiveCouplingData->values(), Eigen::Vector3d(1.0, 2.0, 3.0)));
    BOOST_TEST(receiveCouplingData->getPreviousIterationSize() == 3);
    BOOST_TEST(testing::equals(receiveCouplingData->previousIteration(), Eigen::Vector3d(0.0, 0.0, 0.0)));
    BOOST_TEST(sendCouplingData->getPreviousIterationSize() == 1);
    BOOST_TEST(testing::equals(sendCouplingData->previousIteration()(0), 4.0));
    while (cplScheme.isCouplingOngoing()) {
      if (cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint)) {
        cplScheme.markActionFulfilled(CouplingScheme::Action::WriteCheckpoint);
      }
      BOOST_TEST(cplScheme.getNextTimeStepMaxSize() == timeStepSize);
      sendCouplingData->setSampleAtTime(cplScheme.getTime() + timeStepSize, time::Sample{1, Eigen::VectorXd::Constant(1, 4.0)});
      cplScheme.addComputedTime(timeStepSize);
      cplScheme.firstSynchronization({});
      cplScheme.firstExchange();
      cplScheme.secondSynchronization();
      cplScheme.secondExchange();
      BOOST_TEST(cplScheme.hasDataBeenReceived());
      if (cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint)) {
        cplScheme.markActionFulfilled(CouplingScheme::Action::ReadCheckpoint);
      }
    }
  } else {
    BOOST_TEST(context.isNamed(nameParticipant1));
    BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::InitializeData));
    Eigen::VectorXd v(3);
    v << 1.0, 2.0, 3.0;
    sendCouplingData->setSampleAtTime(0, time::Sample{3, v});
    cplScheme.markActionFulfilled(CouplingScheme::Action::InitializeData);
    BOOST_TEST(testing::equals(receiveCouplingData->values()(0), 0.0));
    BOOST_TEST(receiveCouplingData->values().size() == 1);
    BOOST_TEST(testing::equals(sendCouplingData->values(), Eigen::Vector3d(1.0, 2.0, 3.0)));
    BOOST_TEST(sendCouplingData->values().size() == 3);
    cplScheme.initialize();
    BOOST_TEST(cplScheme.hasDataBeenReceived());
    BOOST_TEST(testing::equals(receiveCouplingData->values()(0), 4.0));
    BOOST_TEST(receiveCouplingData->getPreviousIterationSize() == 1);
    BOOST_TEST(testing::equals(receiveCouplingData->previousIteration()(0), 0.0));
    BOOST_TEST(testing::equals(sendCouplingData->values(), Eigen::Vector3d(1.0, 2.0, 3.0)));
    BOOST_TEST(sendCouplingData->getPreviousIterationSize() == 3);
    BOOST_TEST(testing::equals(sendCouplingData->previousIteration(), Eigen::Vector3d(1.0, 2.0, 3.0)));
    while (cplScheme.isCouplingOngoing()) {
      if (cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint)) {
        cplScheme.markActionFulfilled(CouplingScheme::Action::WriteCheckpoint);
      }
      BOOST_TEST(cplScheme.getNextTimeStepMaxSize() == timeStepSize);
      sendCouplingData->setSampleAtTime(cplScheme.getTime() + timeStepSize, time::Sample{3, v});
      cplScheme.addComputedTime(timeStepSize);
      cplScheme.firstSynchronization({});
      cplScheme.firstExchange();
      cplScheme.secondSynchronization();
      cplScheme.secondExchange();
      BOOST_TEST(cplScheme.hasDataBeenReceived());
      if (cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint)) {
        cplScheme.markActionFulfilled(CouplingScheme::Action::ReadCheckpoint);
      }
    }
  }
  cplScheme.finalize();
}

#endif // not PRECICE_NO_MPI

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

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
#include "cplscheme/impl/MinIterationConvergenceMeasure.hpp"
#include "cplscheme/impl/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/config/M2NConfiguration.hpp"
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

  int dimensions = 3;

  std::string path(_pathToTests + "parallel-implicit-cplscheme-relax-const-config.xml");

  xml::XMLTag          root = xml::getRootTag();
  PtrDataConfiguration dataConfig(new DataConfiguration(root));
  dataConfig->setDimensions(dimensions);
  PtrMeshConfiguration meshConfig(new MeshConfiguration(root, dataConfig));
  meshConfig->setDimensions(dimensions);
  m2n::M2NConfiguration::SharedPointer m2nConfig(
      new m2n::M2NConfiguration(root));
  precice::config::PtrParticipantConfiguration participantConfig(new precice::config::ParticipantConfiguration(root, meshConfig));
  participantConfig->setDimensions(dimensions);
  CouplingSchemeConfiguration cplSchemeConfig(root, meshConfig, m2nConfig, participantConfig);

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
  dataConfig->setDimensions(dimensions);
  dataConfig->addData("Data0", 1);
  dataConfig->addData("Data1", 3);

  mesh::MeshConfiguration meshConfig(root, dataConfig);
  meshConfig.setDimensions(dimensions);
  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 3, testing::nextMeshID()));
  const auto    dataID0 = mesh->createData("Data0", 1, 0_dataID)->getID();
  const auto    dataID1 = mesh->createData("Data1", 3, 1_dataID)->getID();
  mesh->createVertex(Eigen::Vector3d::Zero());
  mesh->allocateDataValues();
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
  int          extrapolationOrder         = 0;
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
  ParallelCouplingScheme cplScheme(
      maxTime, maxTimeWindows, timeWindowSize, 16, nameParticipant0, nameParticipant1,
      context.name, m2n, constants::FIXED_TIME_WINDOW_SIZE, BaseCouplingScheme::Implicit, 100, extrapolationOrder);

  using Fixture = testing::ParallelCouplingSchemeFixture;
  cplScheme.addDataToSend(mesh->data(sendDataIndex), mesh, dataRequiresInitialization);
  CouplingData *sendCouplingData = Fixture::getSendData(cplScheme, sendDataIndex);
  cplScheme.addDataToReceive(mesh->data(receiveDataIndex), mesh, dataRequiresInitialization);
  CouplingData *receiveCouplingData = Fixture::getReceiveData(cplScheme, receiveDataIndex);
  cplScheme.determineInitialDataExchange();

  // Add convergence measures
  int                                    minIterations = 3;
  cplscheme::impl::PtrConvergenceMeasure minIterationConvMeasure1(
      new cplscheme::impl::MinIterationConvergenceMeasure(minIterations));
  cplscheme::impl::PtrConvergenceMeasure minIterationConvMeasure2(
      new cplscheme::impl::MinIterationConvergenceMeasure(minIterations));
  cplScheme.addConvergenceMeasure(dataID1, false, false, minIterationConvMeasure1, true);
  cplScheme.addConvergenceMeasure(dataID0, false, false, minIterationConvMeasure2, true);

  if (context.isNamed(nameParticipant0)) {
    BOOST_TEST(testing::equals(receiveCouplingData->values(), Eigen::Vector3d(0.0, 0.0, 0.0)));
    BOOST_TEST(receiveCouplingData->values().size() == 3);
    BOOST_TEST(receiveCouplingData->getPreviousIterationSize() == 3);
    BOOST_TEST(testing::equals(sendCouplingData->values()(0), 0.0));
    BOOST_TEST(sendCouplingData->values().size() == 1);
    BOOST_TEST(sendCouplingData->getPreviousIterationSize() == 1);
    BOOST_TEST(Fixture::isImplicitCouplingScheme(cplScheme));
    BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::InitializeData));
    sendCouplingData->values() = Eigen::VectorXd::Constant(1, 4.0);
    cplScheme.markActionFulfilled(CouplingScheme::Action::InitializeData);
    cplScheme.initialize(0.0, 0);
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
      if (cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint)) {
        cplScheme.markActionFulfilled(CouplingScheme::Action::ReadCheckpoint);
      }
      cplScheme.addComputedTime(timeStepSize);
      cplScheme.firstSynchronization({});
      cplScheme.firstExchange();
      cplScheme.secondSynchronization();
      cplScheme.secondExchange();
      BOOST_TEST(cplScheme.hasDataBeenReceived());
    }
  } else {
    BOOST_TEST(context.isNamed(nameParticipant1));
    BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::InitializeData));
    Eigen::VectorXd v(3);
    v << 1.0, 2.0, 3.0;
    sendCouplingData->values() = v;
    cplScheme.markActionFulfilled(CouplingScheme::Action::InitializeData);
    BOOST_TEST(testing::equals(receiveCouplingData->values()(0), 0.0));
    BOOST_TEST(receiveCouplingData->values().size() == 1);
    BOOST_TEST(receiveCouplingData->getPreviousIterationSize() == 1);
    BOOST_TEST(testing::equals(sendCouplingData->values(), Eigen::Vector3d(1.0, 2.0, 3.0)));
    BOOST_TEST(sendCouplingData->values().size() == 3);
    BOOST_TEST(sendCouplingData->getPreviousIterationSize() == 3);
    cplScheme.initialize(0.0, 0);
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

BOOST_AUTO_TEST_SUITE(Extrapolation)
BOOST_AUTO_TEST_CASE(FirstOrder)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;

  PtrMesh mesh(new Mesh("MyMesh", 3, testing::nextMeshID()));
  PtrData data   = mesh->createData("MyData", 1, 0_dataID);
  int     dataID = data->getID();
  mesh->createVertex(Eigen::Vector3d::Zero());
  mesh->allocateDataValues();
  BOOST_TEST(data->values().size() == 1);

  const double          maxTime        = CouplingScheme::UNDEFINED_TIME;
  const int             maxTimeWindows = 1;
  const double          timeWindowSize = 1.0;
  std::string           first          = "First";
  std::string           second         = "Second";
  std::string           accessor       = second;
  com::PtrCommunication com(new com::MPIDirectCommunication());
  m2n::PtrM2N           globalCom(new m2n::M2N(com, m2n::DistributedComFactory::SharedPointer()));
  const int             maxIterations      = 1;
  const int             extrapolationOrder = 1;

  // Test first order extrapolation
  ParallelCouplingScheme scheme(maxTime, maxTimeWindows, timeWindowSize, 16, first, second,
                                accessor, globalCom, constants::FIXED_TIME_WINDOW_SIZE,
                                BaseCouplingScheme::Implicit, maxIterations, extrapolationOrder);

  using Fixture = testing::ParallelCouplingSchemeFixture;

  scheme.addDataToSend(data, mesh, true);
  Fixture::initializeStorages(scheme);
  CouplingData *cplData = Fixture::getSendData(scheme, dataID);
  BOOST_CHECK(cplData); // no nullptr
  BOOST_TEST(cplData->getSize() == 1);
  BOOST_TEST(cplData->getPreviousIterationSize() == 1);

  Fixture::moveToNextWindow(scheme);

  // data is uninitialized
  BOOST_TEST(testing::equals(cplData->values()(0), 0.0));
  BOOST_TEST(testing::equals(cplData->previousIteration()(0), 0.0));

  // start first window
  cplData->values()(0) = 1.0; // data provided at end of first window
  Fixture::setTimeWindows(scheme, scheme.getTimeWindows() + 1);
  Fixture::storeExtrapolationData(scheme);
  BOOST_TEST(testing::equals(cplData->values()(0), 1.0));

  // go to second window
  Fixture::moveToNextWindow(scheme); // uses first order extrapolation at end of first window
  BOOST_TEST(testing::equals(cplData->previousIteration()(0), 0.0));
  Fixture::storeIteration(scheme);
  BOOST_TEST(testing::equals(cplData->values()(0), 2.0)); // = 2*1 - 0
  BOOST_TEST(testing::equals(cplData->previousIteration()(0), 2.0));
  cplData->values()(0) = 4.0; // data provided at end of second window
  Fixture::setTimeWindows(scheme, scheme.getTimeWindows() + 1);
  Fixture::storeExtrapolationData(scheme);

  // go to third window
  Fixture::moveToNextWindow(scheme); // uses first order extrapolation (maximum allowed) at end of second window
  BOOST_TEST(testing::equals(cplData->previousIteration()(0), 2.0));
  Fixture::storeIteration(scheme);
  BOOST_TEST(testing::equals(cplData->values()(0), 7.0)); // = 2*4 - 1
  BOOST_TEST(testing::equals(cplData->previousIteration()(0), 7.0));
  cplData->values()(0) = 10.0; // data provided at end of third window
  Fixture::setTimeWindows(scheme, scheme.getTimeWindows() + 1);
  Fixture::storeExtrapolationData(scheme);

  // go to fourth window
  Fixture::moveToNextWindow(scheme); // uses first order extrapolation (maximum allowed) at end of third window
  BOOST_TEST(testing::equals(cplData->previousIteration()(0), 7.0));
  Fixture::storeIteration(scheme);
  BOOST_TEST(testing::equals(cplData->values()(0), 16.0)); // = 2*10 - 4
  BOOST_TEST(testing::equals(cplData->previousIteration()(0), 16.0));
}

/// Test that cplScheme gives correct results when applying extrapolation.
BOOST_AUTO_TEST_CASE(FirstOrderWithAcceleration)
{
  /**
   * Perform first order and constant relaxation acceleration
   *
   * Do two time windows with three iterations each.
   *
   * Each participant writes dummy data to other participant, received data is checked.
   *
   * Make sure that the following happens, if NOT converged (first two iterations):
   * 1. acceleration is performed
   * 2. participants receive correct (accelerated) data
   *
   * Make sure that the following happens, if converged (end of third iteration):
   * 1. old data is stored (we cannot access this from the coupling scheme, but we can deduct this from the extrapolated value)
   * 2. we move to the next window
   * 3. initial guess for first participant is computed via extrapolation from old data
   **/

  PRECICE_TEST("Participant0"_on(1_rank), "Participant1"_on(1_rank), Require::Events);
  testing::ConnectionOptions options;
  options.useOnlyPrimaryCom = true;
  auto m2n                  = context.connectPrimaryRanks("Participant0", "Participant1", options);

  xml::XMLTag root = xml::getRootTag();

  // Create a data configuration, to simplify configuration of data

  mesh::PtrDataConfiguration dataConfig(new mesh::DataConfiguration(root));
  const int                  geometrical_dimensions = 3; // 3d problem
  const int                  data_dimensions        = 1; // only one sample in data
  dataConfig->setDimensions(geometrical_dimensions);
  dataConfig->addData("Data0", data_dimensions);
  dataConfig->addData("Data1", data_dimensions);

  mesh::MeshConfiguration meshConfig(root, dataConfig);
  meshConfig.setDimensions(3);
  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", geometrical_dimensions, testing::nextMeshID()));
  const auto    dataID0 = mesh->createData("Data0", data_dimensions, 0_dataID)->getID();
  const auto    dataID1 = mesh->createData("Data1", data_dimensions, 1_dataID)->getID();
  mesh->createVertex(Eigen::Vector3d::Zero());
  mesh->allocateDataValues();
  meshConfig.addMesh(mesh);

  // Create all parameters necessary to create an ImplicitCouplingScheme object
  const double maxTime            = CouplingScheme::UNDEFINED_TIME;
  const int    maxTimeWindows     = 2;
  const double timeWindowSize     = 0.1;
  const int    maxIterations      = 3;
  const int    extrapolationOrder = 1;
  const double timeStepSize       = timeWindowSize; // solver is not subcycling
  std::string  first("Participant0");
  std::string  second("Participant1");
  int          sendDataIndex        = -1;
  int          receiveDataIndex     = -1;
  int          convergenceDataIndex = -1;

  BOOST_TEST(dataID0 == 0);
  BOOST_TEST(dataID1 == 1);

  if (context.isNamed(first)) {
    sendDataIndex        = dataID0;
    receiveDataIndex     = dataID1;
    convergenceDataIndex = receiveDataIndex;
  } else {
    sendDataIndex        = dataID1;
    receiveDataIndex     = dataID0;
    convergenceDataIndex = sendDataIndex;
  }

  // Create the coupling scheme object
  cplscheme::ParallelCouplingScheme cplScheme(
      maxTime, maxTimeWindows, timeWindowSize, 16, first, second,
      context.name, m2n, constants::FIXED_TIME_WINDOW_SIZE,
      BaseCouplingScheme::Implicit, maxIterations, extrapolationOrder);
  cplScheme.addDataToSend(mesh->data(sendDataIndex), mesh, false);
  cplScheme.addDataToReceive(mesh->data(receiveDataIndex), mesh, false);
  cplScheme.determineInitialDataExchange();

  // Add acceleration
  acceleration::PtrAcceleration ptrAcceleration(new acceleration::ConstantRelaxationAcceleration(0.5, std::vector<int>({sendDataIndex})));
  cplScheme.setAcceleration(ptrAcceleration);

  // Add convergence measures
  const int                              minIterations = maxIterations;
  cplscheme::impl::PtrConvergenceMeasure minIterationConvMeasure1(
      new cplscheme::impl::MinIterationConvergenceMeasure(minIterations));
  cplScheme.addConvergenceMeasure(convergenceDataIndex, false, false, minIterationConvMeasure1, true);

  cplScheme.initialize(0.0, 1);

  Eigen::VectorXd v(1); // buffer for data

  // write data is uninitialized
  BOOST_TEST(mesh->data(sendDataIndex)->values()(0) == 0);

  // first window
  for (int i = 0; i < maxIterations; i++) {
    // first, second and third iteration
    BOOST_TEST(cplScheme.isCouplingOngoing());

    if (context.isNamed(first)) {
      if (i == 0) {
        // data is uninitialized
        BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 0);
      } else if (i == 1) {
        // accelerated data from second participant: 0.5 * 0 + 0.5 * 2 = 1
        BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 1);
      } else if (i == 2) {
        // accelerated data from second participant: 0.5 * 1 + 0.5 * 2 = 1
        BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 1.5);
      }
    } else if (context.isNamed(second)) {
      if (i == 0) {
        // data is uninitialized
        BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 0);
      } else if (i == 1) {
        // accelerated data from first participant: 0.5 * 0 + 0.5 * 1 = 0.5
        BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 0.5);
      } else if (i == 2) {
        // accelerated data from first participant: 0.5 * 0.5 + 0.5 * 1 = 0.75
        BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 0.75);
      }
    }

    if (i == 0) {
      BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint));
      cplScheme.markActionFulfilled(CouplingScheme::Action::WriteCheckpoint);
      BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint));
    } else {
      BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint));
      BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint));
      cplScheme.markActionFulfilled(CouplingScheme::Action::ReadCheckpoint);
    }

    // write data to mesh
    if (context.isNamed(first)) {
      v << 1.0;
    } else if (context.isNamed(second)) {
      v << 2.0;
    }
    mesh->data(sendDataIndex)->values() = v;
    cplScheme.addComputedTime(timeStepSize);

    cplScheme.firstSynchronization({});
    cplScheme.firstExchange();
    cplScheme.secondSynchronization();
    cplScheme.secondExchange();

    if (i < maxIterations - 1) {
      BOOST_TEST(not cplScheme.isTimeWindowComplete());
    } else {
      // window complete since max iterations reached
      BOOST_TEST(cplScheme.isTimeWindowComplete());
    }
  }

  // second window
  for (int i = 0; i < maxIterations; i++) {
    // first, second and third iteration
    BOOST_TEST(cplScheme.isCouplingOngoing());
    if (context.isNamed(first)) {
      if (i == 0) {
        // first order extrapolation
        BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 4); // = 2*2 - 0
      } else if (i == 1) {
        // accelerated data from second participant
        BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 3.5); // = 0.5 * 4 + 0.5 * 3
      } else if (i == 2) {
        // accelerated data from second participant
        BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 3.25); // = 0.5 * 3.5 + 0.5 * 3
      }
    } else if (context.isNamed(second)) {
      if (i == 0) {
        // first order extrapolation
        BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 2); // = 2*1 - 0
      } else if (i == 1) {
        // accelerated data from first participant
        BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 2.5); // = 0.5 * 2 + 0.5 * 3
      } else if (i == 2) {
        // accelerated data from first participant
        BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 2.75); // = 0.5 * 2.5 + 0.5 * 3
      }
    }

    if (i == 0) {
      BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint));
      cplScheme.markActionFulfilled(CouplingScheme::Action::WriteCheckpoint);
      BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint));
    } else {
      BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint));
      BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint));
      cplScheme.markActionFulfilled(CouplingScheme::Action::ReadCheckpoint);
    }

    v << 3.0;
    mesh->data(sendDataIndex)->values() = v;
    cplScheme.addComputedTime(timeStepSize);

    cplScheme.firstSynchronization({});
    cplScheme.firstExchange();
    cplScheme.secondSynchronization();
    cplScheme.secondExchange();

    if (i < maxIterations - 1) {
      BOOST_TEST(not cplScheme.isTimeWindowComplete());
    } else {
      // window complete since max iterations reached
      BOOST_TEST(cplScheme.isTimeWindowComplete());
    }
  }

  // third window
  if (context.isNamed(first)) {
    // first order extrapolation
    BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 4); // = 2*3 - 2
  } else if (context.isNamed(second)) {
    // first order extrapolation
    BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 5); // = 2*3 - 1
  }

  // reached end of simulation, ready to finalize
  BOOST_TEST(not cplScheme.isCouplingOngoing());

  cplScheme.finalize();
}

/// Test that cplScheme gives correct results when applying extrapolation using non-zero initial data.
BOOST_AUTO_TEST_CASE(FirstOrderWithInitializationAndAcceleration)
{
  /**
   * Perform first order extrapolation and use initialization
   *
   * Do two time windows with three iterations each.
   *
   * Each participant writes dummy data to other participant, received data is checked.
   *
   **/

  PRECICE_TEST("Participant0"_on(1_rank), "Participant1"_on(1_rank), Require::Events);
  testing::ConnectionOptions options;
  options.useOnlyPrimaryCom = true;
  auto m2n                  = context.connectPrimaryRanks("Participant0", "Participant1", options);

  xml::XMLTag root = xml::getRootTag();

  // Create a data configuration, to simplify configuration of data

  mesh::PtrDataConfiguration dataConfig(new mesh::DataConfiguration(root));
  const int                  geometrical_dimensions = 3; // 3d problem
  const int                  data_dimensions        = 1; // only one sample in data
  dataConfig->setDimensions(geometrical_dimensions);
  dataConfig->addData("Data0", data_dimensions);
  dataConfig->addData("Data1", data_dimensions);

  mesh::MeshConfiguration meshConfig(root, dataConfig);
  meshConfig.setDimensions(3);
  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", geometrical_dimensions, testing::nextMeshID()));
  const auto    dataID0 = mesh->createData("Data0", data_dimensions, 0_dataID)->getID();
  const auto    dataID1 = mesh->createData("Data1", data_dimensions, 1_dataID)->getID();
  mesh->createVertex(Eigen::Vector3d::Zero());
  mesh->allocateDataValues();
  meshConfig.addMesh(mesh);

  // Create all parameters necessary to create an ImplicitCouplingScheme object
  const double maxTime            = CouplingScheme::UNDEFINED_TIME;
  const int    maxTimeWindows     = 2;
  const double timeWindowSize     = 0.1;
  const int    maxIterations      = 3;
  const int    extrapolationOrder = 1;
  const double timeStepSize       = timeWindowSize; // solver is not subcycling
  std::string  first("Participant0");
  std::string  second("Participant1");
  int          sendDataIndex        = -1;
  int          receiveDataIndex     = -1;
  int          convergenceDataIndex = -1;

  BOOST_TEST(dataID0 == 0);
  BOOST_TEST(dataID1 == 1);

  if (context.isNamed(first)) {
    sendDataIndex        = dataID0;
    receiveDataIndex     = dataID1;
    convergenceDataIndex = receiveDataIndex;
  } else {
    sendDataIndex        = dataID1;
    receiveDataIndex     = dataID0;
    convergenceDataIndex = sendDataIndex;
  }

  // Create the coupling scheme object
  cplscheme::ParallelCouplingScheme cplScheme(
      maxTime, maxTimeWindows, timeWindowSize, 16, first, second,
      context.name, m2n, constants::FIXED_TIME_WINDOW_SIZE,
      BaseCouplingScheme::Implicit, maxIterations, extrapolationOrder);
  cplScheme.addDataToSend(mesh->data(sendDataIndex), mesh, context.isNamed(second));
  cplScheme.addDataToReceive(mesh->data(receiveDataIndex), mesh, context.isNamed(first));
  cplScheme.determineInitialDataExchange();

  // Add acceleration
  acceleration::PtrAcceleration ptrAcceleration(new acceleration::ConstantRelaxationAcceleration(0.5, std::vector<int>({sendDataIndex})));
  cplScheme.setAcceleration(ptrAcceleration);

  // Add convergence measures
  const int                              minIterations = maxIterations;
  cplscheme::impl::PtrConvergenceMeasure minIterationConvMeasure1(
      new cplscheme::impl::MinIterationConvergenceMeasure(minIterations));
  cplScheme.addConvergenceMeasure(convergenceDataIndex, false, false, minIterationConvMeasure1, true);

  Eigen::VectorXd v(1); // buffer for data

  // ensure that data is uninitialized
  BOOST_TEST(mesh->data(receiveDataIndex)->values().size() == 1);
  BOOST_TEST(testing::equals(mesh->data(receiveDataIndex)->values()(0), 0.0));
  BOOST_TEST(mesh->data(sendDataIndex)->values().size() == 1);
  BOOST_TEST(testing::equals(mesh->data(sendDataIndex)->values()(0), 0.0));

  if (context.isNamed(first)) {
    BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::InitializeData));
  } else {
    BOOST_TEST(context.isNamed(second));
    BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::InitializeData));
    v << 4.0;
    mesh->data(sendDataIndex)->values() = v;
    cplScheme.markActionFulfilled(CouplingScheme::Action::InitializeData);
    BOOST_TEST(mesh->data(sendDataIndex)->values().size() == 1);
    BOOST_TEST(testing::equals(mesh->data(sendDataIndex)->values()(0), 4.0));
  }

  if (context.isNamed(first)) {
    // first participant receives initial data = 4 (see above)
    cplScheme.initialize(0.0, 1);
    BOOST_TEST(cplScheme.hasDataBeenReceived());
    BOOST_TEST(mesh->data(receiveDataIndex)->values().size() == 1);
    BOOST_TEST(testing::equals(mesh->data(receiveDataIndex)->values()(0), 4.0));
    // first participant does not send any data here
    BOOST_TEST(mesh->data(sendDataIndex)->values().size() == 1);
    BOOST_TEST(testing::equals(mesh->data(sendDataIndex)->values()(0), 0.0));
  } else {
    // second participant result no initial data = 0
    cplScheme.initialize(0.0, 1);
    BOOST_TEST(not cplScheme.hasDataBeenReceived());
    BOOST_TEST(context.isNamed(second));
    BOOST_TEST(mesh->data(receiveDataIndex)->values().size() == 1);
    BOOST_TEST(testing::equals(mesh->data(receiveDataIndex)->values()(0), 0.0));
    // second participant has send data above (should remain untouched)
    BOOST_TEST(mesh->data(sendDataIndex)->values().size() == 1);
    BOOST_TEST(testing::equals(mesh->data(sendDataIndex)->values()(0), 4.0));
  }

  // first window
  for (int i = 0; i < maxIterations; i++) {
    // first, second and third iteration
    BOOST_TEST(cplScheme.isCouplingOngoing());

    if (context.isNamed(first)) {
      if (i == 0) {
        // data is initialized for first participant
        BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 4);
      } else if (i == 1) {
        // accelerated data from second participant
        BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 3); // = 0.5 * 4 + 0.5 * 2
      } else if (i == 2) {
        // accelerated data from second participant
        BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 2.5); // = 0.5 * 3 + 0.5 * 2
      }
    } else if (context.isNamed(second)) {
      if (i == 0) {
        // data is uninitialized for second participant
        BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 0);
      } else if (i == 1) {
        // accelerated data from first participant
        BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 0.5); // = 0.5 * 0 + 0.5 * 1
      } else if (i == 2) {
        // accelerated data from first participant
        BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 0.75); // = 0.5 * 0.5 + 0.5 * 1
      }
    }

    if (i == 0) {
      BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint));
      cplScheme.markActionFulfilled(CouplingScheme::Action::WriteCheckpoint);
      BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint));
    } else {
      BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint));
      BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint));
      cplScheme.markActionFulfilled(CouplingScheme::Action::ReadCheckpoint);
    }

    // write data to mesh
    if (context.isNamed(first)) {
      v << 1.0;
    } else if (context.isNamed(second)) {
      v << 2.0;
    }
    mesh->data(sendDataIndex)->values() = v;
    cplScheme.addComputedTime(timeStepSize);

    cplScheme.firstSynchronization({});
    cplScheme.firstExchange();
    cplScheme.secondSynchronization();
    cplScheme.secondExchange();

    if (i < maxIterations - 1) {
      BOOST_TEST(not cplScheme.isTimeWindowComplete());
    } else {
      // window complete since max iterations reached
      BOOST_TEST(cplScheme.isTimeWindowComplete());
    }
  }

  // second window
  for (int i = 0; i < maxIterations; i++) {
    // first, second and third iteration
    BOOST_TEST(cplScheme.isCouplingOngoing());
    if (context.isNamed(first)) {
      if (i == 0) {
        // first order extrapolation uses initial data and final value from last window.
        BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 0); // = 2*2 - 4
      } else if (i == 1) {
        // accelerated data from second participant
        BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 1.5); // = 0.5 * 0 + 0.5 * 3
      } else if (i == 2) {
        // accelerated data from second participant
        BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 2.25); // = 0.5 * 1.5 + 0.5 * 3
      }
    } else if (context.isNamed(second)) {
      if (i == 0) {
        // first order extrapolation uses initial data and final value from last window.
        BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 2); // = 2*1 - 0
      } else if (i == 1) {
        // accelerated data from first participant
        BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 2.5); // = 0.5 * 2 + 0.5 * 3
      } else if (i == 2) {
        // accelerated data from first participant
        BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 2.75); // = 0.5 * 2.5 + 0.5 * 3
      }
    }

    if (i == 0) {
      BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint));
      cplScheme.markActionFulfilled(CouplingScheme::Action::WriteCheckpoint);
      BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint));
    } else {
      BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint));
      BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint));
      cplScheme.markActionFulfilled(CouplingScheme::Action::ReadCheckpoint);
    }

    v << 3.0;
    mesh->data(sendDataIndex)->values() = v;
    cplScheme.addComputedTime(timeStepSize);

    cplScheme.firstSynchronization({});
    cplScheme.firstExchange();
    cplScheme.secondSynchronization();
    cplScheme.secondExchange();

    if (i < maxIterations - 1) {
      BOOST_TEST(not cplScheme.isTimeWindowComplete());
    } else {
      // window complete since max iterations reached
      BOOST_TEST(cplScheme.isTimeWindowComplete());
    }
  }

  // third window
  if (context.isNamed(first)) {
    // first order extrapolation
    BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 4); // = 2*3 - 2
  } else if (context.isNamed(second)) {
    // first order extrapolation
    BOOST_TEST(mesh->data(receiveDataIndex)->values()(0) == 5); // = 2*3 - 1
  }

  // reached end of simulation, ready to finalize
  BOOST_TEST(not cplScheme.isCouplingOngoing());

  cplScheme.finalize();
}

BOOST_AUTO_TEST_SUITE_END()

#endif // not PRECICE_NO_MPI

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

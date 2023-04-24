#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include "../CompositionalCouplingScheme.hpp"
#include "../Constants.hpp"
#include "../SharedPointer.hpp"
#include "../config/CouplingSchemeConfiguration.hpp"
#include "DummyCouplingScheme.hpp"
#include "cplscheme/CouplingScheme.hpp"
#include "m2n/M2N.hpp"
#include "m2n/SharedPointer.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "precice/config/ParticipantConfiguration.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "xml/XMLTag.hpp"

using namespace precice;
using namespace precice::cplscheme;

#ifndef PRECICE_NO_MPI

BOOST_AUTO_TEST_SUITE(CplSchemeTests)

struct CompositionalCouplingSchemeFixture : m2n::WhiteboxAccessor {
  std::string _pathToTests;

  CompositionalCouplingSchemeFixture()
  {
    _pathToTests = testing::getPathToSources() + "/cplscheme/tests/";
  }

  void setupAndRunThreeSolverCoupling(const std::string &configFilename, const precice::testing::TestContext &context)
  {
    using namespace mesh;

    std::string nameParticipant0("Participant0");
    std::string nameParticipant1("Participant1");
    std::string nameParticipant2("Participant2");
    int         dimensions = 3;

    xml::XMLTag          root = xml::getRootTag();
    PtrDataConfiguration dataConfig(new DataConfiguration(root));
    dataConfig->setDimensions(dimensions);
    PtrMeshConfiguration meshConfig(new MeshConfiguration(root, dataConfig));
    meshConfig->setDimensions(dimensions);
    m2n::M2NConfiguration::SharedPointer         m2nConfig(new m2n::M2NConfiguration(root));
    precice::config::PtrParticipantConfiguration participantConfig(new precice::config::ParticipantConfiguration(root, meshConfig));
    participantConfig->setDimensions(dimensions);
    CouplingSchemeConfiguration cplSchemeConfig(root, meshConfig, m2nConfig, participantConfig);

    const xml::ConfigurationContext ccontext{context.name, 0, 1};
    xml::configure(root, ccontext, configFilename);

    // some dummy mesh
    meshConfig->meshes().at(0)->createVertex(Eigen::Vector3d(1.0, 1.0, 1.0));
    meshConfig->meshes().at(0)->createVertex(Eigen::Vector3d(2.0, 1.0, -1.0));
    meshConfig->meshes().at(0)->createVertex(Eigen::Vector3d(3.0, 1.0, 1.0));
    meshConfig->meshes().at(0)->createVertex(Eigen::Vector3d(4.0, 1.0, -1.0));

    m2n::PtrM2N m2n0 = m2nConfig->getM2N(nameParticipant0, nameParticipant1);
    m2n::PtrM2N m2n1 = m2nConfig->getM2N(nameParticipant1, nameParticipant2);

    if (context.isNamed(nameParticipant0)) {
      connect(nameParticipant0, nameParticipant1, context.name, m2n0);
    } else if (context.isNamed(nameParticipant1)) {
      connect(nameParticipant0, nameParticipant1, context.name, m2n0);
      connect(nameParticipant1, nameParticipant2, context.name, m2n1);
    } else {
      connect(nameParticipant1, nameParticipant2, context.name, m2n1);
    }

    runThreeSolverCoupling(cplSchemeConfig.getCouplingScheme(context.name),
                           context.name, meshConfig);
  }

  void runThreeSolverCoupling(
      PtrCouplingScheme          cplScheme,
      const std::string &        participantName,
      mesh::PtrMeshConfiguration meshConfig)
  {
    BOOST_TEST(meshConfig->meshes().size() == 1);
    mesh::PtrMesh mesh = meshConfig->meshes().at(0);
    BOOST_TEST(mesh->data().size() == 3);
    BOOST_TEST(mesh->vertices().size() > 0);

    double computedTime      = 0.0;
    int    computedTimesteps = 0;

    if (participantName == std::string("Participant0")) {
      cplScheme->initialize(0.0, 1);
      BOOST_TEST(not cplScheme->hasDataBeenReceived());
      BOOST_TEST(not cplScheme->isTimeWindowComplete());
      BOOST_TEST(cplScheme->isCouplingOngoing());
      while (cplScheme->isCouplingOngoing()) {
        BOOST_REQUIRE(computedTime < 1.1);
        BOOST_TEST(testing::equals(0.1, cplScheme->getNextTimeStepMaxSize()));
        if (cplScheme->isActionRequired(CouplingScheme::Action::WriteCheckpoint)) {
          cplScheme->markActionFulfilled(CouplingScheme::Action::WriteCheckpoint);
        }
        cplScheme->addComputedTime(cplScheme->getNextTimeStepMaxSize());
        cplScheme->firstSynchronization({});
        cplScheme->firstExchange();
        cplScheme->secondSynchronization();
        cplScheme->secondExchange();
        if (cplScheme->isActionRequired(CouplingScheme::Action::ReadCheckpoint)) {
          cplScheme->markActionFulfilled(CouplingScheme::Action::ReadCheckpoint);
        } else {
          BOOST_TEST(cplScheme->isTimeWindowComplete());
          computedTime += cplScheme->getNextTimeStepMaxSize();
          computedTimesteps++;
        }
        BOOST_TEST(testing::equals(computedTime, cplScheme->getTime()));
        BOOST_TEST(computedTimesteps == cplScheme->getTimeWindows() - 1);
        BOOST_TEST(cplScheme->hasDataBeenReceived());
      }
      cplScheme->finalize();
      BOOST_TEST(computedTimesteps == 10);
      BOOST_TEST(cplScheme->isTimeWindowComplete());
      BOOST_TEST(not cplScheme->isCouplingOngoing());
      BOOST_TEST(cplScheme->getNextTimeStepMaxSize() > 0.0); // ??
    } else if (participantName == std::string("Participant1")) {
      cplScheme->initialize(0.0, 1);
      BOOST_TEST(!cplScheme->hasDataBeenReceived());
      cplScheme->receiveResultOfFirstAdvance();
      BOOST_TEST(cplScheme->hasDataBeenReceived());
      BOOST_TEST(not cplScheme->isTimeWindowComplete());
      BOOST_TEST(cplScheme->isCouplingOngoing());
      while (cplScheme->isCouplingOngoing()) {
        BOOST_REQUIRE(computedTime < 1.1);
        BOOST_TEST(testing::equals(0.1, cplScheme->getNextTimeStepMaxSize()));
        if (cplScheme->isActionRequired(CouplingScheme::Action::WriteCheckpoint)) {
          cplScheme->markActionFulfilled(CouplingScheme::Action::WriteCheckpoint);
        }
        cplScheme->addComputedTime(cplScheme->getNextTimeStepMaxSize());
        cplScheme->firstSynchronization({});
        cplScheme->firstExchange();
        cplScheme->secondSynchronization();
        cplScheme->secondExchange();
        if (cplScheme->isActionRequired(CouplingScheme::Action::ReadCheckpoint)) {
          cplScheme->markActionFulfilled(CouplingScheme::Action::ReadCheckpoint);
        } else {
          BOOST_TEST(cplScheme->isTimeWindowComplete());
          computedTime += cplScheme->getNextTimeStepMaxSize();
          computedTimesteps++;
        }
        BOOST_TEST(testing::equals(computedTime, cplScheme->getTime()));
        BOOST_TEST(computedTimesteps == cplScheme->getTimeWindows() - 1);
        BOOST_TEST(cplScheme->hasDataBeenReceived());
      }
      cplScheme->finalize();
      BOOST_TEST(computedTimesteps == 10);
      BOOST_TEST(cplScheme->isTimeWindowComplete());
      BOOST_TEST(not cplScheme->isCouplingOngoing());
      BOOST_TEST(cplScheme->getNextTimeStepMaxSize() > 0.0); // ??
    } else {
      BOOST_TEST(participantName == std::string("Participant2"), participantName);
      cplScheme->initialize(0.0, 1);
      BOOST_TEST(!cplScheme->hasDataBeenReceived());
      cplScheme->receiveResultOfFirstAdvance();
      BOOST_TEST(cplScheme->hasDataBeenReceived());
      BOOST_TEST(not cplScheme->isTimeWindowComplete());
      BOOST_TEST(cplScheme->isCouplingOngoing());
      while (cplScheme->isCouplingOngoing()) {
        BOOST_REQUIRE(computedTime < 1.1);
        BOOST_TEST(testing::equals(0.1, cplScheme->getNextTimeStepMaxSize()));
        if (cplScheme->isActionRequired(CouplingScheme::Action::WriteCheckpoint)) {
          cplScheme->markActionFulfilled(CouplingScheme::Action::WriteCheckpoint);
        }
        cplScheme->addComputedTime(cplScheme->getNextTimeStepMaxSize());
        cplScheme->firstSynchronization({});
        cplScheme->firstExchange();
        cplScheme->secondSynchronization();
        cplScheme->secondExchange();
        if (cplScheme->isActionRequired(CouplingScheme::Action::ReadCheckpoint)) {
          cplScheme->markActionFulfilled(CouplingScheme::Action::ReadCheckpoint);
        } else {
          BOOST_TEST(cplScheme->isTimeWindowComplete());
          computedTime += cplScheme->getNextTimeStepMaxSize();
          computedTimesteps++;
        }
        BOOST_TEST(testing::equals(computedTime, cplScheme->getTime()));
        BOOST_TEST(computedTimesteps == cplScheme->getTimeWindows() - 1);
        if (cplScheme->isCouplingOngoing())
          BOOST_TEST(cplScheme->hasDataBeenReceived());
      }
      cplScheme->finalize();
      BOOST_TEST(computedTimesteps == 10);
      BOOST_TEST(cplScheme->isTimeWindowComplete());
      BOOST_TEST(not cplScheme->isCouplingOngoing());
      BOOST_TEST(cplScheme->getNextTimeStepMaxSize() > 0.0); // ??
    }
  }

  void connect(const std::string &participant0,
               const std::string &participant1,
               const std::string &localParticipant,
               m2n::PtrM2N        communication) const
  {
    BOOST_TEST(communication);
    BOOST_TEST(not communication->isConnected());
    useOnlyPrimaryCom(communication) = true;
    if (participant0 == localParticipant) {
      communication->requestPrimaryRankConnection(participant1, participant0);
    } else {
      BOOST_TEST(participant1 == localParticipant);
      communication->acceptPrimaryRankConnection(participant1, participant0);
    }
  }
};

BOOST_AUTO_TEST_SUITE(DummySchemeCompositionTests)

// Test two explicit dummy coupling schemes
BOOST_AUTO_TEST_CASE(testDummySchemeCompositionExplicit2)
{
  PRECICE_TEST(1_rank, Require::Events);
  int                         numberIterations = 1;
  int                         maxTimeSteps     = 10;
  PtrCouplingScheme           scheme1(new tests::DummyCouplingScheme(numberIterations, maxTimeSteps));
  PtrCouplingScheme           scheme2(new tests::DummyCouplingScheme(numberIterations, maxTimeSteps));
  CompositionalCouplingScheme composition;
  composition.addCouplingScheme(scheme1);
  composition.addCouplingScheme(scheme2);
  composition.initialize(0.0, 1);
  int advances = 0;
  while (composition.isCouplingOngoing()) {
    composition.firstSynchronization({});
    composition.firstExchange();
    composition.secondSynchronization();
    composition.secondExchange();
    advances++;
  }
  composition.finalize();
  BOOST_TEST(advances == 10);
  BOOST_TEST(scheme1->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme2->getTimeWindows() - 1 == 10);
}

// Test three explicit dummy coupling schemes
BOOST_AUTO_TEST_CASE(testDummySchemeCompositionExplicit3)
{
  PRECICE_TEST(1_rank, Require::Events);
  int               numberIterations = 1;
  int               maxTimeSteps     = 10;
  PtrCouplingScheme scheme1(
      new tests::DummyCouplingScheme(numberIterations, maxTimeSteps));
  PtrCouplingScheme scheme2(
      new tests::DummyCouplingScheme(numberIterations, maxTimeSteps));
  PtrCouplingScheme scheme3(
      new tests::DummyCouplingScheme(numberIterations, maxTimeSteps));
  CompositionalCouplingScheme composition;
  composition.addCouplingScheme(scheme1);
  composition.addCouplingScheme(scheme2);
  composition.addCouplingScheme(scheme3);
  composition.initialize(0.0, 1);
  int advances = 0;
  while (composition.isCouplingOngoing()) {
    composition.firstSynchronization({});
    composition.firstExchange();
    composition.secondSynchronization();
    composition.secondExchange();
    advances++;
  }
  composition.finalize();
  BOOST_TEST(advances == 10);
  BOOST_TEST(scheme1->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme2->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme3->getTimeWindows() - 1 == 10);
}

// Test E, I(2)
BOOST_AUTO_TEST_CASE(testDummySchemeCompositionExplicit1Implicit2)
{
  PRECICE_TEST(1_rank, Require::Events);
  int               numberIterations = 1;
  int               maxTimeSteps     = 10;
  PtrCouplingScheme scheme1(new tests::DummyCouplingScheme(numberIterations, maxTimeSteps));
  numberIterations = 2;
  PtrCouplingScheme           scheme2(new tests::DummyCouplingScheme(numberIterations, maxTimeSteps));
  CompositionalCouplingScheme composition;
  composition.addCouplingScheme(scheme1);
  composition.addCouplingScheme(scheme2);
  composition.initialize(0.0, 1);
  int advances = 0;
  BOOST_TEST_MESSAGE("Init Expl " << scheme1->getTimeWindows());
  BOOST_TEST_MESSAGE("Init Impl " << scheme2->getTimeWindows());
  while (composition.isCouplingOngoing()) {
    composition.firstSynchronization({});
    composition.firstExchange();
    composition.secondSynchronization();
    composition.secondExchange();
    advances++;
    // a 1, e 2, i 2
    BOOST_TEST_CONTEXT("Advance Nr " << advances)
    {
      BOOST_TEST_MESSAGE("Expl " << scheme1->getTimeWindows());
      BOOST_TEST_MESSAGE("Impl " << scheme2->getTimeWindows());
      if (advances % 2 == 0) {
        BOOST_TEST(scheme2->isActionRequired(CouplingScheme::Action::WriteCheckpoint));
        BOOST_TEST(scheme1->getTimeWindows() - 1 == advances / 2);
      } else {
        BOOST_TEST(scheme2->isActionRequired(CouplingScheme::Action::ReadCheckpoint));
        BOOST_TEST(scheme1->getTimeWindows() - 1 == (advances + 1) / 2);
      }
    }
  }
  composition.finalize();
  BOOST_TEST(advances == 20);
  BOOST_TEST(scheme1->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme2->getTimeWindows() - 1 == 10);
}

// Test I(2), E
BOOST_AUTO_TEST_CASE(testDummySchemeCompositionImplicit2Explicit1)
{
  PRECICE_TEST(1_rank, Require::Events);
  int               numberIterations = 2;
  int               maxTimeSteps     = 10;
  PtrCouplingScheme scheme1(new tests::DummyCouplingScheme(numberIterations, maxTimeSteps));
  numberIterations = 1;
  PtrCouplingScheme           scheme2(new tests::DummyCouplingScheme(numberIterations, maxTimeSteps));
  CompositionalCouplingScheme composition;
  composition.addCouplingScheme(scheme1);
  composition.addCouplingScheme(scheme2);
  composition.initialize(0.0, 1);
  int advances = 0;
  while (composition.isCouplingOngoing()) {
    composition.firstSynchronization({});
    composition.firstExchange();
    composition.secondSynchronization();
    composition.secondExchange();
    advances++;
    if (advances % 2 == 0) {
      BOOST_TEST(scheme1->isActionRequired(CouplingScheme::Action::WriteCheckpoint));
      BOOST_TEST(scheme1->getTimeWindows() - 1 == advances / 2);
    } else {
      BOOST_TEST(scheme1->isActionRequired(CouplingScheme::Action::ReadCheckpoint));
      BOOST_TEST(scheme1->getTimeWindows() - 1 == (advances - 1) / 2);
    }
  }
  composition.finalize();
  BOOST_TEST(advances == 20);
  BOOST_TEST(scheme1->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme2->getTimeWindows() - 1 == 10);
}

// Test E, I(3)
BOOST_AUTO_TEST_CASE(testDummySchemeCompositionExplicit1Implicit3)
{
  PRECICE_TEST(1_rank, Require::Events);
  int               numberIterations = 1;
  int               maxTimeSteps     = 10;
  PtrCouplingScheme scheme1(
      new tests::DummyCouplingScheme(numberIterations, maxTimeSteps));
  numberIterations = 3;
  PtrCouplingScheme scheme2(
      new tests::DummyCouplingScheme(numberIterations, maxTimeSteps));
  CompositionalCouplingScheme composition;
  composition.addCouplingScheme(scheme1);
  composition.addCouplingScheme(scheme2);
  composition.initialize(0.0, 1);
  int advances = 0;
  while (composition.isCouplingOngoing()) {
    composition.firstSynchronization({});
    composition.firstExchange();
    composition.secondSynchronization();
    composition.secondExchange();
    advances++;
    if (advances % 3 == 0) {
      BOOST_TEST(scheme2->isActionRequired(CouplingScheme::Action::WriteCheckpoint));
      BOOST_TEST(scheme1->getTimeWindows() - 1 == advances / 3);
    } else {
      BOOST_TEST(scheme2->isActionRequired(CouplingScheme::Action::ReadCheckpoint));
      BOOST_TEST(scheme1->getTimeWindows() - 1 == (advances + (3 - advances % 3)) / 3);
    }
  }
  composition.finalize();
  BOOST_TEST(advances == 30);
  BOOST_TEST(scheme1->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme2->getTimeWindows() - 1 == 10);
}

// Test I(3), E
BOOST_AUTO_TEST_CASE(testDummySchemeCompositionImplicit3Explicit1)
{
  PRECICE_TEST(1_rank, Require::Events);
  int               numberIterations = 3;
  int               maxTimeSteps     = 10;
  PtrCouplingScheme scheme1(new tests::DummyCouplingScheme(numberIterations, maxTimeSteps));
  numberIterations = 1;
  PtrCouplingScheme           scheme2(new tests::DummyCouplingScheme(numberIterations, maxTimeSteps));
  CompositionalCouplingScheme composition;
  composition.addCouplingScheme(scheme1);
  composition.addCouplingScheme(scheme2);
  composition.initialize(0.0, 1);
  int advances = 0;
  while (composition.isCouplingOngoing()) {
    composition.firstSynchronization({});
    composition.firstExchange();
    composition.secondSynchronization();
    composition.secondExchange();
    advances++;
    if (advances % 3 == 0) {
      BOOST_TEST(scheme1->isActionRequired(CouplingScheme::Action::WriteCheckpoint));
      BOOST_TEST(scheme1->getTimeWindows() - 1 == advances / 3);
    } else {
      BOOST_TEST(scheme1->isActionRequired(CouplingScheme::Action::ReadCheckpoint));
      BOOST_TEST(scheme1->getTimeWindows() - 1 == (advances - (advances % 3)) / 3);
    }
  }
  composition.finalize();
  BOOST_TEST(advances == 30);
  BOOST_TEST(scheme1->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme2->getTimeWindows() - 1 == 10);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(CompositionalCouplingSchemeTests, CompositionalCouplingSchemeFixture)

/// Test that runs on 3 processors.
BOOST_AUTO_TEST_CASE(testExplicitSchemeComposition1)
{
  PRECICE_TEST("Participant0"_on(1_rank), "Participant1"_on(1_rank), "Participant2"_on(1_rank), Require::Events);

  std::string configPath(_pathToTests + "multi-solver-coupling-1.xml");
  setupAndRunThreeSolverCoupling(configPath, context);
}

/// Test that runs on 3 processors.
BOOST_AUTO_TEST_CASE(testImplicitExplicitSchemeComposition)
{
  PRECICE_TEST("Participant0"_on(1_rank), "Participant1"_on(1_rank), "Participant2"_on(1_rank), Require::Events);

  std::string configPath(_pathToTests + "multi-solver-coupling-3.xml");
  setupAndRunThreeSolverCoupling(configPath, context);
}

/// Test that runs on 3 processors.
BOOST_AUTO_TEST_CASE(testExplicitImplicitSchemeComposition)
{
  PRECICE_TEST("Participant0"_on(1_rank), "Participant1"_on(1_rank), "Participant2"_on(1_rank), Require::Events);

  std::string configPath(_pathToTests + "multi-solver-coupling-4.xml");
  setupAndRunThreeSolverCoupling(configPath, context);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // not PRECICE_NO_MPI

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

    xml::XMLTag          root = xml::getRootTag();
    PtrDataConfiguration dataConfig(new DataConfiguration(root));
    dataConfig->setDimensions(3);
    PtrMeshConfiguration meshConfig(new MeshConfiguration(root, dataConfig));
    meshConfig->setDimensions(3);
    m2n::M2NConfiguration::SharedPointer m2nConfig(new m2n::M2NConfiguration(root));
    CouplingSchemeConfiguration          cplSchemeConfig(root, meshConfig, m2nConfig);

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

    std::string readIterationCheckpoint(constants::actionReadIterationCheckpoint());
    std::string writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());

    double computedTime      = 0.0;
    int    computedTimesteps = 0;

    if (participantName == std::string("Participant0")) {
      cplScheme->initialize(0.0, 1);
      BOOST_TEST(not cplScheme->hasDataBeenReceived());
      BOOST_TEST(not cplScheme->isTimeWindowComplete());
      BOOST_TEST(cplScheme->isCouplingOngoing());
      while (cplScheme->isCouplingOngoing()) {
        BOOST_TEST(testing::equals(0.1, cplScheme->getNextTimestepMaxLength()));
        if (cplScheme->isActionRequired(writeIterationCheckpoint)) {
          cplScheme->markActionFulfilled(writeIterationCheckpoint);
        }
        cplScheme->addComputedTime(cplScheme->getNextTimestepMaxLength());
        cplScheme->advance();
        if (cplScheme->isActionRequired(readIterationCheckpoint)) {
          cplScheme->markActionFulfilled(readIterationCheckpoint);
        } else {
          BOOST_TEST(cplScheme->isTimeWindowComplete());
          computedTime += cplScheme->getNextTimestepMaxLength();
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
      BOOST_TEST(cplScheme->getNextTimestepMaxLength() > 0.0); // ??
    } else if (participantName == std::string("Participant1")) {
      cplScheme->initialize(0.0, 1);
      BOOST_TEST(cplScheme->hasDataBeenReceived());
      BOOST_TEST(not cplScheme->isTimeWindowComplete());
      BOOST_TEST(cplScheme->isCouplingOngoing());
      while (cplScheme->isCouplingOngoing()) {
        BOOST_TEST(testing::equals(0.1, cplScheme->getNextTimestepMaxLength()));
        if (cplScheme->isActionRequired(writeIterationCheckpoint)) {
          cplScheme->markActionFulfilled(writeIterationCheckpoint);
        }
        cplScheme->addComputedTime(cplScheme->getNextTimestepMaxLength());
        cplScheme->advance();
        if (cplScheme->isActionRequired(readIterationCheckpoint)) {
          cplScheme->markActionFulfilled(readIterationCheckpoint);
        } else {
          BOOST_TEST(cplScheme->isTimeWindowComplete());
          computedTime += cplScheme->getNextTimestepMaxLength();
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
      BOOST_TEST(cplScheme->getNextTimestepMaxLength() > 0.0); // ??
    } else {
      BOOST_TEST(participantName == std::string("Participant2"), participantName);
      cplScheme->initialize(0.0, 1);
      BOOST_TEST(cplScheme->hasDataBeenReceived());
      BOOST_TEST(not cplScheme->isTimeWindowComplete());
      BOOST_TEST(cplScheme->isCouplingOngoing());
      while (cplScheme->isCouplingOngoing()) {
        BOOST_TEST(testing::equals(0.1, cplScheme->getNextTimestepMaxLength()));
        if (cplScheme->isActionRequired(writeIterationCheckpoint)) {
          cplScheme->markActionFulfilled(writeIterationCheckpoint);
        }
        cplScheme->addComputedTime(cplScheme->getNextTimestepMaxLength());
        cplScheme->advance();
        if (cplScheme->isActionRequired(readIterationCheckpoint)) {
          cplScheme->markActionFulfilled(readIterationCheckpoint);
        } else {
          BOOST_TEST(cplScheme->isTimeWindowComplete());
          computedTime += cplScheme->getNextTimestepMaxLength();
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
      BOOST_TEST(cplScheme->getNextTimestepMaxLength() > 0.0); // ??
    }
  }

  void connect(const std::string &participant0,
               const std::string &participant1,
               const std::string &localParticipant,
               m2n::PtrM2N        communication) const
  {
    BOOST_TEST(communication);
    BOOST_TEST(not communication->isConnected());
    useOnlyMasterCom(communication) = true;
    if (participant0 == localParticipant) {
      communication->requestMasterConnection(participant1, participant0);
    } else {
      BOOST_TEST(participant1 == localParticipant);
      communication->acceptMasterConnection(participant1, participant0);
    }
  }
};

BOOST_AUTO_TEST_SUITE(DummySchemeCompositionTests)

/// Test one explicit dummy coupling scheme
BOOST_AUTO_TEST_CASE(testDummySchemeCompositionExplicit1)
{
  PRECICE_TEST(1_rank, Require::Events);
  std::string writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterationCheckpoint(constants::actionReadIterationCheckpoint());

  int                         numberIterations = 1;
  int                         maxTimesteps     = 10;
  PtrCouplingScheme           scheme(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  CompositionalCouplingScheme composition;
  composition.addCouplingScheme(scheme);
  composition.initialize(0.0, 1);
  int advances = 0;
  while (composition.isCouplingOngoing()) {
    composition.advance();
    advances++;
  }
  composition.finalize();
  BOOST_TEST(advances == 10);
  BOOST_TEST(scheme->getTimeWindows() - 1 == 10);
}

// Test one implicit dummy coupling scheme
BOOST_AUTO_TEST_CASE(testDummySchemeCompositionImplicit1)
{
  PRECICE_TEST(1_rank, Require::Events);
  std::string                 writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string                 readIterationCheckpoint(constants::actionReadIterationCheckpoint());
  int                         numberIterations = 2;
  int                         maxTimesteps     = 10;
  PtrCouplingScheme           scheme(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  CompositionalCouplingScheme composition;
  composition.addCouplingScheme(scheme);
  composition.initialize(0.0, 1);
  int advances = 0;
  while (composition.isCouplingOngoing()) {
    composition.advance();
    advances++;
  }
  composition.finalize();
  BOOST_TEST(advances == 20);
  BOOST_TEST(scheme->getTimeWindows() - 1 == 10);
}

// Test two explicit dummy coupling schemes
BOOST_AUTO_TEST_CASE(testDummySchemeCompositionExplicit2)
{
  PRECICE_TEST(1_rank, Require::Events);
  std::string writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterationCheckpoint(constants::actionReadIterationCheckpoint());

  int               numberIterations = 1;
  int               maxTimesteps     = 10;
  PtrCouplingScheme scheme1(
      new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  PtrCouplingScheme scheme2(
      new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  CompositionalCouplingScheme composition;
  composition.addCouplingScheme(scheme1);
  composition.addCouplingScheme(scheme2);
  composition.initialize(0.0, 1);
  int advances = 0;
  while (composition.isCouplingOngoing()) {
    composition.advance();
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
  std::string writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterationCheckpoint(constants::actionReadIterationCheckpoint());

  int               numberIterations = 1;
  int               maxTimesteps     = 10;
  PtrCouplingScheme scheme1(
      new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  PtrCouplingScheme scheme2(
      new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  PtrCouplingScheme scheme3(
      new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  CompositionalCouplingScheme composition;
  composition.addCouplingScheme(scheme1);
  composition.addCouplingScheme(scheme2);
  composition.addCouplingScheme(scheme3);
  composition.initialize(0.0, 1);
  int advances = 0;
  while (composition.isCouplingOngoing()) {
    composition.advance();
    advances++;
  }
  composition.finalize();
  BOOST_TEST(advances == 10);
  BOOST_TEST(scheme1->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme2->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme3->getTimeWindows() - 1 == 10);
}

// Test two implicit dummy coupling schemes
BOOST_AUTO_TEST_CASE(testDummySchemeCompositionImplicit2)
{
  PRECICE_TEST(1_rank, Require::Events);
  std::string writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterationCheckpoint(constants::actionReadIterationCheckpoint());

  int               numberIterations = 2;
  int               maxTimesteps     = 10;
  PtrCouplingScheme scheme1(
      new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  PtrCouplingScheme scheme2(
      new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  CompositionalCouplingScheme composition;
  composition.addCouplingScheme(scheme1);
  composition.addCouplingScheme(scheme2);
  composition.initialize(0.0, 1);
  int advances = 0;
  while (composition.isCouplingOngoing()) {
    composition.advance();
    advances++;
    if (advances % 2 == 1) {
      BOOST_TEST(scheme1->isActionRequired(readIterationCheckpoint));
      BOOST_TEST(scheme2->isActionRequired(readIterationCheckpoint));
    } else if (advances % 2 == 0) {
      BOOST_TEST(scheme1->isActionRequired(writeIterationCheckpoint));
      BOOST_TEST(scheme2->isActionRequired(writeIterationCheckpoint));
    }
  }
  composition.finalize();
  BOOST_TEST(advances == 20);
  BOOST_TEST(scheme1->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme2->getTimeWindows() - 1 == 10);
}

// Test two implicit dummy coupling schemes with different iteration number
BOOST_AUTO_TEST_CASE(testDummySchemeCompositionImplicit2DiffIteration)
{
  PRECICE_TEST(1_rank, Require::Events);
  std::string writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterationCheckpoint(constants::actionReadIterationCheckpoint());

  int               numberIterations = 2;
  int               maxTimesteps     = 10;
  PtrCouplingScheme scheme1(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  numberIterations = 3;
  PtrCouplingScheme           scheme2(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  CompositionalCouplingScheme composition;
  composition.addCouplingScheme(scheme1);
  composition.addCouplingScheme(scheme2);
  composition.initialize(0.0, 1);
  int advances = 0;
  while (composition.isCouplingOngoing()) {
    composition.advance();
    advances++;
    if (advances % 3 == 2) {
      BOOST_TEST(scheme1->isActionRequired(writeIterationCheckpoint));
      BOOST_TEST(scheme2->isActionRequired(readIterationCheckpoint));
    } else if (advances % 3 == 1) {
      BOOST_TEST(scheme1->isActionRequired(readIterationCheckpoint));
      BOOST_TEST(scheme2->isActionRequired(readIterationCheckpoint));
    } else if (advances % 3 == 0) {
      BOOST_TEST(scheme1->isActionRequired(writeIterationCheckpoint));
      BOOST_TEST(scheme2->isActionRequired(writeIterationCheckpoint));
    }
  }
  composition.finalize();
  BOOST_TEST(advances == 30);
  BOOST_TEST(scheme1->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme2->getTimeWindows() - 1 == 10);
}

// Test three implicit dummy coupling schemes
BOOST_AUTO_TEST_CASE(testDummySchemeCompositionImplicit3)
{
  PRECICE_TEST(1_rank, Require::Events);
  std::string writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterationCheckpoint(constants::actionReadIterationCheckpoint());

  int                         numberIterations = 2;
  int                         maxTimesteps     = 10;
  PtrCouplingScheme           scheme1(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  PtrCouplingScheme           scheme2(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  PtrCouplingScheme           scheme3(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  CompositionalCouplingScheme composition;
  composition.addCouplingScheme(scheme1);
  composition.addCouplingScheme(scheme2);
  composition.addCouplingScheme(scheme3);
  composition.initialize(0.0, 1);
  int advances = 0;
  while (composition.isCouplingOngoing()) {
    composition.advance();
    advances++;
    if (advances % 2 == 0) {
      BOOST_TEST(scheme1->isActionRequired(writeIterationCheckpoint));
      BOOST_TEST(scheme2->isActionRequired(writeIterationCheckpoint));
      BOOST_TEST(scheme3->isActionRequired(writeIterationCheckpoint));
    } else {
      BOOST_TEST(scheme1->isActionRequired(readIterationCheckpoint));
      BOOST_TEST(scheme2->isActionRequired(readIterationCheckpoint));
      BOOST_TEST(scheme3->isActionRequired(readIterationCheckpoint));
    }
  }
  composition.finalize();
  BOOST_TEST(advances == 20);
  BOOST_TEST(scheme1->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme2->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme3->getTimeWindows() - 1 == 10);
}

// Test three implicit dummy coupling schemes with different iteration number
BOOST_AUTO_TEST_CASE(testDummySchemeCompositionImplicit3DiffIteration)
{
  PRECICE_TEST(1_rank, Require::Events);
  std::string writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterationCheckpoint(constants::actionReadIterationCheckpoint());

  int               numberIterations = 3;
  int               maxTimesteps     = 10;
  PtrCouplingScheme scheme1(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  numberIterations = 4;
  PtrCouplingScheme scheme2(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  numberIterations = 2;
  PtrCouplingScheme           scheme3(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  CompositionalCouplingScheme composition;
  composition.addCouplingScheme(scheme1);
  composition.addCouplingScheme(scheme2);
  composition.addCouplingScheme(scheme3);
  composition.initialize(0.0, 1);
  int advances = 0;
  while (composition.isCouplingOngoing()) {
    composition.advance();
    advances++;
    if (advances % 4 == 0) {
      BOOST_TEST(scheme1->isActionRequired(writeIterationCheckpoint));
      BOOST_TEST(scheme2->isActionRequired(writeIterationCheckpoint));
      BOOST_TEST(scheme3->isActionRequired(writeIterationCheckpoint));
    } else if (advances % 4 == 1) {
      BOOST_TEST(scheme1->isActionRequired(readIterationCheckpoint));
      BOOST_TEST(scheme2->isActionRequired(readIterationCheckpoint));
      BOOST_TEST(scheme3->isActionRequired(readIterationCheckpoint));
    } else if (advances % 4 == 2) {
      BOOST_TEST(scheme1->isActionRequired(readIterationCheckpoint));
      BOOST_TEST(scheme2->isActionRequired(readIterationCheckpoint));
      BOOST_TEST(scheme3->isActionRequired(writeIterationCheckpoint));
    } else if (advances % 4 == 3) {
      BOOST_TEST(scheme1->isActionRequired(writeIterationCheckpoint));
      BOOST_TEST(scheme2->isActionRequired(readIterationCheckpoint));
      BOOST_TEST(scheme3->isActionRequired(writeIterationCheckpoint));
    }
  }
  composition.finalize();
  BOOST_TEST(advances == 40);
  BOOST_TEST(scheme1->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme2->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme3->getTimeWindows() - 1 == 10);
}

// Test E, I(2)
BOOST_AUTO_TEST_CASE(testDummySchemeCompositionExplicit1Implicit2)
{
  PRECICE_TEST(1_rank, Require::Events);
  std::string writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterationCheckpoint(constants::actionReadIterationCheckpoint());

  int               numberIterations = 1;
  int               maxTimesteps     = 10;
  PtrCouplingScheme scheme1(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  numberIterations = 2;
  PtrCouplingScheme           scheme2(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  CompositionalCouplingScheme composition;
  composition.addCouplingScheme(scheme1);
  composition.addCouplingScheme(scheme2);
  composition.initialize(0.0, 1);
  int advances = 0;
  while (composition.isCouplingOngoing()) {
    composition.advance();
    advances++;
    if (advances % 2 == 0) {
      BOOST_TEST(scheme2->isActionRequired(writeIterationCheckpoint));
      BOOST_TEST(scheme1->getTimeWindows() - 1 == advances / 2);
    } else {
      BOOST_TEST(scheme2->isActionRequired(readIterationCheckpoint));
      BOOST_TEST(scheme1->getTimeWindows() - 1 == (advances + 1) / 2);
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
  std::string writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterationCheckpoint(constants::actionReadIterationCheckpoint());

  int               numberIterations = 2;
  int               maxTimesteps     = 10;
  PtrCouplingScheme scheme1(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  numberIterations = 1;
  PtrCouplingScheme           scheme2(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  CompositionalCouplingScheme composition;
  composition.addCouplingScheme(scheme1);
  composition.addCouplingScheme(scheme2);
  composition.initialize(0.0, 1);
  int advances = 0;
  while (composition.isCouplingOngoing()) {
    composition.advance();
    advances++;
    if (advances % 2 == 0) {
      BOOST_TEST(scheme1->isActionRequired(writeIterationCheckpoint));
      BOOST_TEST(scheme1->getTimeWindows() - 1 == advances / 2);
    } else {
      BOOST_TEST(scheme1->isActionRequired(readIterationCheckpoint));
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
  std::string writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterationCheckpoint(constants::actionReadIterationCheckpoint());

  int               numberIterations = 1;
  int               maxTimesteps     = 10;
  PtrCouplingScheme scheme1(
      new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  numberIterations = 3;
  PtrCouplingScheme scheme2(
      new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  CompositionalCouplingScheme composition;
  composition.addCouplingScheme(scheme1);
  composition.addCouplingScheme(scheme2);
  composition.initialize(0.0, 1);
  int advances = 0;
  while (composition.isCouplingOngoing()) {
    composition.advance();
    advances++;
    if (advances % 3 == 0) {
      BOOST_TEST(scheme2->isActionRequired(writeIterationCheckpoint));
      BOOST_TEST(scheme1->getTimeWindows() - 1 == advances / 3);
    } else {
      BOOST_TEST(scheme2->isActionRequired(readIterationCheckpoint));
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
  std::string writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterationCheckpoint(constants::actionReadIterationCheckpoint());

  int               numberIterations = 3;
  int               maxTimesteps     = 10;
  PtrCouplingScheme scheme1(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  numberIterations = 1;
  PtrCouplingScheme           scheme2(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  CompositionalCouplingScheme composition;
  composition.addCouplingScheme(scheme1);
  composition.addCouplingScheme(scheme2);
  composition.initialize(0.0, 1);
  int advances = 0;
  while (composition.isCouplingOngoing()) {
    composition.advance();
    advances++;
    if (advances % 3 == 0) {
      BOOST_TEST(scheme1->isActionRequired(writeIterationCheckpoint));
      BOOST_TEST(scheme1->getTimeWindows() - 1 == advances / 3);
    } else {
      BOOST_TEST(scheme1->isActionRequired(readIterationCheckpoint));
      BOOST_TEST(scheme1->getTimeWindows() - 1 == (advances - (advances % 3)) / 3);
    }
  }
  composition.finalize();
  BOOST_TEST(advances == 30);
  BOOST_TEST(scheme1->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme2->getTimeWindows() - 1 == 10);
}

// Test E, I(2), I(2)
BOOST_AUTO_TEST_CASE(testDummySchemeCompositionExplicit1Implicit2Implicit2)
{
  PRECICE_TEST(1_rank, Require::Events);
  std::string writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterationCheckpoint(constants::actionReadIterationCheckpoint());

  int               numberIterations = 1;
  int               maxTimesteps     = 10;
  PtrCouplingScheme scheme1(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  numberIterations = 2;
  PtrCouplingScheme           scheme2(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  PtrCouplingScheme           scheme3(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  CompositionalCouplingScheme composition;
  composition.addCouplingScheme(scheme1);
  composition.addCouplingScheme(scheme2);
  composition.addCouplingScheme(scheme3);
  composition.initialize(0.0, 1);
  int advances = 0;
  while (composition.isCouplingOngoing()) {
    composition.advance();
    advances++;
    if (advances % 2 == 0) {
      BOOST_TEST(scheme2->isActionRequired(writeIterationCheckpoint));
      BOOST_TEST(scheme3->isActionRequired(writeIterationCheckpoint));
      BOOST_TEST(scheme1->getTimeWindows() - 1 == advances / 2);
    } else {
      BOOST_TEST(scheme2->isActionRequired(readIterationCheckpoint));
      BOOST_TEST(scheme3->isActionRequired(readIterationCheckpoint));
      BOOST_TEST(scheme1->getTimeWindows() - 1 == (advances + 1) / 2);
    }
  }
  composition.finalize();
  BOOST_TEST(advances == 20);
  BOOST_TEST(scheme1->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme2->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme3->getTimeWindows() - 1 == 10);
}

// Test E, I(2), I(3)
BOOST_AUTO_TEST_CASE(testDummySchemeCompositionExplicit1Implicit2Implicit3)
{
  PRECICE_TEST(1_rank, Require::Events);
  std::string writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterationCheckpoint(constants::actionReadIterationCheckpoint());

  int               numberIterations = 1;
  int               maxTimesteps     = 10;
  PtrCouplingScheme scheme1(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  numberIterations = 2;
  PtrCouplingScheme scheme2(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  numberIterations = 3;
  PtrCouplingScheme           scheme3(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  CompositionalCouplingScheme composition;
  composition.addCouplingScheme(scheme1);
  composition.addCouplingScheme(scheme2);
  composition.addCouplingScheme(scheme3);
  composition.initialize(0.0, 1);
  int advances = 0;
  while (composition.isCouplingOngoing()) {
    composition.advance();
    advances++;
    if (advances % 3 == 0) {
      BOOST_TEST(scheme1->getTimeWindows() - 1 == advances / 3);
      BOOST_TEST(scheme2->isActionRequired(writeIterationCheckpoint));
      BOOST_TEST(scheme3->isActionRequired(writeIterationCheckpoint));
    } else if (advances % 3 == 1) {
      BOOST_TEST(scheme2->isActionRequired(readIterationCheckpoint));
      BOOST_TEST(scheme3->isActionRequired(readIterationCheckpoint));
      BOOST_TEST(scheme1->getTimeWindows() - 1 == (advances + 2) / 3);
    } else if (advances % 3 == 2) {
      BOOST_TEST(scheme2->isActionRequired(writeIterationCheckpoint));
      BOOST_TEST(scheme3->isActionRequired(readIterationCheckpoint));
      BOOST_TEST(scheme1->getTimeWindows() - 1 == (advances + 1) / 3);
    }
  }
  composition.finalize();
  BOOST_TEST(advances == 30);
  BOOST_TEST(scheme1->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme2->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme3->getTimeWindows() - 1 == 10);
}

// Test I(2), I(2), E
BOOST_AUTO_TEST_CASE(testDummySchemeCompositionImplicit2Implicit2Explicit1)
{
  PRECICE_TEST(1_rank, Require::Events);
  std::string writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterationCheckpoint(constants::actionReadIterationCheckpoint());

  int               numberIterations = 2;
  int               maxTimesteps     = 10;
  PtrCouplingScheme scheme1(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  PtrCouplingScheme scheme2(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  numberIterations = 1;
  PtrCouplingScheme           scheme3(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  CompositionalCouplingScheme composition;
  composition.addCouplingScheme(scheme1);
  composition.addCouplingScheme(scheme2);
  composition.addCouplingScheme(scheme3);
  composition.initialize(0.0, 1);
  int advances = 0;
  while (composition.isCouplingOngoing()) {
    composition.advance();
    advances++;
    if (advances % 2 == 0) {
      BOOST_TEST(scheme1->isActionRequired(writeIterationCheckpoint));
      BOOST_TEST(scheme2->isActionRequired(writeIterationCheckpoint));
      BOOST_TEST(scheme1->getTimeWindows() - 1 == advances / 2);
      BOOST_TEST(scheme2->getTimeWindows() - 1 == advances / 2);
      BOOST_TEST(scheme3->getTimeWindows() - 1 == advances / 2);
    } else {
      BOOST_TEST(scheme1->isActionRequired(readIterationCheckpoint));
      BOOST_TEST(scheme2->isActionRequired(readIterationCheckpoint));
      BOOST_TEST(scheme1->getTimeWindows() - 1 == (advances - 1) / 2);
      BOOST_TEST(scheme2->getTimeWindows() - 1 == (advances - 1) / 2);
      BOOST_TEST(scheme3->getTimeWindows() - 1 == (advances - 1) / 2);
    }
  }
  composition.finalize();
  BOOST_TEST(advances == 20);
  BOOST_TEST(scheme1->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme2->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme3->getTimeWindows() - 1 == 10);
}

// Test I(2), I(2), E
BOOST_AUTO_TEST_CASE(testDummySchemeCompositionImplicit2Implicit2Explicit1DiffIterations)
{
  PRECICE_TEST(1_rank, Require::Events);
  std::string writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterationCheckpoint(constants::actionReadIterationCheckpoint());

  int               numberIterations = 3;
  int               maxTimesteps     = 10;
  PtrCouplingScheme scheme1(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  numberIterations = 2;
  PtrCouplingScheme scheme2(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  numberIterations = 1;
  PtrCouplingScheme           scheme3(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  CompositionalCouplingScheme composition;
  composition.addCouplingScheme(scheme1);
  composition.addCouplingScheme(scheme2);
  composition.addCouplingScheme(scheme3);
  composition.initialize(0.0, 1);
  int advances = 0;
  while (composition.isCouplingOngoing()) {
    composition.advance();
    advances++;
    if (advances % 3 == 0) {
      BOOST_TEST(scheme1->isActionRequired(writeIterationCheckpoint));
      BOOST_TEST(scheme2->isActionRequired(writeIterationCheckpoint));
      BOOST_TEST(scheme1->getTimeWindows() - 1 == advances / 3);
      BOOST_TEST(scheme2->getTimeWindows() - 1 == advances / 3);
      BOOST_TEST(scheme3->getTimeWindows() - 1 == advances / 3);
    } else if (advances % 3 == 1) {
      BOOST_TEST(scheme1->isActionRequired(readIterationCheckpoint));
      BOOST_TEST(scheme2->isActionRequired(readIterationCheckpoint));
      BOOST_TEST(scheme1->getTimeWindows() - 1 == (advances - 1) / 3);
      BOOST_TEST(scheme2->getTimeWindows() - 1 == (advances - 1) / 3);
      BOOST_TEST(scheme3->getTimeWindows() - 1 == (advances - 1) / 3);
    } else if (advances % 3 == 2) {
      BOOST_TEST(scheme1->isActionRequired(readIterationCheckpoint));
      BOOST_TEST(scheme2->isActionRequired(writeIterationCheckpoint));
      BOOST_TEST(scheme1->getTimeWindows() - 1 == (advances - 2) / 3);
      BOOST_TEST(scheme2->getTimeWindows() - 1 == (advances + 1) / 3);
      BOOST_TEST(scheme3->getTimeWindows() - 1 == (advances - 2) / 3);
    }
  }
  composition.finalize();
  BOOST_TEST(advances == 30);
  BOOST_TEST(scheme1->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme2->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme3->getTimeWindows() - 1 == 10);
}

BOOST_AUTO_TEST_CASE(testDummySchemeCompositionUntitled) /// @todo give a better name, what is this test doing?
{
  PRECICE_TEST(1_rank, Require::Events);
  std::string writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterationCheckpoint(constants::actionReadIterationCheckpoint());

  int               numberIterations = 3;
  int               maxTimesteps     = 10;
  PtrCouplingScheme scheme1(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  numberIterations = 1;
  PtrCouplingScheme scheme2(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  numberIterations = 2;
  PtrCouplingScheme           scheme3(new tests::DummyCouplingScheme(numberIterations, maxTimesteps));
  CompositionalCouplingScheme composition;
  composition.addCouplingScheme(scheme1);
  composition.addCouplingScheme(scheme2);
  composition.addCouplingScheme(scheme3);
  composition.initialize(0.0, 1);
  int advances = 0;
  while (composition.isCouplingOngoing()) {
    composition.advance();
    advances++;
    if (advances % 4 >= 3) {
      BOOST_TEST(scheme1->isActionRequired(writeIterationCheckpoint));
      BOOST_TEST(scheme1->getTimeWindows() - 1 == (advances - (advances % 4) + 4) / 4);
    } else if (advances % 4 != 0) {
      BOOST_TEST(scheme1->isActionRequired(readIterationCheckpoint));
    }
    BOOST_TEST(scheme2->getTimeWindows() - 1 == (advances + 1) / 4);
    BOOST_TEST(scheme2->getTimeWindows() - 1 == (advances + 1) / 4);
  }
  composition.finalize();
  BOOST_TEST(advances == 40);
  BOOST_TEST(scheme1->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme2->getTimeWindows() - 1 == 10);
  BOOST_TEST(scheme3->getTimeWindows() - 1 == 10);
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
BOOST_AUTO_TEST_CASE(testImplicitSchemeComposition)
{
  PRECICE_TEST("Participant0"_on(1_rank), "Participant1"_on(1_rank), "Participant2"_on(1_rank), Require::Events);

  std::string configPath(_pathToTests + "multi-solver-coupling-2.xml");
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

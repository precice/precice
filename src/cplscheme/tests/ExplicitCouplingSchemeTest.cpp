#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include "com/SharedPointer.hpp"
#include "cplscheme/BaseCouplingScheme.hpp"
#include "cplscheme/Constants.hpp"
#include "cplscheme/CouplingScheme.hpp"
#include "cplscheme/SerialCouplingScheme.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"
#include "m2n/M2N.hpp"
#include "m2n/SharedPointer.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "xml/XMLTag.hpp"

using namespace precice;
using namespace precice::cplscheme;

#ifndef PRECICE_NO_MPI

BOOST_AUTO_TEST_SUITE(CplSchemeTests)

void runSimpleExplicitCoupling(
    CouplingScheme &               cplScheme,
    const std::string &            participantName,
    const mesh::MeshConfiguration &meshConfig)
{
  BOOST_TEST(meshConfig.meshes().size() == 1);
  mesh::PtrMesh mesh = meshConfig.meshes().at(0);
  BOOST_TEST(mesh->data().size() == 2);
  auto &dataValues0 = mesh->data().at(0)->values();
  auto &dataValues1 = mesh->data().at(1)->values();
  BOOST_TEST(mesh->vertices().size() > 0);
  mesh::Vertex &  vertex     = mesh->vertices().at(0);
  double          valueData0 = 1.0;
  Eigen::VectorXd valueData1 = Eigen::VectorXd::Constant(3, 1.0);

  double computedTime      = 0.0;
  int    computedTimesteps = 0;

  if (participantName == std::string("Participant0")) {
    cplScheme.initialize(0.0, 1);
    BOOST_TEST(not cplScheme.hasDataBeenReceived());
    BOOST_TEST(not cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));
    BOOST_TEST(not cplScheme.isActionRequired(constants::actionReadIterationCheckpoint()));
    BOOST_TEST(not cplScheme.isTimeWindowComplete());
    BOOST_TEST(cplScheme.isCouplingOngoing());
    while (cplScheme.isCouplingOngoing()) {
      dataValues0(vertex.getID()) = valueData0;
      computedTime += cplScheme.getNextTimestepMaxLength();
      computedTimesteps++;
      cplScheme.addComputedTime(cplScheme.getNextTimestepMaxLength());
      cplScheme.advance();
      BOOST_TEST(cplScheme.isTimeWindowComplete());
      BOOST_TEST(testing::equals(computedTime, cplScheme.getTime()));
      BOOST_TEST(computedTimesteps == cplScheme.getTimeWindows() - 1);
      BOOST_TEST(not cplScheme.isActionRequired("WriteIterationCheckpoint"));
      BOOST_TEST(not cplScheme.isActionRequired("ReadIterationCheckpoint"));
      BOOST_TEST(cplScheme.isTimeWindowComplete());
      if (cplScheme.isCouplingOngoing()) {
        // No receive takes place for the participant that has started the
        // coupled simulation, in the last advance call
        Eigen::VectorXd value = dataValues1.segment(vertex.getID() * 3, 3);
        BOOST_TEST(testing::equals(value, valueData1));
      }
      BOOST_TEST(cplScheme.hasDataBeenReceived());
      // Increment data values, to test if send/receive operations are also
      // correct in following timesteps.
      valueData0 += 1.0;
      valueData1 += Eigen::VectorXd::Constant(3, 1.0);
    }
    cplScheme.finalize();
    // Validate results
    BOOST_TEST(testing::equals(computedTime, 1.0));
    BOOST_TEST(computedTimesteps == 10);
    BOOST_TEST(not cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"));
    BOOST_TEST(not cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"));
    BOOST_TEST(cplScheme.isTimeWindowComplete());
    BOOST_TEST(not cplScheme.isCouplingOngoing());
    BOOST_TEST(cplScheme.getNextTimestepMaxLength() > 0.0);
  } else if (participantName == std::string("Participant1")) {
    cplScheme.initialize(0.0, 1);
    BOOST_TEST(cplScheme.hasDataBeenReceived());
    double value = dataValues0(vertex.getID());
    BOOST_TEST(testing::equals(value, valueData0));
    valueData0 += 1.0;
    BOOST_TEST(not cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"));
    BOOST_TEST(not cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"));
    BOOST_TEST(not cplScheme.isTimeWindowComplete());
    BOOST_TEST(cplScheme.isCouplingOngoing());
    while (cplScheme.isCouplingOngoing()) {
      dataValues1.segment(vertex.getID() * 3, 3) = valueData1;
      computedTime += cplScheme.getNextTimestepMaxLength();
      computedTimesteps++;
      cplScheme.addComputedTime(cplScheme.getNextTimestepMaxLength());
      cplScheme.advance();
      BOOST_TEST(testing::equals(computedTime, cplScheme.getTime()));
      BOOST_TEST(computedTimesteps == cplScheme.getTimeWindows() - 1);
      BOOST_TEST(not cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"));
      BOOST_TEST(not cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"));
      BOOST_TEST(cplScheme.isTimeWindowComplete());
      if (cplScheme.isCouplingOngoing()) {
        // The participant not starting the coupled simulation does neither
        // receive nor send data in the last call to advance
        BOOST_TEST(cplScheme.hasDataBeenReceived());
        double value = dataValues0(vertex.getID());
        BOOST_TEST(testing::equals(value, valueData0));
      }
      valueData0 += 1.0;
      valueData1 += Eigen::VectorXd::Constant(3, 1.0);
    }
    cplScheme.finalize();
    // Validate results
    BOOST_TEST(testing::equals(computedTime, 1.0));
    BOOST_TEST(computedTimesteps == 10);
    BOOST_TEST(not cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"));
    BOOST_TEST(not cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"));
    BOOST_TEST(cplScheme.isTimeWindowComplete());
    BOOST_TEST(not cplScheme.isCouplingOngoing());
    BOOST_TEST(cplScheme.getNextTimestepMaxLength() > 0.0);
  }
}

void runExplicitCouplingWithSubcycling(
    CouplingScheme &               cplScheme,
    const std::string &            participantName,
    const mesh::MeshConfiguration &meshConfig)
{
  BOOST_TEST(meshConfig.meshes().size() == 1);
  mesh::PtrMesh mesh = meshConfig.meshes().at(0);
  BOOST_TEST(mesh->data().size() == 2);
  BOOST_TEST(mesh->vertices().size() > 0);
  mesh::Vertex &  vertex      = mesh->vertices().at(0);
  double          valueData0  = 1.0;
  Eigen::VectorXd valueData1  = Eigen::VectorXd::Constant(3, 1.0);
  auto &          dataValues0 = mesh->data().at(0)->values();
  auto &          dataValues1 = mesh->data().at(1)->values();

  double      computedTime      = 0.0;
  int         computedTimesteps = 0;
  std::string nameParticipant0("Participant0");
  std::string nameParticipant1("Participant1");
  BOOST_TEST(((participantName == nameParticipant0) || (participantName == nameParticipant1)));
  if (participantName == nameParticipant0) {
    cplScheme.initialize(0.0, 1);
    double dtDesired = cplScheme.getNextTimestepMaxLength() / 2.0;
    double dtUsed    = dtDesired;
    BOOST_TEST(not cplScheme.hasDataBeenReceived());
    BOOST_TEST(not cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"));
    BOOST_TEST(not cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"));
    BOOST_TEST(not cplScheme.isTimeWindowComplete());
    BOOST_TEST(cplScheme.isCouplingOngoing());
    while (cplScheme.isCouplingOngoing()) {
      dataValues0(vertex.getID()) = valueData0;
      computedTime += dtUsed;
      computedTimesteps++;
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
      if (computedTimesteps % 2 == 0) {
        // Data exchange takes only place at every second local timestep,
        // since a subcycling of 2 is used.
        BOOST_TEST(cplScheme.isTimeWindowComplete());
        if (cplScheme.isCouplingOngoing()) {
          // No receive takes place for the participant that has started the
          // coupled simulation, in the last advance call.
          Eigen::VectorXd value = dataValues1.segment(vertex.getID() * 3, 3);
          BOOST_TEST(testing::equals(value, valueData1));
        }
        BOOST_TEST(cplScheme.hasDataBeenReceived());
        // Increment data values, to test if send/receive operations are also
        // correct in following timesteps.
        valueData0 += 1.0;
        valueData1 += Eigen::VectorXd::Constant(3, 1.0);
      } else {
        BOOST_TEST(not cplScheme.isTimeWindowComplete());
      }
    }
    cplScheme.finalize();
    BOOST_TEST(testing::equals(computedTime, 1.0));
    BOOST_TEST(computedTimesteps == 20);
    BOOST_TEST(not cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"));
    BOOST_TEST(not cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"));
    BOOST_TEST(cplScheme.isTimeWindowComplete());
    BOOST_TEST(not cplScheme.isCouplingOngoing());
    BOOST_TEST(cplScheme.getNextTimestepMaxLength() > 0.0);
  } else if (participantName == nameParticipant1) {
    // Start coupling
    cplScheme.initialize(0.0, 1);
    // Validate current coupling status
    BOOST_TEST(cplScheme.hasDataBeenReceived());
    BOOST_TEST(testing::equals(dataValues0(vertex.getID()), valueData0));
    valueData0 += 1.0;
    BOOST_TEST(not cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"));
    BOOST_TEST(not cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"));
    BOOST_TEST(not cplScheme.isTimeWindowComplete());
    BOOST_TEST(cplScheme.isCouplingOngoing());
    while (cplScheme.isCouplingOngoing()) {
      dataValues1.segment(vertex.getID() * 3, 3) = valueData1;
      computedTime += cplScheme.getNextTimestepMaxLength();
      computedTimesteps++;
      cplScheme.addComputedTime(cplScheme.getNextTimestepMaxLength());
      cplScheme.advance();
      BOOST_TEST(testing::equals(computedTime, cplScheme.getTime()));
      BOOST_TEST(computedTimesteps == cplScheme.getTimeWindows() - 1);
      BOOST_TEST(not cplScheme.isActionRequired("constants::actionWriteIterationCheckpoint()"));
      BOOST_TEST(not cplScheme.isActionRequired("constants::actionReadIterationCheckpoint()"));
      BOOST_TEST(cplScheme.isTimeWindowComplete());
      if (cplScheme.isCouplingOngoing()) {
        // The participant not starting the coupled simulation does neither
        // receive nor send data in the last call to advance
        BOOST_TEST(cplScheme.hasDataBeenReceived());
        BOOST_TEST(testing::equals(dataValues0(vertex.getID()), valueData0));
        BOOST_TEST(cplScheme.hasDataBeenReceived());
      }
      valueData0 += 1.0;
      valueData1 += Eigen::VectorXd::Constant(3, 1.0);
    }
    cplScheme.finalize();
    BOOST_TEST(testing::equals(computedTime, 1.0));
    BOOST_TEST(computedTimesteps == 10);
    BOOST_TEST(not cplScheme.isActionRequired(constants::actionWriteIterationCheckpoint()));
    BOOST_TEST(not cplScheme.isActionRequired(constants::actionReadIterationCheckpoint()));
    BOOST_TEST(cplScheme.isTimeWindowComplete());
    BOOST_TEST(not cplScheme.isCouplingOngoing());
    BOOST_TEST(cplScheme.getNextTimestepMaxLength() > 0.0);
  }
}

struct ExplicitCouplingSchemeFixture : m2n::WhiteboxAccessor {
  std::string _pathToTests;

  ExplicitCouplingSchemeFixture()
  {
    _pathToTests = testing::getPathToSources() + "/cplscheme/tests/";
  }

  void connect(
      const std::string &participant0,
      const std::string &participant1,
      const std::string &localParticipant,
      m2n::PtrM2N &      communication)
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

BOOST_FIXTURE_TEST_SUITE(ExplicitCouplingSchemeTests, ExplicitCouplingSchemeFixture)

/// Test that runs on 2 processors.

BOOST_AUTO_TEST_CASE(testSimpleExplicitCoupling)
{
  PRECICE_TEST("Participant0"_on(1_rank), "Participant1"_on(1_rank), Require::Events);
  testing::ConnectionOptions options;
  options.useOnlyMasterCom = true;
  auto m2n                 = context.connectMasters("Participant0", "Participant1", options);

  xml::XMLTag                root = xml::getRootTag();
  mesh::PtrDataConfiguration dataConfig(new mesh::DataConfiguration(root));
  dataConfig->addData("Data0", 1);
  dataConfig->addData("Data1", 3);
  mesh::MeshConfiguration meshConfig(root, dataConfig);
  mesh::PtrMesh           mesh(new mesh::Mesh("Mesh", 3, false, testing::nextMeshID()));
  mesh->createData("Data0", 1);
  mesh->createData("Data1", 3);
  mesh->createVertex(Eigen::Vector3d::Zero());
  mesh->allocateDataValues();
  meshConfig.addMesh(mesh);

  double      maxTime        = 1.0;
  int         maxTimesteps   = 10;
  double      timestepLength = 0.1;
  std::string nameParticipant0("Participant0");
  std::string nameParticipant1("Participant1");
  int         sendDataIndex    = -1;
  int         receiveDataIndex = -1;

  if (context.isNamed(nameParticipant0)) {
    sendDataIndex    = 0;
    receiveDataIndex = 1;
  } else {
    sendDataIndex    = 1;
    receiveDataIndex = 0;
  }
  cplscheme::SerialCouplingScheme cplScheme(
      maxTime, maxTimesteps, timestepLength, 12, nameParticipant0,
      nameParticipant1, context.name, m2n, constants::FIXED_TIME_WINDOW_SIZE,
      BaseCouplingScheme::Explicit);
  cplScheme.addDataToSend(mesh->data().at(sendDataIndex), mesh, false);
  cplScheme.addDataToReceive(mesh->data().at(receiveDataIndex), mesh, false);
  runSimpleExplicitCoupling(cplScheme, context.name, meshConfig);
}

/// Test that runs on 2 processors.
BOOST_AUTO_TEST_CASE(testConfiguredSimpleExplicitCoupling)
{
  PRECICE_TEST("Participant0"_on(1_rank), "Participant1"_on(1_rank), Require::Events);

  using namespace mesh;

  std::string configurationPath(_pathToTests + "explicit-coupling-scheme-1.xml");
  std::string nameParticipant0("Participant0");
  std::string nameParticipant1("Participant1");

  xml::XMLTag          root = xml::getRootTag();
  PtrDataConfiguration dataConfig(new DataConfiguration(root));
  dataConfig->setDimensions(3);
  PtrMeshConfiguration meshConfig(new MeshConfiguration(root, dataConfig));
  meshConfig->setDimensions(3);
  m2n::M2NConfiguration::SharedPointer m2nConfig(new m2n::M2NConfiguration(root));
  CouplingSchemeConfiguration          cplSchemeConfig(root, meshConfig, m2nConfig);

  xml::ConfigurationContext ccontext{context.name, 0, 1};
  xml::configure(root, ccontext, configurationPath);
  m2n::PtrM2N m2n = m2nConfig->getM2N(nameParticipant0, nameParticipant1);

  // some dummy mesh
  meshConfig->meshes().at(0)->createVertex(Eigen::Vector3d(1.0, 1.0, 1.0));
  meshConfig->meshes().at(0)->createVertex(Eigen::Vector3d(2.0, 1.0, -1.0));
  meshConfig->meshes().at(0)->createVertex(Eigen::Vector3d(3.0, 1.0, 1.0));
  meshConfig->meshes().at(0)->createVertex(Eigen::Vector3d(4.0, 1.0, -1.0));
  meshConfig->meshes().at(0)->allocateDataValues();

  connect(nameParticipant0, nameParticipant1, context.name, m2n);
  runSimpleExplicitCoupling(*cplSchemeConfig.getCouplingScheme(context.name),
                            context.name, *meshConfig);
}

/// Test that runs on 2 processors.
BOOST_AUTO_TEST_CASE(testExplicitCouplingFirstParticipantSetsDt)
{
  PRECICE_TEST("Participant0"_on(1_rank), "Participant1"_on(1_rank), Require::Events);

  using namespace mesh;
  std::string configurationPath(_pathToTests + "explicit-coupling-scheme-2.xml");
  std::string nameParticipant0("Participant0");
  std::string nameParticipant1("Participant1");

  xml::XMLTag          root = xml::getRootTag();
  PtrDataConfiguration dataConfig(new DataConfiguration(root));
  dataConfig->setDimensions(3);
  PtrMeshConfiguration meshConfig(new MeshConfiguration(root, dataConfig));
  meshConfig->setDimensions(3);
  m2n::M2NConfiguration::SharedPointer m2nConfig(new m2n::M2NConfiguration(root));
  CouplingSchemeConfiguration          cplSchemeConfig(root, meshConfig, m2nConfig);

  xml::ConfigurationContext ccontext{context.name, 0, 1};
  xml::configure(root, ccontext, configurationPath);
  m2n::PtrM2N m2n = m2nConfig->getM2N(nameParticipant0, nameParticipant1);

  // some dummy mesh
  meshConfig->meshes().at(0)->createVertex(Eigen::Vector3d(1.0, 1.0, 1.0));
  meshConfig->meshes().at(0)->createVertex(Eigen::Vector3d(2.0, 1.0, -1.0));
  meshConfig->meshes().at(0)->createVertex(Eigen::Vector3d(3.0, 1.0, 1.0));
  meshConfig->meshes().at(0)->createVertex(Eigen::Vector3d(4.0, 1.0, -1.0));
  meshConfig->meshes().at(0)->allocateDataValues();

  connect(nameParticipant0, nameParticipant1, context.name, m2n);
  CouplingScheme &cplScheme = *cplSchemeConfig.getCouplingScheme(context.name);

  double computedTime      = 0.0;
  int    computedTimesteps = 0;
  if (context.isNamed(nameParticipant0)) {
    double dt = 0.3;
    cplScheme.initialize(0.0, 1);
    BOOST_TEST(not cplScheme.isTimeWindowComplete());
    BOOST_TEST(cplScheme.isCouplingOngoing());
    while (cplScheme.isCouplingOngoing()) {
      computedTime += dt;
      computedTimesteps++;
      cplScheme.addComputedTime(dt);
      cplScheme.advance();
      BOOST_TEST(cplScheme.isTimeWindowComplete());
      BOOST_TEST(testing::equals(computedTime, cplScheme.getTime()));
      BOOST_TEST(computedTimesteps == cplScheme.getTimeWindows() - 1);
      if (cplScheme.isCouplingOngoing()) {
        BOOST_TEST(cplScheme.hasDataBeenReceived());
      }
    }
    cplScheme.finalize();
    BOOST_TEST(testing::equals(computedTime, 1.2));
    BOOST_TEST(computedTimesteps == 4);
    BOOST_TEST(cplScheme.isTimeWindowComplete());
    BOOST_TEST(not cplScheme.isCouplingOngoing());
  } else {
    BOOST_TEST(context.isNamed(nameParticipant1));
    cplScheme.initialize(0.0, 1);
    BOOST_TEST(not cplScheme.isTimeWindowComplete());
    BOOST_TEST(cplScheme.isCouplingOngoing());
    while (cplScheme.isCouplingOngoing()) {
      computedTime += cplScheme.getTimeWindowSize();
      computedTimesteps++;
      cplScheme.addComputedTime(cplScheme.getTimeWindowSize());
      cplScheme.advance();
      BOOST_TEST(cplScheme.isTimeWindowComplete());
      BOOST_TEST(testing::equals(computedTime, cplScheme.getTime()));
      BOOST_TEST(computedTimesteps == cplScheme.getTimeWindows() - 1);
      if (cplScheme.isCouplingOngoing()) {
        BOOST_TEST(cplScheme.hasDataBeenReceived());
      }
    }
    cplScheme.finalize();
    BOOST_TEST(testing::equals(computedTime, 1.2));
    BOOST_TEST(computedTimesteps == 4);
    BOOST_TEST(cplScheme.isTimeWindowComplete());
    BOOST_TEST(not cplScheme.isCouplingOngoing());
  }
}

/// Test that runs on 2 processors.
BOOST_AUTO_TEST_CASE(testSerialDataInitialization)
{
  PRECICE_TEST("Participant0"_on(1_rank), "Participant1"_on(1_rank), Require::Events);

  using namespace mesh;

  std::string configurationPath(_pathToTests + "serial-explicit-coupling-datainit.xml");
  std::string nameParticipant0("Participant0");
  std::string nameParticipant1("Participant1");

  xml::XMLTag          root = xml::getRootTag();
  PtrDataConfiguration dataConfig(new DataConfiguration(root));
  dataConfig->setDimensions(2);
  PtrMeshConfiguration meshConfig(new MeshConfiguration(root, dataConfig));
  meshConfig->setDimensions(2);
  m2n::M2NConfiguration::SharedPointer m2nConfig(new m2n::M2NConfiguration(root));
  CouplingSchemeConfiguration          cplSchemeConfig(root, meshConfig, m2nConfig);

  xml::ConfigurationContext ccontext{context.name, 0, 1};
  xml::configure(root, ccontext, configurationPath);
  m2n::PtrM2N m2n = m2nConfig->getM2N(nameParticipant0, nameParticipant1);

  // some dummy mesh
  meshConfig->meshes().at(0)->createVertex(Eigen::Vector2d(1.0, 1.0));
  meshConfig->meshes().at(0)->createVertex(Eigen::Vector2d(2.0, -1.0));
  meshConfig->meshes().at(0)->createVertex(Eigen::Vector2d(3.0, 1.0));
  meshConfig->meshes().at(0)->createVertex(Eigen::Vector2d(4.0, -1.0));
  meshConfig->meshes().at(0)->allocateDataValues();

  connect(nameParticipant0, nameParticipant1, context.name, m2n);
  CouplingScheme &cplScheme = *cplSchemeConfig.getCouplingScheme(context.name);

  BOOST_TEST(meshConfig->meshes().size() == 1);
  mesh::PtrMesh mesh = meshConfig->meshes().at(0);
  BOOST_TEST(mesh->data().size() == 3);
  auto &dataValues0 = mesh->data().at(0)->values();
  auto &dataValues1 = mesh->data().at(1)->values();
  auto &dataValues2 = mesh->data().at(2)->values();

  if (context.isNamed(nameParticipant0)) {
    cplScheme.initialize(0.0, 1);
    BOOST_TEST(not cplScheme.isActionRequired(constants::actionWriteInitialData()));
    cplScheme.initializeData();
    BOOST_TEST(cplScheme.hasDataBeenReceived());
    BOOST_TEST(testing::equals(dataValues0(0), 0.0));
    BOOST_TEST(testing::equals(dataValues1(0), 1.0));
    dataValues2(0) = 2.0;
    cplScheme.addComputedTime(cplScheme.getNextTimestepMaxLength());
    cplScheme.advance();
    BOOST_TEST(not cplScheme.isCouplingOngoing());
    cplScheme.finalize();
  } else {
    BOOST_TEST(context.isNamed(nameParticipant1));
    cplScheme.initialize(0.0, 1);
    BOOST_TEST(not cplScheme.hasDataBeenReceived());
    BOOST_TEST(cplScheme.isActionRequired(constants::actionWriteInitialData()));
    dataValues1(0) = 1.0;
    cplScheme.markActionFulfilled(constants::actionWriteInitialData());
    cplScheme.initializeData();
    BOOST_TEST(testing::equals(dataValues2(0), 2.0));
    cplScheme.addComputedTime(cplScheme.getNextTimestepMaxLength());
    cplScheme.advance();
    BOOST_TEST(not cplScheme.isCouplingOngoing());
    cplScheme.finalize();
  }
}

/// Test that runs on 2 processors.
BOOST_AUTO_TEST_CASE(testParallelDataInitialization)
{
  PRECICE_TEST("Participant0"_on(1_rank), "Participant1"_on(1_rank), Require::Events);

  using namespace mesh;

  std::string configurationPath(_pathToTests + "parallel-explicit-coupling-datainit.xml");
  std::string nameParticipant0("Participant0");
  std::string nameParticipant1("Participant1");

  xml::XMLTag          root = xml::getRootTag();
  PtrDataConfiguration dataConfig(new DataConfiguration(root));
  dataConfig->setDimensions(2);
  PtrMeshConfiguration meshConfig(new MeshConfiguration(root, dataConfig));
  meshConfig->setDimensions(2);
  m2n::M2NConfiguration::SharedPointer m2nConfig(new m2n::M2NConfiguration(root));
  CouplingSchemeConfiguration          cplSchemeConfig(root, meshConfig, m2nConfig);

  xml::ConfigurationContext ccontext{context.name, 0, 1};
  xml::configure(root, ccontext, configurationPath);
  m2n::PtrM2N m2n = m2nConfig->getM2N(nameParticipant0, nameParticipant1);

  // some dummy mesh
  meshConfig->meshes().at(0)->createVertex(Eigen::Vector2d(1.0, 1.0));
  meshConfig->meshes().at(0)->createVertex(Eigen::Vector2d(2.0, -1.0));
  meshConfig->meshes().at(0)->createVertex(Eigen::Vector2d(3.0, 1.0));
  meshConfig->meshes().at(0)->createVertex(Eigen::Vector2d(4.0, -1.0));
  meshConfig->meshes().at(0)->allocateDataValues();

  connect(nameParticipant0, nameParticipant1, context.name, m2n);
  CouplingScheme &cplScheme = *cplSchemeConfig.getCouplingScheme(context.name);

  BOOST_TEST(meshConfig->meshes().size() == 1);
  mesh::PtrMesh mesh = meshConfig->meshes().at(0);
  BOOST_TEST(mesh->data().size() == 3);
  auto &dataValues0 = mesh->data().at(0)->values();
  auto &dataValues1 = mesh->data().at(1)->values();
  auto &dataValues2 = mesh->data().at(2)->values();

  if (context.isNamed(nameParticipant0)) {
    cplScheme.initialize(0.0, 1);
    BOOST_TEST(cplScheme.isActionRequired(constants::actionWriteInitialData()));
    dataValues2(0) = 3.0;
    cplScheme.markActionFulfilled(constants::actionWriteInitialData());
    cplScheme.initializeData();
    BOOST_TEST(cplScheme.hasDataBeenReceived());
    BOOST_TEST(testing::equals(dataValues0(0), 0.0));
    BOOST_TEST(testing::equals(dataValues1(0), 1.0));
    dataValues2(0) = 2.0;
    cplScheme.addComputedTime(cplScheme.getNextTimestepMaxLength());
    cplScheme.advance();
    BOOST_TEST(testing::equals(dataValues0(0), 4.0));
    BOOST_TEST(not cplScheme.isCouplingOngoing());
    cplScheme.finalize();
  } else {
    BOOST_TEST(context.isNamed(nameParticipant1));
    cplScheme.initialize(0.0, 1);
    BOOST_TEST(not cplScheme.hasDataBeenReceived());
    BOOST_TEST(cplScheme.isActionRequired(constants::actionWriteInitialData()));
    dataValues1(0) = 1.0;
    cplScheme.markActionFulfilled(constants::actionWriteInitialData());
    cplScheme.initializeData();
    BOOST_TEST(cplScheme.hasDataBeenReceived());
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
BOOST_AUTO_TEST_CASE(testExplicitCouplingWithSubcycling)
{
  PRECICE_TEST("Participant0"_on(1_rank), "Participant1"_on(1_rank), Require::Events);
  testing::ConnectionOptions options;
  options.useOnlyMasterCom = true;
  auto m2n                 = context.connectMasters("Participant0", "Participant1", options);

  xml::XMLTag                root = xml::getRootTag();
  mesh::PtrDataConfiguration dataConfig(new mesh::DataConfiguration(root));
  dataConfig->setDimensions(3);
  dataConfig->addData("Data0", 1);
  dataConfig->addData("Data1", 3);
  mesh::MeshConfiguration meshConfig(root, dataConfig);
  meshConfig.setDimensions(3);
  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 3, false, testing::nextMeshID()));
  mesh->createData("Data0", 1);
  mesh->createData("Data1", 3);
  mesh->createVertex(Eigen::Vector3d::Zero());
  mesh->allocateDataValues();
  meshConfig.addMesh(mesh);

  double      maxTime        = 1.0;
  int         maxTimesteps   = 10;
  double      timestepLength = 0.1;
  std::string nameParticipant0("Participant0");
  std::string nameParticipant1("Participant1");
  int         sendDataIndex    = -1;
  int         receiveDataIndex = -1;
  if (context.isNamed(nameParticipant0)) {
    sendDataIndex    = 0;
    receiveDataIndex = 1;
  } else {
    sendDataIndex    = 1;
    receiveDataIndex = 0;
  }
  cplscheme::SerialCouplingScheme cplScheme(
      maxTime, maxTimesteps, timestepLength, 12, nameParticipant0,
      nameParticipant1, context.name, m2n, constants::FIXED_TIME_WINDOW_SIZE,
      BaseCouplingScheme::Explicit);
  cplScheme.addDataToSend(mesh->data().at(sendDataIndex), mesh, false);
  cplScheme.addDataToReceive(mesh->data().at(receiveDataIndex), mesh, false);
  runExplicitCouplingWithSubcycling(cplScheme, context.name, meshConfig);
}

/// Test that runs on 2 processors.
BOOST_AUTO_TEST_CASE(testConfiguredExplicitCouplingWithSubcycling)
{
  PRECICE_TEST("Participant0"_on(1_rank), "Participant1"_on(1_rank), Require::Events);

  using namespace mesh;

  std::string configurationPath(_pathToTests + "explicit-coupling-scheme-1.xml");
  std::string nameParticipant0("Participant0");
  std::string nameParticipant1("Participant1");

  xml::XMLTag          root = xml::getRootTag();
  PtrDataConfiguration dataConfig(new DataConfiguration(root));
  dataConfig->setDimensions(3);
  PtrMeshConfiguration meshConfig(new MeshConfiguration(root, dataConfig));
  meshConfig->setDimensions(3);
  m2n::M2NConfiguration::SharedPointer m2nConfig(new m2n::M2NConfiguration(root));
  CouplingSchemeConfiguration          cplSchemeConfig(root, meshConfig, m2nConfig);

  xml::ConfigurationContext ccontext{context.name, 0, 1};
  xml::configure(root, ccontext, configurationPath);
  m2n::PtrM2N m2n = m2nConfig->getM2N(nameParticipant0, nameParticipant1);
  // some dummy mesh
  meshConfig->meshes().at(0)->createVertex(Eigen::Vector3d(1.0, 1.0, 1.0));
  meshConfig->meshes().at(0)->createVertex(Eigen::Vector3d(2.0, -1.0, 1.0));
  meshConfig->meshes().at(0)->createVertex(Eigen::Vector3d(3.0, 1.0, 1.0));
  meshConfig->meshes().at(0)->createVertex(Eigen::Vector3d(4.0, -1.0, 1.0));
  meshConfig->meshes().at(0)->allocateDataValues();

  connect(nameParticipant0, nameParticipant1, context.name, m2n);
  runExplicitCouplingWithSubcycling(
      *cplSchemeConfig.getCouplingScheme(context.name), context.name,
      *meshConfig);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // not PRECICE_NO_MPI

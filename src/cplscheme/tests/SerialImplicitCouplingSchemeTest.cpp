#include <Eigen/Core>
#include <algorithm>
#include <boost/test/tools/old/interface.hpp>
#include <iterator>
#include <memory>
#include <string>
#include <vector>
#include "acceleration/ConstantRelaxationAcceleration.hpp"
#include "acceleration/SharedPointer.hpp"
#include "acceleration/config/AccelerationConfiguration.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "com/SharedPointer.hpp"
#include "cplscheme/BaseCouplingScheme.hpp"
#include "cplscheme/Constants.hpp"
#include "cplscheme/CouplingData.hpp"
#include "cplscheme/CouplingScheme.hpp"
#include "cplscheme/SerialCouplingScheme.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"
#include "cplscheme/impl/AbsoluteConvergenceMeasure.hpp"
#include "cplscheme/impl/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/DistributedComFactory.hpp"
#include "m2n/M2N.hpp"
#include "m2n/SharedPointer.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "math/differences.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "precice/config/ParticipantConfiguration.hpp"
#include "testing/SerialCouplingSchemeFixture.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "xml/XMLTag.hpp"

using namespace precice;
using namespace precice::cplscheme;

#ifndef PRECICE_NO_MPI

BOOST_AUTO_TEST_SUITE(CplSchemeTests)

void runCoupling(
    CouplingScheme                &cplScheme,
    const std::string             &nameParticipant,
    const mesh::MeshConfiguration &meshConfig,
    const std::vector<int>        &validIterations)
{
  BOOST_REQUIRE(meshConfig.meshes().size() == 1);
  mesh::PtrMesh mesh = meshConfig.meshes().at(0);
  BOOST_REQUIRE(mesh->data().size() == 2);
  BOOST_REQUIRE(!mesh->empty());
  BOOST_REQUIRE(!validIterations.empty());

  mesh::Vertex   &vertex               = mesh->vertex(0);
  int             index                = vertex.getID();
  auto           &dataValues0          = mesh->data(0)->values();
  auto           &dataValues1          = mesh->data(1)->values();
  double          initialStepsizeData0 = 5.0;
  double          stepsizeData0        = 5.0;
  Eigen::VectorXd initialStepsizeData1 = Eigen::VectorXd::Constant(3, 5.0);
  Eigen::VectorXd stepsizeData1        = Eigen::VectorXd::Constant(3, 5.0);
  double          computedTime         = 0.0;
  int             computedTimesteps    = 0;
  std::string     nameParticipant0("Participant0");
  std::string     nameParticipant1("Participant1");
  BOOST_TEST(((nameParticipant == nameParticipant0) || (nameParticipant == nameParticipant1)));
  int                              iterationCount      = 0;
  std::vector<int>::const_iterator iterValidIterations = validIterations.begin();

  if (nameParticipant == nameParticipant0) {
    mesh->data(0)->setSampleAtTime(0, time::Sample{mesh->data(0)->getDimensions(), mesh->data(0)->values()});
    cplScheme.initialize();
    BOOST_TEST(not cplScheme.isTimeWindowComplete());
    BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint));
    BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint));
    BOOST_TEST(not cplScheme.hasDataBeenReceived());

    // Tells coupling scheme, that a checkpoint has been created.
    // All required actions have to be performed before calling advance().
    cplScheme.markActionFulfilled(CouplingScheme::Action::WriteCheckpoint);
    BOOST_TEST(cplScheme.isActionFulfilled(CouplingScheme::Action::WriteCheckpoint));

    while (cplScheme.isCouplingOngoing()) {
      dataValues0(index) += stepsizeData0;
      // The max time step size is required to be obeyed.
      double maxTimeStepSize = cplScheme.getNextTimeStepMaxSize();
      cplScheme.addComputedTime(maxTimeStepSize);
      mesh->data(0)->setSampleAtTime(cplScheme.getTime(), time::Sample{mesh->data(0)->getDimensions(), mesh->data(0)->values()});
      cplScheme.firstSynchronization({});
      cplScheme.firstExchange();
      cplScheme.secondSynchronization();
      cplScheme.secondExchange();
      iterationCount++;
      // A coupling timestep is complete, when the coupling iterations are
      // globally converged and if subcycling steps have filled one global
      // timestep.
      if (cplScheme.isTimeWindowComplete()) {
        // Advance participant time and timestep
        computedTime += maxTimeStepSize;
        computedTimesteps++;
        BOOST_TEST(testing::equals(computedTime, cplScheme.getTime()));
        BOOST_TEST(computedTimesteps == cplScheme.getTimeWindows() - 1);
        // The iteration number is enforced by the controlled decrease of the
        // change of data written
        BOOST_TEST(iterationCount == *iterValidIterations);
        if (cplScheme.isCouplingOngoing()) {
          BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint));
          cplScheme.markActionFulfilled(CouplingScheme::Action::WriteCheckpoint);
          BOOST_TEST(cplScheme.isActionFulfilled(CouplingScheme::Action::WriteCheckpoint));
        } else {
          BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint));
          BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint));
        }
        iterationCount = 0;
        iterValidIterations++;
        if (iterValidIterations == validIterations.end()) {
          BOOST_REQUIRE(not cplScheme.isCouplingOngoing());
        }
        // Reset data values, to simulate same convergence behavior of
        // interface values in next timestep.
        stepsizeData0 = initialStepsizeData0;
      } else { // coupling timestep is not yet complete
        BOOST_TEST(cplScheme.isCouplingOngoing());
        BOOST_TEST(iterationCount < *iterValidIterations);
        BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint));
        cplScheme.markActionFulfilled(CouplingScheme::Action::ReadCheckpoint);
        BOOST_TEST(cplScheme.isActionFulfilled(CouplingScheme::Action::ReadCheckpoint));
        // The written data value is decreased in a regular manner, in order
        // to achieve a predictable convergence.
        stepsizeData0 -= 1.0;
      }
      // the first participant always receives new data
      // if(cplScheme.isCouplingOngoing())
      BOOST_TEST(cplScheme.hasDataBeenReceived());
    }
    cplScheme.finalize(); // Ends the coupling scheme
    BOOST_TEST(testing::equals(computedTime, 0.3));
    BOOST_TEST(computedTimesteps == 3);
  } else if (nameParticipant == nameParticipant1) {
    mesh->data(1)->setSampleAtTime(0, time::Sample{mesh->data(1)->getDimensions(), mesh->data(1)->values()});
    cplScheme.initialize();
    BOOST_TEST(not cplScheme.isTimeWindowComplete());
    BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint));
    BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint));
    BOOST_TEST(cplScheme.hasDataBeenReceived());

    // Tells coupling scheme, that a checkpoint has been created.
    // All required actions have to be performed before calling advance().
    cplScheme.markActionFulfilled(CouplingScheme::Action::WriteCheckpoint);
    BOOST_TEST(cplScheme.isActionFulfilled(CouplingScheme::Action::WriteCheckpoint));

    while (cplScheme.isCouplingOngoing()) {
      Eigen::VectorXd currentData(3);
      currentData = dataValues1.segment(index * 3, 3);
      currentData += stepsizeData1;
      dataValues1.segment(index * 3, 3) = currentData;
      // The max time step size is required to be obeyed.
      double maxTimeStepSize = cplScheme.getNextTimeStepMaxSize();
      cplScheme.addComputedTime(maxTimeStepSize);
      mesh->data(1)->setSampleAtTime(cplScheme.getTime(), time::Sample{mesh->data(1)->getDimensions(), mesh->data(1)->values()});
      cplScheme.firstSynchronization({});
      cplScheme.firstExchange();
      cplScheme.secondSynchronization();
      cplScheme.secondExchange();
      iterationCount++;
      // A coupling timestep is complete, when the coupling iterations are
      // globally converged and if subcycling steps have filled one global
      // timestep.
      if (cplScheme.isTimeWindowComplete()) {
        // Advance participant time and timestep
        computedTime += maxTimeStepSize;
        computedTimesteps++;
        BOOST_TEST(testing::equals(computedTime, cplScheme.getTime()));
        BOOST_TEST(computedTimesteps == cplScheme.getTimeWindows() - 1);
        // The iterations are enforced by the controlled decrease of the
        // change of data written
        BOOST_TEST(iterationCount == *iterValidIterations);
        if (cplScheme.isCouplingOngoing()) {
          BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint));
          cplScheme.markActionFulfilled(CouplingScheme::Action::WriteCheckpoint);
          BOOST_TEST(cplScheme.isActionFulfilled(CouplingScheme::Action::WriteCheckpoint));
        } else {
          BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint));
          BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint));
        }
        iterationCount = 0;
        iterValidIterations++;
        if (iterValidIterations == validIterations.end()) {
          BOOST_REQUIRE(not cplScheme.isCouplingOngoing());
        }
        // Reset data values, to simulate same convergence behavior of
        // interface values in next timestep.
        stepsizeData1 = initialStepsizeData1;
      } else { // coupling timestep is not yet complete
        BOOST_TEST(cplScheme.isCouplingOngoing());
        BOOST_TEST(iterationCount < *iterValidIterations);
        BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint));
        // The load checkpoint action requires to fallback to the cplScheme of the
        // first implicit iteration of the current timestep/time.
        cplScheme.markActionFulfilled(CouplingScheme::Action::ReadCheckpoint);
        BOOST_TEST(cplScheme.isActionFulfilled(CouplingScheme::Action::ReadCheckpoint));
        // The written data value is decreased in a regular manner, in order
        // to achieve a predictable convergence.
        // stepsizeData1 -= 1.0;
        stepsizeData1 -= Eigen::Vector3d::Constant(1.0);
      }
      // only check if data is received
      if (cplScheme.isCouplingOngoing())
        BOOST_TEST(cplScheme.hasDataBeenReceived());
    }
    cplScheme.finalize(); // Ends the coupling scheme
    BOOST_TEST(testing::equals(computedTime, 0.3));
    BOOST_TEST(computedTimesteps == 3);
  }
}

void runCouplingWithSubcycling(
    CouplingScheme                &cplScheme,
    const std::string             &nameParticipant,
    const mesh::MeshConfiguration &meshConfig,
    const std::vector<int>        &validIterations)
{
  BOOST_REQUIRE(meshConfig.meshes().size() == 1);
  mesh::PtrMesh mesh = meshConfig.meshes().at(0);
  BOOST_REQUIRE(mesh->data().size() == 2);
  BOOST_REQUIRE(!mesh->empty());
  BOOST_REQUIRE(!validIterations.empty());

  double          initialStepsizeData0 = 5.0;
  double          stepsizeData0        = 5.0;
  Eigen::Vector3d initialStepsizeData1 = Eigen::Vector3d::Constant(5.0);
  Eigen::Vector3d stepsizeData1        = Eigen::Vector3d::Constant(5.0);
  double          computedTime         = 0.0;
  int             computedTimesteps    = 0;
  std::string     nameParticipant0("Participant0");
  std::string     nameParticipant1("Participant1");
  BOOST_TEST(((nameParticipant == nameParticipant0) || (nameParticipant == nameParticipant1)));
  int                              iterationCount = 0;
  std::vector<int>::const_iterator iterValidIterations =
      validIterations.begin();

  if (nameParticipant == nameParticipant0) {
    iterationCount++; // different handling due to subcycling
    mesh->data(0)->setSampleAtTime(0, time::Sample{mesh->data(0)->getDimensions(), mesh->data(0)->values()});
    cplScheme.initialize();
    BOOST_TEST(not cplScheme.isTimeWindowComplete());
    BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint));
    BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint));
    BOOST_TEST(not cplScheme.hasDataBeenReceived());

    // Tells coupling scheme, that a checkpoint has been created.
    // All required actions have to be performed before calling advance().
    cplScheme.markActionFulfilled(CouplingScheme::Action::WriteCheckpoint);
    BOOST_TEST(cplScheme.isActionFulfilled(CouplingScheme::Action::WriteCheckpoint));

    double maxTimeStepSize      = cplScheme.getNextTimeStepMaxSize();
    double computedTimeStepSize = maxTimeStepSize / 2.0;
    int    subcyclingStep       = 0;

    // Clear data for iteration.
    mesh->data(0)->timeStepsStorage().trim();

    // Main coupling loop
    while (cplScheme.isCouplingOngoing()) {
      cplScheme.addComputedTime(computedTimeStepSize);
      mesh->data(0)->setSampleAtTime(cplScheme.getTime(), time::Sample{mesh->data(0)->getDimensions(), mesh->data(0)->values()});
      cplScheme.firstSynchronization({});
      cplScheme.firstExchange();
      cplScheme.secondSynchronization();
      cplScheme.secondExchange();
      // A coupling timestep is complete, when the coupling iterations are
      // globally converged and if subcycling steps have filled one global
      // timestep.
      if (cplScheme.isTimeWindowComplete()) {
        // Advance participant time and timestep
        mesh->data(0)->timeStepsStorage().trim();
        computedTime += maxTimeStepSize;
        computedTimesteps++;
        BOOST_TEST(testing::equals(computedTime, cplScheme.getTime()));
        BOOST_TEST(computedTimesteps == cplScheme.getTimeWindows() - 1);
        // The iteration number is enforced by the controlled decrease of the
        // change of data written
        BOOST_TEST(iterationCount == *iterValidIterations);
        if (cplScheme.isCouplingOngoing()) {
          BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint));
          cplScheme.markActionFulfilled(CouplingScheme::Action::WriteCheckpoint);
          BOOST_TEST(cplScheme.isActionFulfilled(CouplingScheme::Action::WriteCheckpoint));
        } else {
          BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint));
          BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint));
        }
        iterationCount = 1;
        iterValidIterations++;
        if (iterValidIterations == validIterations.end()) {
          BOOST_REQUIRE(not cplScheme.isCouplingOngoing());
        }
        // Reset data values, to simulate same convergence behavior of
        // interface values in next time step.
        stepsizeData0 = initialStepsizeData0;
        BOOST_TEST(subcyclingStep == 1);
        subcyclingStep = 0;
      } else { // coupling timestep is not yet complete
        BOOST_TEST(cplScheme.isCouplingOngoing());
        // If end of time window is reached
        if (cplScheme.hasDataBeenReceived()) {
          BOOST_TEST(iterationCount <= *iterValidIterations);
          BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint));
          BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint));
          cplScheme.markActionFulfilled(CouplingScheme::Action::ReadCheckpoint);
          BOOST_TEST(cplScheme.isActionFulfilled(CouplingScheme::Action::ReadCheckpoint));
          // The written data value is decreased in a regular manner, in order
          // to achieve a predictable convergence.
          stepsizeData0 -= 1.0;
          subcyclingStep = 0; // Subcycling steps
          iterationCount++;   // Implicit coupling iterations
        } else {              // If subcycling
          BOOST_TEST(iterationCount <= *iterValidIterations);
          BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint));
          BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint));
          BOOST_TEST(subcyclingStep < 2);
          subcyclingStep++;
        }
      }
    }
    cplScheme.finalize(); // Ends the coupling scheme
    BOOST_TEST(testing::equals(computedTime, 0.3));
    BOOST_TEST(computedTimesteps == 3);
    BOOST_TEST(stepsizeData0 == 5.0);
  }

  else if (nameParticipant == nameParticipant1) {
    iterationCount++; // different handling due to subcycling
    mesh->data(1)->setSampleAtTime(0, time::Sample{mesh->data(1)->getDimensions(), mesh->data(1)->values()});
    cplScheme.initialize();
    BOOST_TEST(not cplScheme.isTimeWindowComplete());
    BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint));
    BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint));
    BOOST_TEST(cplScheme.hasDataBeenReceived());

    // Tells coupling scheme, that a checkpoint has been created.
    // All required actions have to be performed before calling advance().
    cplScheme.markActionFulfilled(CouplingScheme::Action::WriteCheckpoint);
    BOOST_TEST(cplScheme.isActionFulfilled(CouplingScheme::Action::WriteCheckpoint));

    double maxTimeStepSize       = cplScheme.getNextTimeStepMaxSize();
    double preferredTimeStepSize = maxTimeStepSize / 2.5;
    double computedTimeStepSize  = preferredTimeStepSize;
    int    subcyclingStep        = 0;

    // Clear data for iteration.
    mesh->data(1)->timeStepsStorage().trim();

    // Main coupling loop
    while (cplScheme.isCouplingOngoing()) {
      cplScheme.addComputedTime(computedTimeStepSize);
      mesh->data(1)->setSampleAtTime(cplScheme.getTime(), time::Sample{mesh->data(1)->getDimensions(), mesh->data(1)->values()});
      cplScheme.firstSynchronization({});
      cplScheme.firstExchange();
      cplScheme.secondSynchronization();
      cplScheme.secondExchange();
      computedTimeStepSize =
          cplScheme.getNextTimeStepMaxSize() < preferredTimeStepSize
              ? cplScheme.getNextTimeStepMaxSize()
              : preferredTimeStepSize;
      // A coupling time step is complete, when the coupling iterations are
      // globally converged and if subcycling steps have filled one global
      // time step.
      if (cplScheme.isTimeWindowComplete()) {
        mesh->data(1)->timeStepsStorage().trim();
        // Advance participant time and time step
        computedTime += maxTimeStepSize;
        computedTimesteps++;
        BOOST_TEST(testing::equals(computedTime, cplScheme.getTime()));
        BOOST_TEST(computedTimesteps == cplScheme.getTimeWindows() - 1);
        // The iteration number is enforced by the controlled decrease of the
        // change of data written
        BOOST_TEST(iterationCount == *iterValidIterations);
        if (cplScheme.isCouplingOngoing()) {
          BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint));
          cplScheme.markActionFulfilled(CouplingScheme::Action::WriteCheckpoint);
          BOOST_TEST(cplScheme.isActionFulfilled(CouplingScheme::Action::WriteCheckpoint));
        } else {
          BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint));
          BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint));
        }
        iterationCount = 1;
        iterValidIterations++;
        if (iterValidIterations == validIterations.end()) {
          BOOST_REQUIRE(not cplScheme.isCouplingOngoing());
        }
        // Reset data values, to simulate same convergence behavior of
        // interface values in next time step.
        stepsizeData1 = initialStepsizeData1;
        BOOST_TEST(subcyclingStep == 2);
        subcyclingStep = 0;
      } else { // coupling timestep is not yet complete
        BOOST_TEST(cplScheme.isCouplingOngoing());
        // If end of time window is reached
        if (cplScheme.hasDataBeenReceived()) {
          BOOST_TEST(iterationCount <= *iterValidIterations);
          BOOST_TEST(cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint));
          BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint));
          cplScheme.markActionFulfilled(CouplingScheme::Action::ReadCheckpoint);
          BOOST_TEST(cplScheme.isActionFulfilled(CouplingScheme::Action::ReadCheckpoint));
          // The written data value is decreased in a regular manner, in order
          // to achieve a predictable convergence.
          stepsizeData1.array() -= 1.0;
          subcyclingStep = 0; // Subcycling steps
          iterationCount++;   // Implicit coupling iterations
        } else {              // If subcycling
          BOOST_TEST(iterationCount <= *iterValidIterations);
          BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint));
          BOOST_TEST(not cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint));
          BOOST_TEST(subcyclingStep < 3);
          subcyclingStep++;
        }
      }
    }
    cplScheme.finalize(); // Ends the coupling scheme
    BOOST_TEST(testing::equals(computedTime, 0.3));
    BOOST_TEST(computedTimesteps == 3);
  }
}

struct SerialImplicitCouplingSchemeFixture : m2n::WhiteboxAccessor {
  std::string _pathToTests;

  SerialImplicitCouplingSchemeFixture()
  {
    _pathToTests = testing::getPathToSources() + "/cplscheme/tests/";
  }
};

BOOST_FIXTURE_TEST_SUITE(SerialImplicitCouplingSchemeTests, SerialImplicitCouplingSchemeFixture)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(testParseConfigurationWithRelaxation)
{
  PRECICE_TEST();
  using namespace mesh;

  std::string path(_pathToTests + "serial-implicit-cplscheme-relax-const-config.xml");

  xml::XMLTag                          root = xml::getRootTag();
  PtrDataConfiguration                 dataConfig(new DataConfiguration(root));
  PtrMeshConfiguration                 meshConfig(new MeshConfiguration(root, dataConfig));
  m2n::M2NConfiguration::SharedPointer m2nConfig(
      new m2n::M2NConfiguration(root));
  precice::config::PtrParticipantConfiguration participantConfig(new precice::config::ParticipantConfiguration(root, meshConfig));
  CouplingSchemeConfiguration                  cplSchemeConfig(root, meshConfig, m2nConfig, participantConfig);

  xml::configure(root, xml::ConfigurationContext{}, path);
  BOOST_CHECK(cplSchemeConfig.getData("Data0", "Mesh") != cplSchemeConfig.getData("Data1", "Mesh"));
  BOOST_CHECK(cplSchemeConfig.findDataByID(cplSchemeConfig.getData("Data0", "Mesh")->getID()) != cplSchemeConfig.findDataByID(cplSchemeConfig.getData("Data1", "Mesh")->getID()));
  BOOST_CHECK(cplSchemeConfig.getData("Data0", "Mesh") == cplSchemeConfig.findDataByID(cplSchemeConfig.getData("Data0", "Mesh")->getID()));
  BOOST_CHECK(cplSchemeConfig.getData("Data1", "Mesh") == cplSchemeConfig.findDataByID(cplSchemeConfig.getData("Data1", "Mesh")->getID()));
  BOOST_CHECK(cplSchemeConfig.findDataByID(2) == nullptr);                   // nullptr, there are only two pieces of data.
  BOOST_CHECK(cplSchemeConfig._accelerationConfig->getAcceleration().get()); // no nullptr
}

/// Test that runs on 2 processors.
PRECICE_TEST_SETUP("Participant0"_on(1_rank), "Participant1"_on(1_rank), Require::Events)
BOOST_AUTO_TEST_CASE(testAbsConvergenceMeasureSynchronized)
{
  PRECICE_TEST();
  testing::ConnectionOptions options;
  options.useOnlyPrimaryCom = true;
  auto m2n                  = context.connectPrimaryRanks("Participant0", "Participant1", options);

  using namespace mesh;

  xml::XMLTag root = xml::getRootTag();
  // Create a data configuration, to simplify configuration of data
  PtrDataConfiguration dataConfig(new DataConfiguration(root));
  dataConfig->addData("data0", mesh::Data::typeName::SCALAR);
  dataConfig->addData("data1", mesh::Data::typeName::VECTOR);

  MeshConfiguration meshConfig(root, dataConfig);
  mesh::PtrMesh     mesh(new Mesh("Mesh", 3, testing::nextMeshID()));
  mesh->createData("data0", 1, 0_dataID);
  mesh->createData("data1", 3, 1_dataID);
  mesh->createVertex(Eigen::Vector3d::Zero());
  mesh->allocateDataValues();
  meshConfig.insertMeshToMeshDimensionsMap(mesh->getName(), mesh->getDimensions());
  meshConfig.addMesh(mesh);

  // Create all parameters necessary to create an ImplicitCouplingScheme object
  const double maxTime        = 1.0;
  const int    maxTimeWindows = 3;
  const double timeWindowSize = 0.1;
  std::string  nameParticipant0("Participant0");
  std::string  nameParticipant1("Participant1");
  int          sendDataIndex        = -1;
  int          receiveDataIndex     = -1;
  int          convergenceDataIndex = -1;
  if (context.isNamed(nameParticipant0)) {
    sendDataIndex        = 0;
    receiveDataIndex     = 1;
    convergenceDataIndex = receiveDataIndex;
  } else {
    sendDataIndex        = 1;
    receiveDataIndex     = 0;
    convergenceDataIndex = sendDataIndex;
  }

  // Create the coupling scheme object
  const int                       minIterations = 1;
  const int                       maxIterations = 100;
  cplscheme::SerialCouplingScheme cplScheme(maxTime, maxTimeWindows, timeWindowSize, nameParticipant0, nameParticipant1, context.name, m2n, constants::FIXED_TIME_WINDOW_SIZE, BaseCouplingScheme::Implicit, minIterations, maxIterations);
  cplScheme.addDataToSend(mesh->data(sendDataIndex), mesh, false, false, true);
  cplScheme.addDataToReceive(mesh->data(receiveDataIndex), mesh, false, false, true);
  cplScheme.determineInitialDataExchange();

  double                                 convergenceLimit1 = sqrt(3.0); // when diff_vector = (1.0, 1.0, 1.0)
  cplscheme::impl::PtrConvergenceMeasure absoluteConvMeasure1(
      new cplscheme::impl::AbsoluteConvergenceMeasure(convergenceLimit1));
  cplScheme.addConvergenceMeasure(convergenceDataIndex, false, false, absoluteConvMeasure1);

  // Expected iterations per implicit timesptep
  std::vector<int> validIterations = {5, 5, 5};
  runCoupling(cplScheme, context.name, meshConfig, validIterations);
}

PRECICE_TEST_SETUP("Participant0"_on(1_rank), "Participant1"_on(1_rank), Require::Events)
BOOST_AUTO_TEST_CASE(testConfiguredAbsConvergenceMeasureSynchronized)
{
  PRECICE_TEST();

  using namespace mesh;

  std::string configurationPath(
      _pathToTests + "serial-implicit-cplscheme-absolute-config.xml");

  xml::XMLTag                                  root = xml::getRootTag();
  PtrDataConfiguration                         dataConfig(new DataConfiguration(root));
  PtrMeshConfiguration                         meshConfig(new MeshConfiguration(root, dataConfig));
  m2n::M2NConfiguration::SharedPointer         m2nConfig(new m2n::M2NConfiguration(root));
  precice::config::PtrParticipantConfiguration participantConfig(new precice::config::ParticipantConfiguration(root, meshConfig));
  CouplingSchemeConfiguration                  cplSchemeConfig(root, meshConfig, m2nConfig, participantConfig);

  xml::configure(root, xml::ConfigurationContext{}, configurationPath);
  m2n::PtrM2N m2n        = m2nConfig->getM2N("Participant0", "Participant1");
  useOnlyPrimaryCom(m2n) = true;

  // some dummy mesh
  meshConfig->meshes().at(0)->createVertex(Eigen::Vector3d(1.0, 1.0, 1.0));
  meshConfig->meshes().at(0)->createVertex(Eigen::Vector3d(2.0, 1.0, -1.0));
  meshConfig->meshes().at(0)->createVertex(Eigen::Vector3d(3.0, 1.0, 1.0));
  meshConfig->meshes().at(0)->createVertex(Eigen::Vector3d(4.0, 1.0, -1.0));
  meshConfig->meshes().at(0)->allocateDataValues();

  std::vector<int> validIterations = {5, 5, 5};

  if (context.isNamed("Participant0")) {
    m2n->requestPrimaryRankConnection("Participant1", "Participant0", "");
  } else {
    m2n->acceptPrimaryRankConnection("Participant1", "Participant0", "");
  }

  runCoupling(*cplSchemeConfig.getCouplingScheme(context.name),
              context.name, *meshConfig, validIterations);
}

PRECICE_TEST_SETUP("Participant0"_on(1_rank), "Participant1"_on(1_rank), Require::Events)
BOOST_AUTO_TEST_CASE(testMinIterConvergenceMeasureSynchronized)
{
  PRECICE_TEST();
  testing::ConnectionOptions options;
  options.useOnlyPrimaryCom = true;
  auto m2n                  = context.connectPrimaryRanks("Participant0", "Participant1", options);

  xml::XMLTag root = xml::getRootTag();
  // Create a data configuration, to simplify configuration of data
  mesh::PtrDataConfiguration dataConfig(new mesh::DataConfiguration(root));
  dataConfig->addData("data0", mesh::Data::typeName::SCALAR);
  dataConfig->addData("data1", mesh::Data::typeName::VECTOR);

  mesh::MeshConfiguration meshConfig(root, dataConfig);
  mesh::PtrMesh           mesh(new mesh::Mesh("Mesh", 3, testing::nextMeshID()));
  mesh->createData("data0", 1, 0_dataID);
  mesh->createData("data1", 3, 1_dataID);
  mesh->createVertex(Eigen::Vector3d::Zero());
  mesh->allocateDataValues();
  meshConfig.insertMeshToMeshDimensionsMap(mesh->getName(), mesh->getDimensions());
  meshConfig.addMesh(mesh);

  // Create all parameters necessary to create an ImplicitCouplingScheme object
  const double maxTime        = 1.0;
  const int    maxTimeWindows = 3;
  const double timeWindowSize = 0.1;
  std::string  nameParticipant0("Participant0");
  std::string  nameParticipant1("Participant1");
  int          sendDataIndex    = -1;
  int          receiveDataIndex = -1;
  if (context.isNamed(nameParticipant0)) {
    sendDataIndex    = 0;
    receiveDataIndex = 1;
  } else {
    sendDataIndex    = 1;
    receiveDataIndex = 0;
  }

  // Create the coupling scheme object
  const int                       minIterations = 1;
  const int                       maxIterations = 3;
  cplscheme::SerialCouplingScheme cplScheme(maxTime, maxTimeWindows, timeWindowSize, nameParticipant0, nameParticipant1, context.name, m2n, constants::FIXED_TIME_WINDOW_SIZE, BaseCouplingScheme::Implicit, minIterations, maxIterations);
  cplScheme.addDataToSend(mesh->data(sendDataIndex), mesh, false, false, true);
  cplScheme.addDataToReceive(mesh->data(receiveDataIndex), mesh, false, false, true);
  cplScheme.determineInitialDataExchange();

  // Expected iterations per implicit timesptep
  std::vector<int> validIterations = {3, 3, 3};
  runCoupling(cplScheme, context.name, meshConfig, validIterations);
}

PRECICE_TEST_SETUP("Participant0"_on(1_rank), "Participant1"_on(1_rank), Require::Events)
BOOST_AUTO_TEST_CASE(testMinIterConvergenceMeasureSynchronizedWithSubcycling)
{
  PRECICE_TEST();
  testing::ConnectionOptions options;
  options.useOnlyPrimaryCom = true;
  auto m2n                  = context.connectPrimaryRanks("Participant0", "Participant1", options);

  xml::XMLTag root = xml::getRootTag();
  // Create a data configuration, to simplify configuration of data
  mesh::PtrDataConfiguration dataConfig(new mesh::DataConfiguration(root));
  dataConfig->addData("data0", mesh::Data::typeName::SCALAR);
  dataConfig->addData("data1", mesh::Data::typeName::VECTOR);

  mesh::MeshConfiguration meshConfig(root, dataConfig);
  mesh::PtrMesh           mesh(new mesh::Mesh("Mesh", 3, testing::nextMeshID()));
  mesh->createData("data0", 1, 0_dataID);
  mesh->createData("data1", 3, 1_dataID);
  mesh->createVertex(Eigen::Vector3d::Zero());
  mesh->allocateDataValues();
  meshConfig.insertMeshToMeshDimensionsMap(mesh->getName(), mesh->getDimensions());
  meshConfig.addMesh(mesh);

  // Create all parameters necessary to create an ImplicitCouplingScheme object
  double           maxTime        = 1.0;
  int              maxTimeWindows = 3;
  double           timeWindowSize = 0.1;
  std::string      nameParticipant0("Participant0");
  std::string      nameParticipant1("Participant1");
  int              sendDataIndex    = -1;
  int              receiveDataIndex = -1;
  std::vector<int> validIterations;
  if (context.isNamed(nameParticipant0)) {
    sendDataIndex    = 0;
    receiveDataIndex = 1;
    validIterations  = {3, 3, 3};
  } else {
    sendDataIndex    = 1;
    receiveDataIndex = 0;
    validIterations  = {3, 3, 3};
  }

  // Create the coupling scheme object
  const int                       minIterations = 1;
  const int                       maxIterations = 3;
  cplscheme::SerialCouplingScheme cplScheme(maxTime, maxTimeWindows, timeWindowSize, nameParticipant0, nameParticipant1, context.name, m2n, constants::FIXED_TIME_WINDOW_SIZE, BaseCouplingScheme::Implicit, minIterations, maxIterations);
  cplScheme.addDataToSend(mesh->data(sendDataIndex), mesh, false, false, true);
  cplScheme.addDataToReceive(mesh->data(receiveDataIndex), mesh, false, false, true);
  cplScheme.determineInitialDataExchange();

  runCouplingWithSubcycling(
      cplScheme, context.name, meshConfig, validIterations);
}

PRECICE_TEST_SETUP("Participant0"_on(1_rank), "Participant1"_on(1_rank), Require::Events)
BOOST_AUTO_TEST_CASE(testInitializeData)
{
  PRECICE_TEST();
  testing::ConnectionOptions options;
  options.useOnlyPrimaryCom = true;
  auto m2n                  = context.connectPrimaryRanks("Participant0", "Participant1", options);

  xml::XMLTag root = xml::getRootTag();

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

  // Create all parameters necessary to create an ImplicitCouplingScheme object
  const double maxTime        = 1.0;
  const int    maxTimeWindows = 3;
  const double timeWindowSize = 0.1;
  const double timeStepSize   = timeWindowSize; // solver is not subcycling
  std::string  nameParticipant0("Participant0");
  std::string  nameParticipant1("Participant1");
  int          sendDataIndex              = -1;
  int          receiveDataIndex           = -1;
  bool         dataRequiresInitialization = false;
  if (context.isNamed(nameParticipant0)) {
    sendDataIndex    = dataID0;
    receiveDataIndex = dataID1;
  } else {
    sendDataIndex              = dataID1;
    receiveDataIndex           = dataID0;
    dataRequiresInitialization = true;
  }

  // Create the coupling scheme object
  const int                       minIterations = 1;
  const int                       maxIterations = 3;
  cplscheme::SerialCouplingScheme cplScheme(maxTime, maxTimeWindows, timeWindowSize, nameParticipant0, nameParticipant1, context.name, m2n, constants::FIXED_TIME_WINDOW_SIZE, BaseCouplingScheme::Implicit, minIterations, maxIterations);
  using Fixture = testing::SerialCouplingSchemeFixture;

  cplScheme.addDataToSend(mesh->data(sendDataIndex), mesh, dataRequiresInitialization, false, true);
  CouplingData *sendCouplingData = Fixture::getSendData(cplScheme, sendDataIndex);
  cplScheme.addDataToReceive(mesh->data(receiveDataIndex), mesh, not dataRequiresInitialization, false, true);
  CouplingData *receiveCouplingData = Fixture::getReceiveData(cplScheme, receiveDataIndex);
  cplScheme.determineInitialDataExchange();

  if (context.isNamed(nameParticipant0)) {
    // ensure that read data is uninitialized
    BOOST_TEST(receiveCouplingData->getSize() == 3);
    BOOST_TEST(testing::equals(receiveCouplingData->values(), Eigen::Vector3d(0.0, 0.0, 0.0)));
    // ensure that write data is uninitialized
    BOOST_TEST(sendCouplingData->getSize() == 1);
    BOOST_TEST(testing::equals(sendCouplingData->values()(0), 0.0));

    BOOST_TEST(Fixture::isImplicitCouplingScheme(cplScheme));
    sendCouplingData->setSampleAtTime(0, time::Sample{sendCouplingData->getDimensions(), sendCouplingData->values()});
    cplScheme.initialize();
    BOOST_TEST(cplScheme.hasDataBeenReceived());
    // ensure that initial data was read
    BOOST_TEST(receiveCouplingData->getSize() == 3);
    BOOST_TEST(testing::equals(receiveCouplingData->values(), Eigen::Vector3d(1.0, 2.0, 3.0)));
    BOOST_TEST(receiveCouplingData->getPreviousIterationSize() == 3);
    BOOST_TEST(testing::equals(receiveCouplingData->previousIteration(), Eigen::Vector3d(0.0, 0.0, 0.0)));
    // ensure that write data is still uninitialized
    BOOST_TEST(sendCouplingData->getSize() == 1);
    BOOST_TEST(testing::equals(sendCouplingData->values()(0), 0.0));
    BOOST_TEST(sendCouplingData->getPreviousIterationSize() == 1);
    BOOST_TEST(testing::equals(sendCouplingData->previousIteration()(0), 0.0));
    BOOST_TEST(sendCouplingData->getPreviousIterationSize() == 1);
    // set write data
    while (cplScheme.isCouplingOngoing()) {
      if (cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint)) {
        cplScheme.markActionFulfilled(CouplingScheme::Action::WriteCheckpoint);
      }
      BOOST_TEST(cplScheme.getNextTimeStepMaxSize() == timeStepSize);
      sendCouplingData->setSampleAtTime(cplScheme.getTime() + timeStepSize, time::Sample{sendCouplingData->getDimensions(), Eigen::VectorXd::Constant(sendCouplingData->getSize(), 4.0)});
      cplScheme.addComputedTime(timeStepSize);
      sendCouplingData->setSampleAtTime(cplScheme.getTime(), time::Sample{sendCouplingData->getDimensions(), Eigen::VectorXd::Constant(sendCouplingData->getSize(), 4.0)});
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
    sendCouplingData->setSampleAtTime(0, time::Sample{sendCouplingData->getDimensions(), v});
    cplScheme.markActionFulfilled(CouplingScheme::Action::InitializeData);
    BOOST_TEST(receiveCouplingData->getSize() == 1);
    BOOST_TEST(testing::equals(receiveCouplingData->values()(0), 0.0));
    BOOST_TEST(sendCouplingData->getSize() == 3);
    BOOST_TEST(testing::equals(sendCouplingData->values(), Eigen::Vector3d(1.0, 2.0, 3.0)));
    cplScheme.initialize();
    BOOST_TEST(cplScheme.hasDataBeenReceived());
    BOOST_TEST(receiveCouplingData->getSize() == 1);
    BOOST_TEST(testing::equals(receiveCouplingData->values()(0), 4.0));
    BOOST_TEST(receiveCouplingData->getPreviousIterationSize() == 1);
    BOOST_TEST(testing::equals(receiveCouplingData->previousIteration()(0), 0.0));
    BOOST_TEST(sendCouplingData->getSize() == 3);
    BOOST_TEST(testing::equals(sendCouplingData->values(), Eigen::Vector3d(1.0, 2.0, 3.0)));
    BOOST_TEST(sendCouplingData->getPreviousIterationSize() == 3);
    BOOST_TEST(testing::equals(sendCouplingData->previousIteration(), Eigen::Vector3d(1.0, 2.0, 3.0)));
    while (cplScheme.isCouplingOngoing()) {
      if (cplScheme.isActionRequired(CouplingScheme::Action::WriteCheckpoint)) {
        cplScheme.markActionFulfilled(CouplingScheme::Action::WriteCheckpoint);
      }
      BOOST_TEST(cplScheme.getNextTimeStepMaxSize() == timeStepSize);
      sendCouplingData->setSampleAtTime(cplScheme.getTime() + timeStepSize, time::Sample{sendCouplingData->getDimensions(), v});
      cplScheme.addComputedTime(timeStepSize);
      cplScheme.firstSynchronization({});
      cplScheme.firstExchange();
      cplScheme.secondSynchronization();
      cplScheme.secondExchange();
      if (cplScheme.isCouplingOngoing()) {
        BOOST_TEST(cplScheme.hasDataBeenReceived());
      }
      if (cplScheme.isActionRequired(CouplingScheme::Action::ReadCheckpoint)) {
        cplScheme.markActionFulfilled(CouplingScheme::Action::ReadCheckpoint);
      }
    }
  }
  cplScheme.finalize();
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // not PRECICE_NO_MPI

#ifndef PRECICE_NO_MPI

#include <Eigen/Core>
#include <algorithm>
#include <deque>
#include <fstream>
#include <istream>
#include <iterator>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "action/RecorderAction.hpp"
#include "logging/LogMacros.hpp"
#include "math/constants.hpp"
#include "math/geometry.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Utils.hpp"
#include "mesh/Vertex.hpp"
#include "precice/SolverInterface.hpp"
#include "precice/impl/MeshContext.hpp"
#include "precice/impl/Participant.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"
#include "precice/types.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using precice::testing::TestContext;

struct SerialTestFixture : testing::WhiteboxAccessor {

  std::string _pathToTests;

  void reset()
  {
    mesh::Data::resetDataCount();
  }

  SerialTestFixture()
  {
    _pathToTests = testing::getPathToSources() + "/precice/tests/";
    reset();
  }
};

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_FIXTURE_TEST_SUITE(Serial, SerialTestFixture)
BOOST_AUTO_TEST_SUITE(InitializeData)

/**
 * @brief helper function for a simple test with data initialization
 */
void testDataInitialization(precice::testing::TestContext context, std::string config)
{
  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, config, 0, 1);
  if (context.isNamed("SolverOne")) {
    int      meshOneID = cplInterface.getMeshID("MeshOne");
    Vector3d pos       = Vector3d::Zero();
    cplInterface.setMeshVertex(meshOneID, pos.data());
    double maxDt      = cplInterface.initialize();
    int    dataID     = cplInterface.getDataID("Data", meshOneID);
    double valueDataB = 0.0;
    cplInterface.initializeData();
    cplInterface.readScalarData(dataID, 0, valueDataB);
    BOOST_TEST(2.0 == valueDataB);
    cplInterface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    int      meshTwoID = cplInterface.getMeshID("MeshTwo");
    Vector3d pos       = Vector3d::Zero();
    cplInterface.setMeshVertex(meshTwoID, pos.data());
    double maxDt  = cplInterface.initialize();
    int    dataID = cplInterface.getDataID("Data", meshTwoID);
    cplInterface.writeScalarData(dataID, 0, 2.0);
    //tell preCICE that data has been written and call initializeData
    cplInterface.markActionFulfilled(precice::constants::actionWriteInitialData());
    cplInterface.initializeData();
    cplInterface.finalize();
  }
}

/**
 * @brief The second solver initializes the data of the first. Use write mapping for data.
 */
BOOST_AUTO_TEST_CASE(testDataInitializationWriteMapping)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  testDataInitialization(context, _pathToTests + "oneway-data-init-write-mapping.xml");
}

/**
 * @brief The second solver initializes the data of the first. Use read mapping for data.
 */
BOOST_AUTO_TEST_CASE(testDataInitializationReadMapping)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  testDataInitialization(context, _pathToTests + "oneway-data-init-read-mapping.xml");
}

/**
 * @brief The second solver initializes the data of the first.
 *
 * A mapping is employed for the second solver, i.e., at the end of
 * initializeData(), the mapping needs to be invoked.
 */
BOOST_AUTO_TEST_CASE(testExplicitWithDataInitialization)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, _pathToTests + "explicit-data-init.xml", 0, 1);
  if (context.isNamed("SolverOne")) {
    int meshOneID = cplInterface.getMeshID("MeshOne");
    cplInterface.setMeshVertex(meshOneID, Vector3d(1.0, 2.0, 3.0).data());
    double maxDt      = cplInterface.initialize();
    int    dataAID    = cplInterface.getDataID("DataOne", meshOneID);
    int    dataBID    = cplInterface.getDataID("DataTwo", meshOneID);
    double valueDataB = 0.0;
    cplInterface.initializeData();
    cplInterface.readScalarData(dataBID, 0, valueDataB);
    BOOST_TEST(2.0 == valueDataB);
    while (cplInterface.isCouplingOngoing()) {
      Vector3d valueDataA(1.0, 1.0, 1.0);
      cplInterface.writeVectorData(dataAID, 0, valueDataA.data());
      maxDt = cplInterface.advance(maxDt);
      cplInterface.readScalarData(dataBID, 0, valueDataB);
      BOOST_TEST(2.5 == valueDataB);
    }
    cplInterface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    int      meshTwoID = cplInterface.getMeshID("MeshTwo");
    Vector3d pos       = Vector3d::Zero();
    cplInterface.setMeshVertex(meshTwoID, pos.data());
    double maxDt   = cplInterface.initialize();
    int    dataAID = cplInterface.getDataID("DataOne", meshTwoID);
    int    dataBID = cplInterface.getDataID("DataTwo", meshTwoID);
    cplInterface.writeScalarData(dataBID, 0, 2.0);
    //tell preCICE that data has been written and call initializeData
    cplInterface.markActionFulfilled(precice::constants::actionWriteInitialData());
    cplInterface.initializeData();
    Vector3d valueDataA;
    cplInterface.readVectorData(dataAID, 0, valueDataA.data());
    Vector3d expected(1.0, 1.0, 1.0);
    BOOST_TEST(valueDataA == expected);
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(dataBID, 0, 2.5);
      maxDt = cplInterface.advance(maxDt);
      cplInterface.readVectorData(dataAID, 0, valueDataA.data());
      BOOST_TEST(valueDataA == expected);
    }
    cplInterface.finalize();
  }
}

/// Test simple coupled simulation with iterations, data initialization and without acceleration
BOOST_AUTO_TEST_CASE(testImplicitWithDataInitialization)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  using namespace precice::constants;

  SolverInterface couplingInterface(context.name, _pathToTests + "implicit-data-init.xml", 0, 1);

  int         dimensions = couplingInterface.getDimensions();
  std::string meshName;
  std::string writeDataName;
  std::string readDataName;
  double      writeValue, expectedReadValue;

  if (context.isNamed("SolverOne")) {
    meshName          = "MeshOne";
    writeDataName     = "Forces";
    readDataName      = "Velocities";
    writeValue        = 1;
    expectedReadValue = 2;
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    meshName          = "MeshTwo";
    writeDataName     = "Velocities";
    readDataName      = "Forces";
    writeValue        = 2;
    expectedReadValue = 1;
  }
  int                 meshID      = couplingInterface.getMeshID(meshName);
  int                 writeDataID = couplingInterface.getDataID(writeDataName, meshID);
  int                 readDataID  = couplingInterface.getDataID(readDataName, meshID);
  std::vector<double> vertex(dimensions, 0);
  int                 vertexID = couplingInterface.setMeshVertex(meshID, vertex.data());

  double dt = 0;
  dt        = couplingInterface.initialize();
  std::vector<double> writeData(dimensions, writeValue);
  std::vector<double> readData(dimensions, -1);
  const std::string & cowid = actionWriteInitialData();

  if (couplingInterface.isActionRequired(cowid)) {
    BOOST_TEST(context.isNamed("SolverTwo"));
    couplingInterface.writeVectorData(writeDataID, vertexID, writeData.data());
    couplingInterface.markActionFulfilled(cowid);
  }

  couplingInterface.initializeData();

  while (couplingInterface.isCouplingOngoing()) {
    if (couplingInterface.isActionRequired(actionWriteIterationCheckpoint())) {
      couplingInterface.markActionFulfilled(actionWriteIterationCheckpoint());
    }
    couplingInterface.readVectorData(readDataID, vertexID, readData.data());
    BOOST_TEST(expectedReadValue == readData.at(0));
    BOOST_TEST(expectedReadValue == readData.at(1));
    couplingInterface.writeVectorData(writeDataID, vertexID, writeData.data());
    dt = couplingInterface.advance(dt);
    if (couplingInterface.isActionRequired(actionReadIterationCheckpoint())) {
      couplingInterface.markActionFulfilled(actionReadIterationCheckpoint());
    }
  }
  couplingInterface.finalize();
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
#endif // PRECICE_NO_MPI

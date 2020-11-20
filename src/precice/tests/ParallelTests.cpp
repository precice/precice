#ifndef PRECICE_NO_MPI
#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include <mpi.h>
#include <string>
#include <vector>
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "math/constants.hpp"
#include "mesh/Mesh.hpp"
#include "precice/SolverInterface.hpp"
#include "precice/config/Configuration.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/Parallel.hpp"

using namespace precice;
using testing::TestContext;

struct ParallelTestFixture : testing::WhiteboxAccessor {

  std::string _pathToTests;

  ParallelTestFixture()
  {
    _pathToTests = testing::getPathToSources() + "/precice/tests/";
  }
};

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_FIXTURE_TEST_SUITE(Parallel, ParallelTestFixture)

BOOST_AUTO_TEST_CASE(TestMasterSlaveSetup)
{
  PRECICE_TEST("SolverOne"_on(4_ranks));
  std::string     configFilename = _pathToTests + "config1.xml";
  SolverInterface interface(context.name, configFilename, context.rank, context.size);
  BOOST_TEST(interface.getDimensions() == 3);

  if (context.isMaster()) {
    BOOST_TEST(utils::MasterSlave::isMaster() == true);
    BOOST_TEST(utils::MasterSlave::isSlave() == false);
  } else {
    BOOST_TEST(utils::MasterSlave::isMaster() == false);
    BOOST_TEST(utils::MasterSlave::isSlave() == true);
  }

  BOOST_TEST(utils::MasterSlave::getRank() == context.rank);
  BOOST_TEST(utils::MasterSlave::getSize() == context.size);
  BOOST_TEST(utils::MasterSlave::_communication.use_count() > 0);
  BOOST_TEST(utils::MasterSlave::_communication->isConnected());
}

BOOST_AUTO_TEST_CASE(TestFinalize)
{
  PRECICE_TEST("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks));
  std::string configFilename = _pathToTests + "config1.xml";
  if (context.isNamed("SolverOne")) {
    SolverInterface interface(context.name, configFilename, context.rank, context.size);
    int             meshID = interface.getMeshID("MeshOne");
    double          xCoord = 0.0 + context.rank;
    interface.setMeshVertex(meshID, Eigen::Vector3d(xCoord, 0.0, 0.0).data());
    interface.initialize();
    BOOST_TEST(impl(interface).mesh("MeshOne").vertices().size() == 1);
    BOOST_TEST(impl(interface).mesh("MeshTwo").vertices().size() == 1);
    interface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    SolverInterface interface(context.name, configFilename, context.rank, context.size);
    int             meshID = interface.getMeshID("MeshTwo");
    double          xCoord = 0.0 + context.rank;
    interface.setMeshVertex(meshID, Eigen::Vector3d(xCoord, 0.0, 0.0).data());
    interface.initialize();
    BOOST_TEST(impl(interface).mesh("MeshTwo").vertices().size() == 1);
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE(Lifecycle)

// Test representing the full explicit lifecycle of a SolverInterface
BOOST_AUTO_TEST_CASE(Full)
{
  PRECICE_TEST("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks));
  std::string config = _pathToTests + "lifecycle.xml";

  SolverInterface interface(context.name, config, context.rank, context.size);

  constexpr double y{0};
  constexpr double z{0};
  constexpr double x1{1};
  constexpr double dx{1};

  if (context.isNamed("SolverOne")) {
    auto   meshid   = interface.getMeshID("MeshOne");
    double coords[] = {x1 + dx * context.rank, y, z};
    auto   vertexid = interface.setMeshVertex(meshid, coords);

    auto   dataid = interface.getDataID("DataOne", meshid);
    double data[] = {3.4, 4.5, 5.6};
    interface.writeVectorData(dataid, vertexid, data);
  } else {
    auto   meshid   = interface.getMeshID("MeshTwo");
    double coords[] = {x1 + dx * context.rank, y, z};
    auto   vertexid = interface.setMeshVertex(meshid, coords);

    auto dataid = interface.getDataID("DataTwo", meshid);
    interface.writeScalarData(dataid, vertexid, 7.8);
  }
  interface.initialize();
  BOOST_TEST(interface.isCouplingOngoing());
  interface.finalize();
}

// Test representing the full lifecycle of a SolverInterface
// Finalize is not called explicitly here.
// The destructor has to cleanup.
BOOST_AUTO_TEST_CASE(ImplicitFinalize)
{
  PRECICE_TEST("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks));
  std::string config = _pathToTests + "lifecycle.xml";

  SolverInterface interface(context.name, config, context.rank, context.size);

  constexpr double y{0};
  constexpr double z{0};
  constexpr double x1{1};
  constexpr double dx{1};

  if (context.isNamed("SolverOne")) {
    auto   meshid   = interface.getMeshID("MeshOne");
    double coords[] = {x1 + dx * context.rank, y, z};
    auto   vertexid = interface.setMeshVertex(meshid, coords);

    auto   dataid = interface.getDataID("DataOne", meshid);
    double data[] = {3.4, 4.5, 5.6};
    interface.writeVectorData(dataid, vertexid, data);
  } else {
    auto   meshid   = interface.getMeshID("MeshTwo");
    double coords[] = {x1 + dx * context.rank, y, z};
    auto   vertexid = interface.setMeshVertex(meshid, coords);

    auto dataid = interface.getDataID("DataTwo", meshid);
    interface.writeScalarData(dataid, vertexid, 7.8);
  }
  interface.initialize();
  BOOST_TEST(interface.isCouplingOngoing());
}

// Test representing the minimal lifecylce, which consists out of construction only.
// The destructor has to cleanup correctly.
BOOST_AUTO_TEST_CASE(ConstructOnly)
{
  PRECICE_TEST("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks));
  std::string config = _pathToTests + "lifecycle.xml";

  SolverInterface interface(context.name, config, context.rank, context.size);
}

// Test representing the minimal lifecylce with explicit finalization.
// This shows how to manually finalize MPI etc without using the SolverInterface.
BOOST_AUTO_TEST_CASE(ConstructAndExplicitFinalize)
{
  PRECICE_TEST("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks));
  std::string config = _pathToTests + "lifecycle.xml";

  SolverInterface interface(context.name, config, context.rank, context.size);

  interface.finalize();
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_CASE(GlobalRBFPartitioning)
{
  PRECICE_TEST("SolverOne"_on(3_ranks), "SolverTwo"_on(1_rank));
  std::string configFilename = _pathToTests + "globalRBFPartitioning.xml";

  if (context.isNamed("SolverOne")) {
    SolverInterface interface(context.name, configFilename, context.rank, context.size);
    int             meshID = interface.getMeshID("MeshOne");
    int             dataID = interface.getDataID("Data2", meshID);

    int    vertexIDs[2];
    double xCoord       = context.rank * 0.4;
    double positions[4] = {xCoord, 0.0, xCoord + 0.2, 0.0};
    interface.setMeshVertices(meshID, 2, positions, vertexIDs);
    interface.initialize();
    double values[2];
    interface.advance(1.0);
    interface.readBlockScalarData(dataID, 2, vertexIDs, values);
    //    std::cout << context.rank <<": " << values << '\n';
    interface.finalize();
  } else {
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    SolverInterface interface(context.name, configFilename, context.rank, context.size);
    int             meshID = interface.getMeshID("MeshTwo");
    int             vertexIDs[6];
    double          positions[12] = {0.0, 0.0, 0.2, 0.0, 0.4, 0.0, 0.6, 0.0, 0.8, 0.0, 1.0, 0.0};
    interface.setMeshVertices(meshID, 6, positions, vertexIDs);
    interface.initialize();
    int    dataID    = interface.getDataID("Data2", meshID);
    double values[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    interface.writeBlockScalarData(dataID, 6, vertexIDs, values);
    interface.advance(1.0);
    interface.finalize();
  }
}

BOOST_AUTO_TEST_CASE(LocalRBFPartitioning)
{
  PRECICE_TEST("SolverOne"_on(3_ranks), "SolverTwo"_on(1_rank));
  std::string configFilename = _pathToTests + "localRBFPartitioning.xml";

  if (context.name == "SolverOne") {
    SolverInterface interface(context.name, configFilename, context.rank, context.size);
    int             meshID = interface.getMeshID("MeshOne");
    int             dataID = interface.getDataID("Data2", meshID);

    int    vertexIDs[2];
    double xCoord       = context.rank * 0.4;
    double positions[4] = {xCoord, 0.0, xCoord + 0.2, 0.0};
    interface.setMeshVertices(meshID, 2, positions, vertexIDs);
    interface.initialize();
    double values[2];
    interface.advance(1.0);
    interface.readBlockScalarData(dataID, 2, vertexIDs, values);
    interface.finalize();
  } else {
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    SolverInterface interface(context.name, configFilename, context.rank, context.size);
    int             meshID = interface.getMeshID("MeshTwo");
    int             vertexIDs[6];
    double          positions[12] = {0.0, 0.0, 0.2, 0.0, 0.4, 0.0, 0.6, 0.0, 0.8, 0.0, 1.0, 0.0};
    interface.setMeshVertices(meshID, 6, positions, vertexIDs);
    interface.initialize();
    int    dataID    = interface.getDataID("Data2", meshID);
    double values[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    interface.writeBlockScalarData(dataID, 6, vertexIDs, values);
    interface.advance(1.0);
    interface.finalize();
  }
}

/// This testcase is based on a bug reported by Thorsten for acoustic FASTEST-Ateles coupling
BOOST_AUTO_TEST_CASE(CouplingOnLine)
{
  PRECICE_TEST("Ateles"_on(3_ranks), "FASTEST"_on(1_rank));
  std::string configFilename = _pathToTests + "line-coupling.xml";

  if (context.isNamed("Ateles")) {
    SolverInterface interface(context.name, configFilename, context.rank, context.size);
    int             meshID = interface.getMeshID("Ateles_Mesh");

    int    vertexIDs[4];
    double offset        = context.rank * 0.4;
    double xCoord        = 0.0;
    double yCoord        = 1.0;
    double positions[12] = {xCoord, yCoord, 0.1 + offset,
                            xCoord, yCoord, 0.2 + offset,
                            xCoord, yCoord, 0.3 + offset,
                            xCoord, yCoord, 0.4 + offset};
    interface.setMeshVertices(meshID, 4, positions, vertexIDs);
    interface.initialize();
    interface.advance(1.0);
    interface.finalize();
  } else {
    SolverInterface interface(context.name, configFilename, context.rank, context.size);
    int             meshID = interface.getMeshID("FASTEST_Mesh");
    int             vertexIDs[10];
    double          xCoord        = -0.0001;
    double          yCoord        = 1.00001;
    double          positions[30] = {xCoord, yCoord, 0.12,
                            xCoord, yCoord, 0.24,
                            xCoord, yCoord, 0.36,
                            xCoord, yCoord, 0.48,
                            xCoord, yCoord, 0.60,
                            xCoord, yCoord, 0.72,
                            xCoord, yCoord, 0.84,
                            xCoord, yCoord, 0.96,
                            xCoord, yCoord, 1.08,
                            xCoord, yCoord, 1.2};
    interface.setMeshVertices(meshID, 10, positions, vertexIDs);
    interface.initialize();
    interface.advance(1.0);
    interface.finalize();
  }
}

/// tests for different QN settings if correct fixed point is reached
void runTestQN(std::string const &config, TestContext const &context)
{
  std::string meshName, writeDataName, readDataName;

  if (context.isNamed("SolverOne")) {
    meshName      = "MeshOne";
    writeDataName = "Data1";
    readDataName  = "Data2";
  } else {
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    meshName      = "MeshTwo";
    writeDataName = "Data2";
    readDataName  = "Data1";
  }

  SolverInterface interface(context.name, config, context.rank, context.size);
  int             meshID      = interface.getMeshID(meshName);
  int             writeDataID = interface.getDataID(writeDataName, meshID);
  int             readDataID  = interface.getDataID(readDataName, meshID);

  int vertexIDs[4];

  // meshes for rank 0 and rank 1, we use matching meshes for both participants
  double positions0[8] = {1.0, 0.0, 1.0, 0.5, 1.0, 1.0, 1.0, 1.5};
  double positions1[8] = {2.0, 0.0, 2.0, 0.5, 2.0, 1.0, 2.0, 1.5};

  if (context.isNamed("SolverOne")) {
    if (context.isMaster()) {
      interface.setMeshVertices(meshID, 4, positions0, vertexIDs);
    } else {
      interface.setMeshVertices(meshID, 4, positions1, vertexIDs);
    }
  } else {
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    if (context.isMaster()) {
      interface.setMeshVertices(meshID, 4, positions0, vertexIDs);
    } else {
      interface.setMeshVertices(meshID, 4, positions1, vertexIDs);
    }
  }

  interface.initialize();
  double inValues[4]  = {0.0, 0.0, 0.0, 0.0};
  double outValues[4] = {0.0, 0.0, 0.0, 0.0};

  int iterations = 0;

  while (interface.isCouplingOngoing()) {
    if (interface.isActionRequired(precice::constants::actionWriteIterationCheckpoint())) {
      interface.markActionFulfilled(precice::constants::actionWriteIterationCheckpoint());
    }

    if (interface.isReadDataAvailable())
      interface.readBlockScalarData(readDataID, 4, vertexIDs, inValues);

    /*
      Solves the following non-linear equations, which are extended to a fixed-point equation (simply +x)
      2 * x_1^2 - x_2 * x_3 - 8 = 0
      x_1^2 * x_2 + 2 * x_1 * x_2 * x_3 + x_2 * x_3^2 + x_2 = 0
      x_3^2 - 4 = 0
      x_4^2 - 4 = 0

      Analyical solutions are (+/-2, 0, +/-2, +/-2).
      Assumably due to the initial relaxation the iteration always converges to the solution in the negative quadrant.
    */

    if (context.isNamed("SolverOne")) {
      for (int i = 0; i < 4; i++) {
        outValues[i] = inValues[i]; //only pushes solution through
      }
    } else {
      outValues[0] = 2 * inValues[0] * inValues[0] - inValues[1] * inValues[2] - 8.0 + inValues[0];
      outValues[1] = inValues[0] * inValues[0] * inValues[1] + 2.0 * inValues[0] * inValues[1] * inValues[2] + inValues[1] * inValues[2] * inValues[2] + inValues[1];
      outValues[2] = inValues[2] * inValues[2] - 4.0 + inValues[2];
      outValues[3] = inValues[3] * inValues[3] - 4.0 + inValues[3];
    }

    interface.writeBlockScalarData(writeDataID, 4, vertexIDs, outValues);
    interface.advance(1.0);

    if (interface.isActionRequired(precice::constants::actionReadIterationCheckpoint())) {
      interface.markActionFulfilled(precice::constants::actionReadIterationCheckpoint());
    }
    iterations++;
  }

  interface.finalize();

  //relative residual in config is 1e-7, so 2 orders of magnitude less strict
  BOOST_TEST(outValues[0] == -2.0, boost::test_tools::tolerance(1e-5));
  BOOST_TEST(outValues[1] == 0.0, boost::test_tools::tolerance(1e-5));
  BOOST_TEST(outValues[2] == -2.0, boost::test_tools::tolerance(1e-5));
  BOOST_TEST(outValues[3] == -2.0, boost::test_tools::tolerance(1e-5));

  // to exclude false or no convergence
  BOOST_TEST(iterations <= 20);
  BOOST_TEST(iterations >= 5);
}

/// tests for different QN settings if correct fixed point is reached mesh with empty partition
void runTestQNEmptyPartition(std::string const &config, TestContext const &context)
{
  std::string meshName, writeDataName, readDataName;

  if (context.isNamed("SolverOne")) {
    meshName      = "MeshOne";
    writeDataName = "Data1";
    readDataName  = "Data2";
  } else {
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    meshName      = "MeshTwo";
    writeDataName = "Data2";
    readDataName  = "Data1";
  }

  SolverInterface interface(context.name, config, context.rank, context.size);
  int             meshID      = interface.getMeshID(meshName);
  int             writeDataID = interface.getDataID(writeDataName, meshID);
  int             readDataID  = interface.getDataID(readDataName, meshID);

  int vertexIDs[4];

  // meshes for rank 0 and rank 1, we use matching meshes for both participants
  double positions0[8] = {1.0, 0.0, 1.0, 0.5, 1.0, 1.0, 1.0, 1.5};

  if (context.isNamed("SolverOne")) {
    // All mesh is on Master
    if (context.isMaster()) {
      interface.setMeshVertices(meshID, 4, positions0, vertexIDs);
    }
  } else {
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    // All mesh is on Slave
    if (not context.isMaster()) {
      interface.setMeshVertices(meshID, 4, positions0, vertexIDs);
    }
  }

  interface.initialize();
  double inValues[4]  = {0.0, 0.0, 0.0, 0.0};
  double outValues[4] = {0.0, 0.0, 0.0, 0.0};

  int iterations = 0;

  while (interface.isCouplingOngoing()) {
    if (interface.isActionRequired(precice::constants::actionWriteIterationCheckpoint())) {
      interface.markActionFulfilled(precice::constants::actionWriteIterationCheckpoint());
    }

    if (interface.isReadDataAvailable())
      if ((context.isNamed("SolverOne") and context.isMaster()) or
          (context.isNamed("SolverTwo") and (not context.isMaster()))) {
        interface.readBlockScalarData(readDataID, 4, vertexIDs, inValues);
      }

    /*
      Solves the following non-linear equations, which are extended to a fixed-point equation (simply +x)
      2 * x_1^2 - x_2 * x_3 - 8 = 0
      x_1^2 * x_2 + 2 * x_1 * x_2 * x_3 + x_2 * x_3^2 + x_2 = 0
      x_3^2 - 4 = 0
      x_4^2 - 4 = 0

      Analyical solutions are (+/-2, 0, +/-2, +/-2).
      Assumably due to the initial relaxation the iteration always converges to the solution in the negative quadrant.
    */

    if (context.isNamed("SolverOne")) {
      for (int i = 0; i < 4; i++) {
        outValues[i] = inValues[i]; //only pushes solution through
      }
    } else {
      outValues[0] = 2 * inValues[0] * inValues[0] - inValues[1] * inValues[2] - 8.0 + inValues[0];
      outValues[1] = inValues[0] * inValues[0] * inValues[1] + 2.0 * inValues[0] * inValues[1] * inValues[2] + inValues[1] * inValues[2] * inValues[2] + inValues[1];
      outValues[2] = inValues[2] * inValues[2] - 4.0 + inValues[2];
      outValues[3] = inValues[3] * inValues[3] - 4.0 + inValues[3];
    }

    if ((context.isNamed("SolverOne") and context.isMaster()) or
        (context.isNamed("SolverTwo") and (not context.isMaster()))) {
      interface.writeBlockScalarData(writeDataID, 4, vertexIDs, outValues);
    }
    interface.advance(1.0);

    if (interface.isActionRequired(precice::constants::actionReadIterationCheckpoint())) {
      interface.markActionFulfilled(precice::constants::actionReadIterationCheckpoint());
    }
    iterations++;
  }

  interface.finalize();

  //relative residual in config is 1e-7, so 2 orders of magnitude less strict
  if ((context.isNamed("SolverOne") and context.isMaster()) or
      (context.isNamed("SolverTwo") and (not context.isMaster()))) {
    BOOST_TEST(outValues[0] == -2.0, boost::test_tools::tolerance(1e-5));
    BOOST_TEST(outValues[1] == 0.0, boost::test_tools::tolerance(1e-5));
    BOOST_TEST(outValues[2] == -2.0, boost::test_tools::tolerance(1e-5));
    BOOST_TEST(outValues[3] == -2.0, boost::test_tools::tolerance(1e-5));

    // to exclude false or no convergence
    BOOST_TEST(iterations <= 20);
    BOOST_TEST(iterations >= 5);
  }
}

BOOST_AUTO_TEST_CASE(TestQN1)
{
  PRECICE_TEST("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks));
  // serial coupling, IQN-ILS, strict QR2 filter
  std::string config = _pathToTests + "QN1.xml";
  runTestQN(config, context);
}

BOOST_AUTO_TEST_CASE(TestQN1EmptyPartition)
{
  PRECICE_TEST("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks));
  // serial coupling, IQN-ILS, strict QR2 filter
  std::string config = _pathToTests + "QN1.xml";
  runTestQNEmptyPartition(config, context);
}

BOOST_AUTO_TEST_CASE(TestQN2)
{
  PRECICE_TEST("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks));
  // parallel coupling, IQN-ILS, strict QR2 filter
  std::string config = _pathToTests + "QN2.xml";
  runTestQN(config, context);
}

BOOST_AUTO_TEST_CASE(TestQN2EmptyPartition)
{
  PRECICE_TEST("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks));
  // parallel coupling, IQN-ILS, strict QR2 filter
  std::string config = _pathToTests + "QN2.xml";
  runTestQNEmptyPartition(config, context);
}

BOOST_AUTO_TEST_CASE(TestQN3)
{
  PRECICE_TEST("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks));
  // serial coupling, IQN-IMVJ (which is identical to IQN-ILS as only first timestep is considered), relaxed QR2 filter
  std::string config = _pathToTests + "QN3.xml";
  runTestQN(config, context);
}

BOOST_AUTO_TEST_CASE(TestQN3EmptyPartition)
{
  PRECICE_TEST("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks));
  // parallel coupling, IQN-ILS, strict QR2 filter
  std::string config = _pathToTests + "QN3.xml";
  runTestQNEmptyPartition(config, context);
}

/// Tests various distributed communication schemes.
void runTestDistributedCommunication(std::string const &config, TestContext const &context)
{
  std::string meshName;
  int         i1 = -1, i2 = -1; //indices for data and positions

  std::vector<Eigen::VectorXd> positions;
  std::vector<Eigen::VectorXd> data;
  std::vector<Eigen::VectorXd> expectedData;

  Eigen::Vector3d position;
  Eigen::Vector3d datum;

  for (int i = 0; i < 4; i++) {
    position[0] = i * 1.0;
    position[1] = 0.0;
    position[2] = 0.0;
    positions.push_back(position);
    datum[0] = i * 1.0;
    datum[1] = i * 1.0;
    datum[2] = 0.0;
    data.push_back(datum);
    datum[0] = i * 2.0 + 1.0;
    datum[1] = i * 2.0 + 1.0;
    datum[2] = 1.0;
    expectedData.push_back(datum);
  }

  if (context.isNamed("Fluid")) {
    meshName = "FluidMesh";
    if (context.isMaster()) {
      i1 = 0;
      i2 = 2;
    } else {
      i1 = 2;
      i2 = 4;
    }
  } else {
    meshName = "StructureMesh";
    if (context.isMaster()) {
      i1 = 0;
      i2 = 1;
    } else {
      i1 = 1;
      i2 = 4;
    }
  }

  SolverInterface precice(context.name, config, context.rank, context.size);
  int             meshID   = precice.getMeshID(meshName);
  int             forcesID = precice.getDataID("Forces", meshID);
  int             velocID  = precice.getDataID("Velocities", meshID);

  std::vector<int> vertexIDs;
  for (int i = i1; i < i2; i++) {
    int vertexID = precice.setMeshVertex(meshID, positions[i].data());
    vertexIDs.push_back(vertexID);
  }

  precice.initialize();

  if (context.isNamed("Fluid")) { //Fluid
    for (size_t i = 0; i < vertexIDs.size(); i++) {
      precice.writeVectorData(forcesID, vertexIDs[i], data[i + i1].data());
    }
  } else {
    BOOST_TEST(context.isNamed("Structure"));
    for (size_t i = 0; i < vertexIDs.size(); i++) {
      precice.readVectorData(forcesID, vertexIDs[i], data[i].data());
      data[i] = (data[i] * 2).array() + 1.0;
      precice.writeVectorData(velocID, vertexIDs[i], data[i].data());
    }
  }

  precice.advance(1.0);

  if (context.isNamed("Fluid")) { //Fluid
    for (size_t i = 0; i < vertexIDs.size(); i++) {
      precice.readVectorData(velocID, vertexIDs[i], data[i + i1].data());
      for (size_t d = 0; d < 3; d++) {
        BOOST_TEST(expectedData[i + i1][d] == data[i + i1][d]);
      }
    }
  }

  precice.finalize();
}

BOOST_AUTO_TEST_CASE(TestDistributedCommunicationP2PSockets)
{
  PRECICE_TEST("Fluid"_on(2_ranks), "Structure"_on(2_ranks));
  std::string config = _pathToTests + "point-to-point-sockets.xml";
  runTestDistributedCommunication(config, context);
}

BOOST_AUTO_TEST_CASE(TestDistributedCommunicationP2PMPI)
{
  PRECICE_TEST("Fluid"_on(2_ranks), "Structure"_on(2_ranks));
  std::string config = _pathToTests + "point-to-point-mpi.xml";
  runTestDistributedCommunication(config, context);
}

BOOST_AUTO_TEST_CASE(TestDistributedCommunicationGatherScatterMPI)
{
  PRECICE_TEST("Fluid"_on(2_ranks), "Structure"_on(2_ranks));
  std::string config = _pathToTests + "gather-scatter-mpi.xml";
  runTestDistributedCommunication(config, context);
}

BOOST_AUTO_TEST_CASE(TestBoundingBoxInitialization)
{
  PRECICE_TEST("Fluid"_on(2_ranks), "Structure"_on(2_ranks));

  std::vector<Eigen::Vector3d> positions;
  std::vector<Eigen::Vector3d> data;
  std::vector<Eigen::Vector3d> expectedData;

  Eigen::Vector3d position;
  Eigen::Vector3d datum;

  for (int i = 0; i < 4; i++) {
    position[0] = i * 1.0;
    position[1] = i * 0.1;
    position[2] = -i * 10.0;
    positions.push_back(position);
    datum[0] = i * 1.0;
    datum[1] = i * 2.0;
    datum[2] = i * 3.0;
    data.push_back(datum);
    datum[0] = i * 1.0;
    datum[1] = i * 2.0;
    datum[2] = i * 3.0;
    expectedData.push_back(datum);
  }

  int i1 = -1, i2 = -1; //indices for data and positions

  if (context.isNamed("Fluid")) {
    if (context.isMaster()) {
      i1 = 2;
      i2 = 4;
    } else {
      i1 = 0;
      i2 = 2;
    }
  } else {
    BOOST_TEST(context.isNamed("Structure"));
    // This partiticipant starts with negated data
    for (int i = 0; i < 4; i++) {
      data[i] = -data[i];
    }
    if (context.isMaster()) {
      i1 = 0;
      i2 = 2;
    } else {
      i1 = 2;
      i2 = 4;
    }
  }
  BOOST_REQUIRE(i1 >= 0);
  BOOST_REQUIRE(i2 >= 0);

  std::string     configFilename = _pathToTests + "BB-sockets-explicit-oneway.xml";
  SolverInterface precice(context.name, configFilename, context.rank, context.size);

  int meshID   = precice.getMeshID(context.name + "Mesh");
  int forcesID = precice.getDataID("Forces", meshID);

  std::vector<int> vertexIDs;
  for (int i = i1; i < i2; i++) {
    int vertexID = precice.setMeshVertex(meshID, positions[i].data());
    vertexIDs.push_back(vertexID);
  }

  precice.initialize();

  if (context.isNamed("Fluid")) {
    for (size_t i = 0; i < vertexIDs.size(); i++) {
      precice.writeVectorData(forcesID, vertexIDs[i], data[i + i1].data());
    }
  }

  precice.advance(1.0);

  if (context.isNamed("Structure")) {
    for (size_t i = 0; i < vertexIDs.size(); i++) {
      precice.readVectorData(forcesID, vertexIDs[i], data[i + i1].data());
      for (size_t d = 0; d < 3; d++) {
        BOOST_TEST(expectedData[i + i1][d] == data[i + i1][d]);
      }
    }
  }

  precice.finalize();
}

BOOST_AUTO_TEST_CASE(TestBoundingBoxInitializationTwoWay)
{
  PRECICE_TEST("Fluid"_on(2_ranks), "Structure"_on(2_ranks));

  std::vector<Eigen::Vector3d> positions;
  std::vector<Eigen::Vector3d> data;
  std::vector<Eigen::Vector3d> expectedData;

  Eigen::Vector3d position;
  Eigen::Vector3d datum;

  int i1 = -1, i2 = -1; //indices for data and positions

  if (context.isNamed("Fluid")) {
    if (context.isMaster()) {
      i1 = 0;
      i2 = 2;
    } else {
      i1 = 2;
      i2 = 4;
    }
    for (int i = 0; i < 4; i++) {
      position[0] = i * 1.0;
      position[1] = 0;
      position[2] = 0;
      positions.push_back(position);
      datum[0] = i * 1.0;
      datum[1] = i * 2.0;
      datum[2] = i * 3.0;
      data.push_back(datum);
      datum[0] = -i * 1.0;
      datum[1] = -i * 2.0;
      datum[2] = -i * 3.0;
      expectedData.push_back(datum);
    }
  }

  if (context.isNamed("Structure")) {
    if (context.isMaster()) {
      i1 = 2;
      i2 = 4;
    } else {
      i1 = 0;
      i2 = 2;
    }

    for (int i = 0; i < 4; i++) {
      position[0] = i * 1.0;
      position[1] = 0.0;
      position[2] = 0.0;
      positions.push_back(position);
      datum[0] = -1.0;
      datum[1] = -1.0;
      datum[2] = -1.0;
      data.push_back(datum);
      datum[0] = i * 1.0;
      datum[1] = i * 2.0;
      datum[2] = i * 3.0;
      expectedData.push_back(datum);
    }
  }

  std::string     configFilename = _pathToTests + "BB-sockets-explicit-twoway.xml";
  SolverInterface precice(context.name, configFilename, context.rank, context.size);

  int meshID       = precice.getMeshID(context.name + "Mesh");
  int forcesID     = precice.getDataID("Forces", meshID);
  int velocitiesID = precice.getDataID("Velocities", meshID);

  std::vector<int> vertexIDs;
  for (int i = i1; i < i2; i++) {
    int vertexID = precice.setMeshVertex(meshID, positions[i].data());
    vertexIDs.push_back(vertexID);
  }

  precice.initialize();

  if (context.isNamed("Fluid")) {
    for (size_t i = 0; i < vertexIDs.size(); i++) {
      precice.writeVectorData(forcesID, vertexIDs[i], data[i + i1].data());
    }
  }

  if (context.isNamed("Structure")) {
    for (size_t i = 0; i < vertexIDs.size(); i++) {
      precice.readVectorData(forcesID, vertexIDs[i], data[i + i1].data());
      for (size_t d = 0; d < 3; d++) {
        BOOST_TEST(expectedData[i + i1][d] == data[i + i1][d]);
      }
    }

    for (size_t j = 0; j < 4; j++) {
      data[j] = -data[j].array();
    }

    for (size_t i = 0; i < vertexIDs.size(); i++) {
      precice.writeVectorData(velocitiesID, vertexIDs[i], data[i + i1].data());
    }
  }

  precice.advance(1.0);

  if (context.isNamed("Fluid")) {
    for (size_t i = 0; i < vertexIDs.size(); i++) {
      precice.readVectorData(velocitiesID, vertexIDs[i], data[i + i1].data());
      for (size_t d = 0; d < 3; d++) {
        BOOST_TEST(expectedData[i + i1][d] == data[i + i1][d]);
      }
    }
  }

  precice.finalize();
}

/// This testcase is based on a bug documented in issue #371
BOOST_AUTO_TEST_CASE(NearestProjectionRePartitioning)
{
  PRECICE_TEST("FluidSolver"_on(3_ranks), "SolidSolver"_on(1_rank));
  std::string configFilename = _pathToTests + "np-repartitioning.xml";

  if (context.isNamed("FluidSolver")) {
    SolverInterface interface(context.name, configFilename, context.rank, context.size);

    if (context.isMaster()) {
      interface.initialize();
      interface.advance(1.0);
      interface.finalize();
    } else {
      const int meshID     = interface.getMeshID("CellCenters");
      const int dimensions = 3;
      BOOST_TEST(interface.getDimensions() == dimensions);

      const int                 numberOfVertices = 65;
      const double              yCoord           = 0.0;
      const double              zCoord           = 0.005;
      const std::vector<double> positions{
          0.00124795, yCoord, zCoord,
          0.00375646, yCoord, zCoord,
          0.00629033, yCoord, zCoord,
          0.00884982, yCoord, zCoord,
          0.0114352, yCoord, zCoord,
          0.0140467, yCoord, zCoord,
          0.0166846, yCoord, zCoord,
          0.0193492, yCoord, zCoord,
          0.0220407, yCoord, zCoord,
          0.0247594, yCoord, zCoord,
          0.0275056, yCoord, zCoord,
          0.0302796, yCoord, zCoord,
          0.0330816, yCoord, zCoord,
          0.0359119, yCoord, zCoord,
          0.0387709, yCoord, zCoord,
          0.0416588, yCoord, zCoord,
          0.0445758, yCoord, zCoord,
          0.0475224, yCoord, zCoord,
          0.0504987, yCoord, zCoord,
          0.0535051, yCoord, zCoord,
          0.0565419, yCoord, zCoord,
          0.0596095, yCoord, zCoord,
          0.062708, yCoord, zCoord,
          0.0658378, yCoord, zCoord,
          0.0689993, yCoord, zCoord,
          0.0721928, yCoord, zCoord,
          0.0754186, yCoord, zCoord,
          0.0786769, yCoord, zCoord,
          0.0819682, yCoord, zCoord,
          0.0852928, yCoord, zCoord,
          0.088651, yCoord, zCoord,
          0.0920431, yCoord, zCoord,
          0.0954695, yCoord, zCoord,
          0.0989306, yCoord, zCoord,
          0.102427, yCoord, zCoord,
          0.105958, yCoord, zCoord,
          0.109525, yCoord, zCoord,
          0.113128, yCoord, zCoord,
          0.116768, yCoord, zCoord,
          0.120444, yCoord, zCoord,
          0.124158, yCoord, zCoord,
          0.127909, yCoord, zCoord,
          0.131698, yCoord, zCoord,
          0.135525, yCoord, zCoord,
          0.139391, yCoord, zCoord,
          0.143296, yCoord, zCoord,
          0.147241, yCoord, zCoord,
          0.151226, yCoord, zCoord,
          0.15525, yCoord, zCoord,
          0.159316, yCoord, zCoord,
          0.163422, yCoord, zCoord,
          0.16757, yCoord, zCoord,
          0.17176, yCoord, zCoord,
          0.175993, yCoord, zCoord,
          0.180268, yCoord, zCoord,
          0.184586, yCoord, zCoord,
          0.188948, yCoord, zCoord,
          0.193354, yCoord, zCoord,
          0.197805, yCoord, zCoord,
          0.202301, yCoord, zCoord,
          0.206842, yCoord, zCoord,
          0.211429, yCoord, zCoord,
          0.216062, yCoord, zCoord,
          0.220742, yCoord, zCoord,
          0.22547, yCoord, zCoord};
      BOOST_TEST(numberOfVertices * dimensions == positions.size());
      std::vector<int> vertexIDs(numberOfVertices);
      interface.setMeshVertices(meshID, numberOfVertices, positions.data(), vertexIDs.data());
      interface.initialize();
      BOOST_TEST(impl(interface).mesh("Nodes").triangles().size() == 15);
      interface.advance(1.0);
      interface.finalize();
    }
  } else {
    BOOST_TEST(context.isNamed("SolidSolver"));
    SolverInterface interface(context.name, configFilename, context.rank, context.size);
    const int       meshID     = interface.getMeshID("Nodes");
    const int       dimensions = 3;
    BOOST_TEST(interface.getDimensions() == dimensions);
    const int                 numberOfVertices = 34;
    const double              yCoord           = 0.0;
    const double              zCoord1          = 0.0;
    const double              zCoord2          = 0.01;
    const std::vector<double> positions{
        0.0, yCoord, zCoord2,
        0.0, yCoord, zCoord1,
        0.03125, yCoord, zCoord2,
        0.03125, yCoord, zCoord1,
        0.0625, yCoord, zCoord2,
        0.0625, yCoord, zCoord1,
        0.09375, yCoord, zCoord2,
        0.09375, yCoord, zCoord1,
        0.125, yCoord, zCoord2,
        0.125, yCoord, zCoord1,
        0.15625, yCoord, zCoord2,
        0.15625, yCoord, zCoord1,
        0.1875, yCoord, zCoord2,
        0.1875, yCoord, zCoord1,
        0.21875, yCoord, zCoord2,
        0.21875, yCoord, zCoord1,
        0.25, yCoord, zCoord2,
        0.25, yCoord, zCoord1,
        0.28125, yCoord, zCoord2,
        0.28125, yCoord, zCoord1,
        0.3125, yCoord, zCoord2,
        0.3125, yCoord, zCoord1,
        0.34375, yCoord, zCoord2,
        0.34375, yCoord, zCoord1,
        0.375, yCoord, zCoord2,
        0.375, yCoord, zCoord1,
        0.40625, yCoord, zCoord2,
        0.40625, yCoord, zCoord1,
        0.4375, yCoord, zCoord2,
        0.4375, yCoord, zCoord1,
        0.46875, yCoord, zCoord2,
        0.46875, yCoord, zCoord1,
        0.5, yCoord, zCoord2,
        0.5, yCoord, zCoord1};
    BOOST_TEST(numberOfVertices * dimensions == positions.size());
    std::vector<int> vertexIDs(numberOfVertices);
    interface.setMeshVertices(meshID, numberOfVertices, positions.data(), vertexIDs.data());

    const int        numberOfCells = numberOfVertices / 2 - 1;
    const int        numberOfEdges = numberOfCells * 4 + 1;
    std::vector<int> edgeIDs(numberOfEdges);

    for (int i = 0; i < numberOfCells; i++) {
      edgeIDs.at(4 * i)     = interface.setMeshEdge(meshID, vertexIDs.at(i * 2), vertexIDs.at(i * 2 + 1));     //left
      edgeIDs.at(4 * i + 1) = interface.setMeshEdge(meshID, vertexIDs.at(i * 2), vertexIDs.at(i * 2 + 2));     //top
      edgeIDs.at(4 * i + 2) = interface.setMeshEdge(meshID, vertexIDs.at(i * 2 + 1), vertexIDs.at(i * 2 + 3)); //bottom
      edgeIDs.at(4 * i + 3) = interface.setMeshEdge(meshID, vertexIDs.at(i * 2), vertexIDs.at(i * 2 + 3));     //diagonal
    }
    edgeIDs.at(numberOfEdges - 1) = interface.setMeshEdge(meshID, vertexIDs.at(numberOfVertices - 2), vertexIDs.at(numberOfVertices - 1)); //very right

    for (int i = 0; i < numberOfCells; i++) {
      interface.setMeshTriangle(meshID, edgeIDs.at(4 * i), edgeIDs.at(4 * i + 3), edgeIDs.at(4 * i + 2));     //left-diag-bottom
      interface.setMeshTriangle(meshID, edgeIDs.at(4 * i + 1), edgeIDs.at(4 * i + 3), edgeIDs.at(4 * i + 4)); //top-diag-right
    }

    interface.initialize();
    interface.advance(1.0);
    interface.finalize();
  }
}

BOOST_AUTO_TEST_CASE(MasterSockets)
{
  PRECICE_TEST("ParallelSolver"_on(3_ranks), "SerialSolver"_on(1_rank));
  std::string configFilename = _pathToTests + "master-sockets.xml";
  std::string myMeshName;
  if (context.isNamed("ParallelSolver")) {
    myMeshName = "ParallelMesh";
  } else {
    myMeshName = "SerialMesh";
  }
  SolverInterface interface(context.name, configFilename, context.rank, context.size);
  int             meshID      = interface.getMeshID(myMeshName);
  double          position[2] = {0, 0};
  interface.setMeshVertex(meshID, position);
  interface.initialize();
  interface.advance(1.0);
  interface.finalize();
}

// Tests SolverInterface() with a user-defined MPI communicator.
BOOST_AUTO_TEST_CASE(UserDefinedMPICommunicator)
{
  PRECICE_TEST("SolverOne"_on(3_ranks), "SolverTwo"_on(1_rank));
  std::string configFilename = _pathToTests + "userDefinedMPICommunicator.xml";

  if (context.isNamed("SolverOne")) {
    MPI_Comm        myComm = utils::Parallel::current()->comm;
    SolverInterface interface(context.name, configFilename, context.rank, context.size, &myComm);
    int             meshID = interface.getMeshID("MeshOne");

    int    vertexIDs[2];
    double xCoord       = context.rank * 0.4;
    double positions[4] = {xCoord, 0.0, xCoord + 0.2, 0.0};
    interface.setMeshVertices(meshID, 2, positions, vertexIDs);
    interface.initialize();
    interface.finalize();
  } else {
    SolverInterface interface(context.name, configFilename, context.rank, context.size);
    int             meshID = interface.getMeshID("MeshTwo");
    int             vertexIDs[6];
    double          positions[12] = {0.0, 0.0, 0.2, 0.0, 0.4, 0.0, 0.6, 0.0, 0.8, 0.0, 1.0, 0.0};
    interface.setMeshVertices(meshID, 6, positions, vertexIDs);
    interface.initialize();
    interface.finalize();
  }
}

#ifndef PRECICE_NO_PETSC

// Tests SolverInterface() with a user-defined MPI communicator.
// Since PETSc also uses MPI, we use petrbf mapping here.
BOOST_AUTO_TEST_CASE(UserDefinedMPICommunicatorPetRBF)
{
  PRECICE_TEST("SolverOne"_on(3_ranks), "SolverTwo"_on(1_rank));
  std::string           configFilename = _pathToTests + "userDefinedMPICommunicatorPetRBF.xml";
  config::Configuration config;

  if (context.isNamed("SolverOne")) {
    MPI_Comm myComm = utils::Parallel::current()->comm;

    SolverInterface interface(context.name, configFilename, context.rank, context.size, &myComm);
    int             meshID = interface.getMeshID("MeshOne");

    int    vertexIDs[2];
    double xCoord       = context.rank * 0.4;
    double positions[4] = {xCoord, 0.0, xCoord + 0.2, 0.0};
    interface.setMeshVertices(meshID, 2, positions, vertexIDs);
    interface.initialize();
    interface.finalize();
  } else {
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    SolverInterface interface(context.name, configFilename, context.rank, context.size);
    int             meshID = interface.getMeshID("MeshTwo");
    int             vertexIDs[6];
    double          positions[12] = {0.0, 0.0, 0.2, 0.0, 0.4, 0.0, 0.6, 0.0, 0.8, 0.0, 1.0, 0.0};
    interface.setMeshVertices(meshID, 6, positions, vertexIDs);
    interface.initialize();
    interface.finalize();
  }
}
#endif // PRECICE_NO_PETSC

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
#endif // PRECICE_NO_MPI

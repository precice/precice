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

std::string pathToTests = testing::getPathToSources() + "/precice/tests/gradient-tests/";

BOOST_AUTO_TEST_SUITE(PreciceTests)

BOOST_AUTO_TEST_SUITE(GradientMappingTests)

BOOST_AUTO_TEST_SUITE(SerialGradientMappingTests)

// Unidirectional Nearest Neighbor Gradient Read Mapping
BOOST_AUTO_TEST_CASE(NNG_Unidirectional_Read_Only)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank))
  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, pathToTests + "nng-unidirectional.xml", 0, 1);
  if (context.isNamed("A")) {

    int      meshOneID = cplInterface.getMeshID("MeshA");
    Vector3d posOne    = Vector3d::Constant(0.0);
    cplInterface.setMeshVertex(meshOneID, posOne.data());
    int dataID = cplInterface.getDataID("DataA", meshOneID);

    // Initialize, thus sending the mesh.
    double maxDt = cplInterface.initialize();
    BOOST_TEST(cplInterface.isCouplingOngoing(), "Sending participant should have to advance once!");

    double valueA = 1.0;
    cplInterface.writeScalarData(dataID, 0, valueA);
    cplInterface.writeScalarGradientData(dataID, 0, 1.0, 1.0, 1.0);

    // Participant must make move after writing
    maxDt = cplInterface.advance(maxDt);

    BOOST_TEST(!cplInterface.isCouplingOngoing(), "Sending participant should have to advance once!");
    cplInterface.finalize();

  } else {
    BOOST_TEST(context.isNamed("B"));
    int      meshTwoID = cplInterface.getMeshID("MeshB");
    Vector3d pos       = Vector3d::Constant(0.1);
    cplInterface.setMeshVertex(meshTwoID, pos.data());

    double maxDt = cplInterface.initialize();
    BOOST_TEST(cplInterface.isCouplingOngoing(), "Receiving participant should have to advance once!");

    int    dataID = cplInterface.getDataID("DataA", meshTwoID);
    double valueData;
    cplInterface.readScalarData(dataID, 0, valueData);
    BOOST_TEST(valueData == 1.3);

    cplInterface.advance(maxDt);
    BOOST_TEST(!cplInterface.isCouplingOngoing(), "Receiving participant should have to advance once!");

    cplInterface.finalize();
  }
}

// Bidirectional test : Read: Vector & NN - Write: Scalar & NNG (Serial coupling)
BOOST_AUTO_TEST_CASE(NNG_Bidirectional_Write_Scalar)
{

  //precice.isActionRequired(precice::constants::actionWriteInitialData()
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, pathToTests + "nng-write-scalar.xml", 0, 1);
  if (context.isNamed("SolverOne")) {
    int      meshOneID = cplInterface.getMeshID("MeshOne");
    Vector3d vec1      = Vector3d::Constant(0.1);
    cplInterface.setMeshVertex(meshOneID, vec1.data());
    double maxDt   = cplInterface.initialize();
    int    dataAID = cplInterface.getDataID("DataOne", meshOneID);
    int    dataBID = cplInterface.getDataID("DataTwo", meshOneID);

    double valueDataB = 0.0;
    cplInterface.initializeData();
    cplInterface.readScalarData(dataBID, 0, valueDataB);
    BOOST_TEST(1.3 == valueDataB);

    while (cplInterface.isCouplingOngoing()) {
      Vector3d valueDataA(1.0, 1.0, 1.0);
      cplInterface.writeVectorData(dataAID, 0, valueDataA.data());
      maxDt = cplInterface.advance(maxDt);

      cplInterface.readScalarData(dataBID, 0, valueDataB);
      BOOST_TEST(1.8 == valueDataB);
    }
    cplInterface.finalize();

  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    int      meshTwoID = cplInterface.getMeshID("MeshTwo");
    Vector3d vec2      = Vector3d::Constant(0.0);
    cplInterface.setMeshVertex(meshTwoID, vec2.data());

    double maxDt   = cplInterface.initialize();
    int    dataAID = cplInterface.getDataID("DataOne", meshTwoID);
    int    dataBID = cplInterface.getDataID("DataTwo", meshTwoID);

    double valueDataB = 1.0;
    cplInterface.writeScalarData(dataBID, 0, valueDataB);
    cplInterface.writeScalarGradientData(dataBID, 0, 1.0, 1.0, 1.0);

    //tell preCICE that data has been written and call initializeData
    cplInterface.markActionFulfilled(precice::constants::actionWriteInitialData());
    cplInterface.initializeData();

    Vector3d valueDataA;
    cplInterface.readVectorData(dataAID, 0, valueDataA.data());
    Vector3d expected(1.0, 1.0, 1.0);
    BOOST_TEST(valueDataA == expected);

    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(dataBID, 0, 1.5);
      cplInterface.writeScalarGradientData(dataBID, 0, 1.0, 1.0, 1.0);

      maxDt = cplInterface.advance(maxDt);
      cplInterface.readVectorData(dataAID, 0, valueDataA.data());
      BOOST_TEST(valueDataA == expected);
    }
    cplInterface.finalize();
  }
}

// Read : NN & Vector - Write : NNG & Vector (Parallel Coupling Scheme)
BOOST_AUTO_TEST_CASE(NNG_Bidirectional_Write_Vector)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))

  using Eigen::Vector2d;
  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, pathToTests + "nng-write-vector.xml", 0, 1);
  if (context.isNamed("SolverOne")) {
    int      meshOneID = cplInterface.getMeshID("MeshOne");
    Vector3d posOne    = Vector3d::Constant(0.0);
    cplInterface.setMeshVertex(meshOneID, posOne.data());
    double maxDt   = cplInterface.initialize();
    int    dataAID = cplInterface.getDataID("DataOne", meshOneID);
    int    dataBID = cplInterface.getDataID("DataTwo", meshOneID);

    Vector3d valueDataA(1.0, 1.0, 1.0);
    cplInterface.writeVectorData(dataAID, 0, valueDataA.data());

    cplInterface.markActionFulfilled(precice::constants::actionWriteInitialData());
    cplInterface.initializeData();

    Vector2d valueDataB;
    cplInterface.readVectorData(dataBID, 0, valueDataB.data());
    Vector2d expected(-1.0, 0.0);
    BOOST_TEST(valueDataB == expected);

    while (cplInterface.isCouplingOngoing()) {
      Vector3d valueDataA(2.0, 2.0, 2.0);
      cplInterface.writeVectorData(dataAID, 0, valueDataA.data());

      maxDt = cplInterface.advance(maxDt);

      cplInterface.readVectorData(dataBID, 0, valueDataB.data());
      expected << -0.5, 0.5;
      BOOST_TEST(valueDataB == expected);
    }
    cplInterface.finalize();

  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    int      meshTwoID = cplInterface.getMeshID("MeshTwo");
    Vector3d pos       = Vector3d::Constant(1.0);
    cplInterface.setMeshVertex(meshTwoID, pos.data());

    double maxDt   = cplInterface.initialize();
    int    dataAID = cplInterface.getDataID("DataOne", meshTwoID);
    int    dataBID = cplInterface.getDataID("DataTwo", meshTwoID);

    Vector2d valueDataB(2.0, 3.0);
    Vector2d gradient(1.0, 1.0);
    cplInterface.writeVectorData(dataBID, 0, valueDataB.data());
    cplInterface.writeVectorGradientData(dataBID, 0, gradient.data(), gradient.data(), gradient.data());

    //tell preCICE that data has been written and call initializeData
    cplInterface.markActionFulfilled(precice::constants::actionWriteInitialData());
    cplInterface.initializeData();

    Vector3d valueDataA;
    cplInterface.readVectorData(dataAID, 0, valueDataA.data());
    Vector3d expected(1.0, 1.0, 1.0);
    BOOST_TEST(valueDataA == expected);

    while (cplInterface.isCouplingOngoing()) {

      valueDataB << 2.5, 3.5;
      cplInterface.writeVectorData(dataBID, 0, valueDataB.data());
      cplInterface.writeVectorGradientData(dataBID, 0, gradient.data(), gradient.data(), gradient.data());

      maxDt = cplInterface.advance(maxDt);
      cplInterface.readVectorData(dataAID, 0, valueDataA.data());
      expected << 2.0, 2.0, 2.0;
      BOOST_TEST(valueDataA == expected);
    }
    cplInterface.finalize();
  }
}

// Bidirectional test : Read: Vector & NNG - Write: Vector & NN (Serial coupling)
BOOST_AUTO_TEST_CASE(NNG_Bidirectional_Read_Vector)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))

  using Eigen::Vector2d;
  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, pathToTests + "nng-read-vector.xml", 0, 1);
  if (context.isNamed("SolverOne")) {
    int      meshOneID = cplInterface.getMeshID("MeshOne");
    Vector3d posOne    = Vector3d::Constant(0.0);
    cplInterface.setMeshVertex(meshOneID, posOne.data());
    double maxDt   = cplInterface.initialize();
    int    dataAID = cplInterface.getDataID("DataOne", meshOneID);
    int    dataBID = cplInterface.getDataID("DataTwo", meshOneID);

    Vector2d valueDataB;

    cplInterface.markActionFulfilled(precice::constants::actionWriteInitialData());
    cplInterface.initializeData();
    cplInterface.readVectorData(dataBID, 0, valueDataB.data());
    Vector2d expected(2.0, 3.0);
    BOOST_TEST(valueDataB == expected);

    while (cplInterface.isCouplingOngoing()) {
      Vector3d valueDataA(1.0, 1.0, 1.0);
      Vector3d gradient(1.0, 1.0, 1.0);
      cplInterface.writeVectorData(dataAID, 0, valueDataA.data());
      cplInterface.writeVectorGradientData(dataAID, 0, gradient.data(), gradient.data(), gradient.data());

      maxDt = cplInterface.advance(maxDt);

      cplInterface.readVectorData(dataBID, 0, valueDataB.data());
      expected << 2.5, 3.5;
      BOOST_TEST(valueDataB == expected);
    }
    cplInterface.finalize();

  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    int      meshTwoID = cplInterface.getMeshID("MeshTwo");
    Vector3d pos       = Vector3d::Constant(1.0);
    cplInterface.setMeshVertex(meshTwoID, pos.data());

    double maxDt   = cplInterface.initialize();
    int    dataAID = cplInterface.getDataID("DataOne", meshTwoID);
    int    dataBID = cplInterface.getDataID("DataTwo", meshTwoID);

    Vector2d valueDataB(2.0, 3.0);
    cplInterface.writeVectorData(dataBID, 0, valueDataB.data());

    //tell preCICE that data has been written and call initializeData
    cplInterface.markActionFulfilled(precice::constants::actionWriteInitialData());
    cplInterface.initializeData();

    Vector3d valueDataA;
    cplInterface.readVectorData(dataAID, 0, valueDataA.data());
    Vector3d expected(4.0, 4.0, 4.0);
    BOOST_TEST(valueDataA == expected);

    while (cplInterface.isCouplingOngoing()) {

      valueDataB << 2.5, 3.5;
      cplInterface.writeVectorData(dataBID, 0, valueDataB.data());

      maxDt = cplInterface.advance(maxDt);
      cplInterface.readVectorData(dataAID, 0, valueDataA.data());
      BOOST_TEST(valueDataA == expected);
    }
    cplInterface.finalize();
  }
}

// Bidirectional test : Read: Vector & NN - Write: Scalar & NNG (Serial coupling)
BOOST_AUTO_TEST_CASE(NNG_Bidirectional_Read_Scalar)
{

  //precice.isActionRequired(precice::constants::actionWriteInitialData()
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, pathToTests + "nng-read-scalar.xml", 0, 1);
  if (context.isNamed("SolverOne")) {
    int      meshOneID = cplInterface.getMeshID("MeshOne");
    Vector3d vec1      = Vector3d::Constant(0.1);
    cplInterface.setMeshVertex(meshOneID, vec1.data());
    double maxDt   = cplInterface.initialize();
    int    dataAID = cplInterface.getDataID("DataOne", meshOneID);
    int    dataBID = cplInterface.getDataID("DataTwo", meshOneID);

    double valueDataB = 0.0;
    cplInterface.initializeData();
    cplInterface.readScalarData(dataBID, 0, valueDataB);
    BOOST_TEST(1.0 == valueDataB);

    while (cplInterface.isCouplingOngoing()) {

      cplInterface.writeScalarData(dataAID, 0, 3.0);
      cplInterface.writeScalarGradientData(dataAID, 0, 1.0, 2.0, 3.0);

      maxDt = cplInterface.advance(maxDt);

      cplInterface.readScalarData(dataBID, 0, valueDataB);
      BOOST_TEST(1.5 == valueDataB);
    }
    cplInterface.finalize();

  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    int      meshTwoID = cplInterface.getMeshID("MeshTwo");
    Vector3d vec2      = Vector3d::Constant(0.0);
    cplInterface.setMeshVertex(meshTwoID, vec2.data());

    double maxDt   = cplInterface.initialize();
    int    dataAID = cplInterface.getDataID("DataOne", meshTwoID);
    int    dataBID = cplInterface.getDataID("DataTwo", meshTwoID);

    double valueDataB = 1.0;
    cplInterface.writeScalarData(dataBID, 0, valueDataB);

    //tell preCICE that data has been written and call initializeData
    cplInterface.markActionFulfilled(precice::constants::actionWriteInitialData());
    cplInterface.initializeData();

    double valueDataA;
    cplInterface.readScalarData(dataAID, 0, valueDataA);
    BOOST_TEST(valueDataA == 2.4);

    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(dataBID, 0, 1.5);

      maxDt = cplInterface.advance(maxDt);
    }
    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(ParallelGradientMappingTests)

// In order to test enforced gather scatter communication with an empty master rank (see below)
void runTestEnforceGatherScatter(std::vector<double> masterPartition, std::string configFile)
{
  PRECICE_TEST("ParallelSolver"_on(2_ranks), "SerialSolver"_on(1_rank));
  std::string configFilename = configFile;

  if (context.isNamed("ParallelSolver")) {
    // Get mesh and data IDs
    SolverInterface interface(context.name, configFilename, context.rank, context.size);
    const int       meshID      = interface.getMeshID("ParallelMesh");
    const int       writeDataID = interface.getDataID("MyData1", meshID);
    const int       readDataID  = interface.getDataID("MyData2", meshID);
    const int       dim         = interface.getDimensions();
    BOOST_TEST(dim == 2);

    // Set coordinates, master according to input argument
    const std::vector<double> coordinates = context.isMaster() ? masterPartition : std::vector<double>{0.0, 0.5, 0.0, 3.5, 0.0, 5.0};
    const unsigned int        size        = coordinates.size() / dim;
    std::vector<int>          ids(size, 0);

    // Set mesh vertices
    interface.setMeshVertices(meshID, size, coordinates.data(), ids.data());

    // Initialize the solverinterface
    double dt = interface.initialize();

    // Create some dummy writeData
    std::vector<double> writeData;
    std::vector<double> writeGradientData;
    for (unsigned int i = 0; i < size; ++i) {
      writeData.emplace_back(i + 1);
      writeGradientData.emplace_back(0.0);
    }

    // Allocate memory for readData
    std::vector<double> readData(size);
    while (interface.isCouplingOngoing()) {
      // Write data, advance the solverinterface and readData
      interface.writeBlockScalarData(writeDataID, size,
                                     ids.data(), writeData.data());
      interface.writeBlockScalarGradientData(writeDataID, size,
                                             ids.data(), writeGradientData.data(), writeGradientData.data());

      dt = interface.advance(dt);
      interface.readBlockScalarData(readDataID, size,
                                    ids.data(), readData.data());
      // The received data on the slave rank is always the same
      if (!context.isMaster()) {
        BOOST_TEST(readData == std::vector<double>({3.4, 5.7, 4.0}));
      }
    }
  } else {
    // The serial participant
    BOOST_REQUIRE(context.isNamed("SerialSolver"));
    SolverInterface interface(context.name, configFilename, context.rank, context.size);
    // Get IDs
    const MeshID meshID      = interface.getMeshID("SerialMesh");
    const int    writeDataID = interface.getDataID("MyData2", meshID);
    const int    readDataID  = interface.getDataID("MyData1", meshID);
    const int    dim         = interface.getDimensions();
    BOOST_TEST(interface.getDimensions() == 2);

    // Define the interface
    const std::vector<double> coordinates{0.0, 0.5, 0.0, 3.5, 0.0, 5.0};
    const unsigned int        size = coordinates.size() / dim;
    std::vector<int>          ids(size);

    // Set vertices
    interface.setMeshVertices(meshID, size, coordinates.data(), ids.data());

    // Initialize the solverinterface
    double dt = interface.initialize();

    // Somce arbitrary write data
    std::vector<double> writeData{3.4, 5.7, 4.0};
    std::vector<double> writeGradientData{0.0, 0.0, 0.0};
    std::vector<double> readData(size);

    // Start the time loop
    while (interface.isCouplingOngoing()) {
      // Write data, advance solverinterface and read data
      interface.writeBlockScalarData(writeDataID, size,
                                     ids.data(), writeData.data());
      interface.writeBlockScalarGradientData(writeDataID, size,
                                             ids.data(), writeGradientData.data(), writeGradientData.data());
      dt = interface.advance(dt);
      interface.readBlockScalarData(readDataID, size,
                                    ids.data(), readData.data());
      // The received data is always the same
      if (!context.isMaster()) {
        BOOST_TEST(readData == std::vector<double>({1, 2, 3}));
      }
    }
  }
}

// Test case for an enforced gather scatter communication, where the partition
// on the master rank is empty (recieved and provided). See issue #1013 for details.
BOOST_AUTO_TEST_CASE(EnforceGatherScatterEmptyMaster)
{
  // Provided master partition is empty and received master partition is empty
  runTestEnforceGatherScatter(std::vector<double>{}, pathToTests + "nng-enforce-gather-scatter.xml");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
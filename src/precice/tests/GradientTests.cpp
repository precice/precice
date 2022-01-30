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

std::string pathToTests = testing::getPathToSources() + "/precice/tests/";

BOOST_AUTO_TEST_SUITE(PreciceTests)

BOOST_AUTO_TEST_SUITE(GradientMappingTests)

BOOST_AUTO_TEST_SUITE(SerialGradientMappingTests)

// Check for correct mesh requirements

BOOST_AUTO_TEST_CASE(NNG_Unidirectional_Parallel_A)
{
  PRECICE_TEST(1_rank);
  std::string filename = pathToTests + "nng-unidirectional-parallel.xml";

  SolverInterface interfaceA("A", filename, 0, 1);
  auto            meshIDA = interfaceA.getMeshID("MeshA");
  BOOST_TEST(interfaceA.isGradientRequired(meshIDA));
}

BOOST_AUTO_TEST_CASE(NNG_Unidirectional_Parallel_B)
{
  PRECICE_TEST(1_rank);
  std::string     filename = pathToTests + "nng-unidirectional-parallel.xml";
  SolverInterface interfaceB("B", filename, 0, 1);
  auto            meshIDB = interfaceB.getMeshID("MeshB");
  BOOST_TEST(!interfaceB.isGradientRequired(meshIDB));
}

BOOST_AUTO_TEST_CASE(NNG_Unidirectional_Serial_A)
{
  PRECICE_TEST(1_rank);
  std::string filename = pathToTests + "nng-unidirectional-serial.xml";

  SolverInterface interfaceA("A", filename, 0, 1);
  auto            meshIDA = interfaceA.getMeshID("MeshA");
  BOOST_TEST(interfaceA.isGradientRequired(meshIDA));
}

BOOST_AUTO_TEST_CASE(NNG_Unidirectional_Serial_B)
{
  PRECICE_TEST(1_rank);
  std::string     filename = pathToTests + "nng-unidirectional-serial.xml";
  SolverInterface interfaceB("B", filename, 0, 1);
  auto            meshIDB = interfaceB.getMeshID("MeshB");
  BOOST_TEST(!interfaceB.isGradientRequired(meshIDB));
}

// Unidirectional Nearest Neighbor Gradient Read Mapping
BOOST_AUTO_TEST_CASE(NNG_Unidirectional_Serial_Read_Only)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank))
  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, pathToTests + "nng-unidirectional-serial.xml", 0, 1);
  if (context.isNamed("A")) {

    int      meshOneID = cplInterface.getMeshID("MeshA");
    Vector3d posOne    = Vector3d::Constant(0.1);
    cplInterface.setMeshVertex(meshOneID, posOne.data());
    int dataID = cplInterface.getDataID("DataA", meshOneID);

    // Initialize, thus sending the mesh.
    double maxDt = cplInterface.initialize();
    BOOST_TEST(cplInterface.isCouplingOngoing(), "Sending participant should have to advance once!");

    BOOST_TEST(cplInterface.isGradientRequired(meshOneID));

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
    Vector3d pos       = Vector3d::Constant(0.0);
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

// Serial coupling, Bidirectional test : Read: Vector & NN - Write: Vector & NNG
BOOST_AUTO_TEST_CASE(NNG_Bidirectional_Serial_Write_Vector)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))

  using Eigen::Vector2d;
  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, pathToTests + "nng-write-serial-vector.xml", 0, 1);
  if (context.isNamed("SolverOne")) {
    int      meshOneID = cplInterface.getMeshID("MeshOne");
    Vector3d posOne    = Vector3d::Constant(1.0);
    cplInterface.setMeshVertex(meshOneID, posOne.data());
    double maxDt   = cplInterface.initialize();
    int    dataAID = cplInterface.getDataID("DataOne", meshOneID);
    int    dataBID = cplInterface.getDataID("DataTwo", meshOneID);

    Vector2d valueDataB;

    cplInterface.markActionFulfilled(precice::constants::actionWriteInitialData());
    cplInterface.initializeData();
    cplInterface.readVectorData(dataBID, 0, valueDataB.data());
    Vector2d expected(-1.0, 0.0);
    BOOST_TEST(valueDataB == expected);

    while (cplInterface.isCouplingOngoing()) {
      Vector3d valueDataA(1.0, 1.0, 1.0);
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
    Vector3d pos       = Vector3d::Constant(0.0);
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
      BOOST_TEST(valueDataA == expected);
    }
    cplInterface.finalize();
  }
}

// Serial coupling, Bidirectional test : Read: Vector & NN - Write: Scalar & NNG
BOOST_AUTO_TEST_CASE(NNG_Bidirectional_Serial_Write_Scalar)
{

  //precice.isActionRequired(precice::constants::actionWriteInitialData()
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, pathToTests + "nng-write-serial-scalar.xml", 0, 1);
  if (context.isNamed("SolverOne")) {
    int      meshOneID = cplInterface.getMeshID("MeshOne");
    Vector3d vec1      = Vector3d::Constant(0.0);
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
    Vector3d vec2      = Vector3d::Constant(0.1);
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

// Serial coupling, Bidirectional test : Read: Vector & NNG - Write: Vector & NN
BOOST_AUTO_TEST_CASE(NNG_Bidirectional_Serial_Read_Vector)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))

  using Eigen::Vector2d;
  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, pathToTests + "nng-read-serial-vector.xml", 0, 1);
  if (context.isNamed("SolverOne")) {
    int      meshOneID = cplInterface.getMeshID("MeshOne");
    Vector3d posOne    = Vector3d::Constant(1.0);
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
    Vector3d pos       = Vector3d::Constant(0.0);
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

// Serial coupling, Bidirectional test : Read: Vector & NN - Write: Scalar & NNG
BOOST_AUTO_TEST_CASE(NNG_Bidirectional_Serial_Read_Scalar)
{

  //precice.isActionRequired(precice::constants::actionWriteInitialData()
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, pathToTests + "nng-read-serial-scalar.xml", 0, 1);
  if (context.isNamed("SolverOne")) {
    int      meshOneID = cplInterface.getMeshID("MeshOne");
    Vector3d vec1      = Vector3d::Constant(0.0);
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
    Vector3d vec2      = Vector3d::Constant(0.1);
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

// Parallel Coupling Scheme : Read : NN & Vector - Write : NNG & Vector
BOOST_AUTO_TEST_CASE(NNG_Bidirectional_Parallel_Explicit_Vector)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))

  using Eigen::Vector2d;
  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, pathToTests + "nng-write-parallel-vector.xml", 0, 1);
  if (context.isNamed("SolverOne")) {
    int      meshOneID = cplInterface.getMeshID("MeshOne");
    Vector3d posOne    = Vector3d::Constant(1.0);
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
    Vector3d pos       = Vector3d::Constant(0.0);
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

// Parallel Coupling Scheme : Read : NN & Vector - Write : NNG & Scalar
BOOST_AUTO_TEST_CASE(NNG_Bidirectional_Parallel_Explicit_Scalar)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, pathToTests + "nng-write-parallel-scalar.xml", 0, 1);
  if (context.isNamed("SolverOne")) {
    int      meshOneID = cplInterface.getMeshID("MeshOne");
    Vector3d vec1      = Vector3d::Constant(0.0);
    cplInterface.setMeshVertex(meshOneID, vec1.data());
    double maxDt   = cplInterface.initialize();
    int    dataAID = cplInterface.getDataID("DataOne", meshOneID);
    int    dataBID = cplInterface.getDataID("DataTwo", meshOneID);

    Vector3d valueDataA(1.0, 1.0, 1.0);
    cplInterface.writeVectorData(dataAID, 0, valueDataA.data());

    cplInterface.markActionFulfilled(precice::constants::actionWriteInitialData());
    cplInterface.initializeData();

    double valueDataB = 0.0;
    cplInterface.readScalarData(dataBID, 0, valueDataB);
    BOOST_TEST(1.3 == valueDataB);

    while (cplInterface.isCouplingOngoing()) {
      Vector3d valueDataA(2.0, 2.0, 2.0);
      cplInterface.writeVectorData(dataAID, 0, valueDataA.data());

      maxDt = cplInterface.advance(maxDt);

      cplInterface.readScalarData(dataBID, 0, valueDataB);
      BOOST_TEST(1.8 == valueDataB);
    }
    cplInterface.finalize();

  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    int      meshTwoID = cplInterface.getMeshID("MeshTwo");
    Vector3d vec2      = Vector3d::Constant(0.1);
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
    Vector3d expected(1.0, 1.0, 1.0); //MUST ADD INITIALIZE TO EXCHANGE DATA (OTHERWISE DOESNT WORK)
    BOOST_TEST(valueDataA == expected);

    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(dataBID, 0, 1.5);
      cplInterface.writeScalarGradientData(dataBID, 0, 1.0, 1.0, 1.0);

      maxDt = cplInterface.advance(maxDt);
      cplInterface.readVectorData(dataAID, 0, valueDataA.data());
      expected << 2.0, 2.0, 2.0;
      BOOST_TEST(valueDataA == expected);
    }
    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(ParallelGradientMappingTests)

//TODO: FIX COMMUNICATION !!

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

/*
// ! No gradient data send is zero !


void runTestDistributedCommunication(std::string const &config, TestContext const &context)
{
  std::string meshName;
  int         i1 = -1, i2 = -1; //indices for data and positions

  std::vector<Eigen::VectorXd> positions;
  std::vector<Eigen::VectorXd> data;
  std::vector<Eigen::VectorXd> expectedData;
  std::vector<Eigen::VectorXd> gradientData;

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

    //gradient
    //datum[0] = 0.0;
    //datum[1] = 0.0;
    //datum[2] = 0.0;
    //gradientData.push_back(datum);
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
    VertexID vertexID = precice.setMeshVertex(meshID, positions[i].data());
    vertexIDs.push_back(vertexID);
  }

  precice.initialize();

  if (context.isNamed("Fluid")) { //Fluid
    for (size_t i = 0; i < vertexIDs.size(); i++) {
      precice.writeVectorData(forcesID, vertexIDs[i], data[i + i1].data());
      precice.writeVectorGradientData(forcesID, vertexIDs[i], gradientData[i + i1].data(), gradientData[i + i1].data(), gradientData[i + i1].data());
    }
  } else {
    BOOST_TEST(context.isNamed("Structure"));
    for (size_t i = 0; i < vertexIDs.size(); i++) {
      precice.readVectorData(forcesID, vertexIDs[i], data[i].data());
      data[i] = (data[i] * 2).array() + 1.0;
      precice.writeVectorData(velocID, vertexIDs[i], data[i].data());
      //precice.writeVectorGradientData(velocID, vertexIDs[i], gradientData[i].data(), gradientData[i].data(), gradientData[i].data());
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

BOOST_AUTO_TEST_CASE(TestDistributedCommunicationP2PMPI)
{
  PRECICE_TEST("Fluid"_on(2_ranks), "Structure"_on(2_ranks));
  std::string config = pathToTests + "nng-parallel-read-write-point-to-point-mpi.xml";
  runTestDistributedCommunication(config, context);
}
 */

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
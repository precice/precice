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

std::string pathToTests = testing::getPathToSources() + "/tests/serial/parallel-serial-mapping-nearest-neighbor-gradient/";

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(SerialGradientMappingTests)

// Unidirectional Nearest Neighbor Gradient Read Mapping
BOOST_AUTO_TEST_CASE(NNG_Unidirectional_Read_Only)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank))
  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, pathToTests + "nng-unidirectional-scalar.xml", 0, 1);
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

    if (cplInterface.isDataGradientRequired(dataID)) {
      BOOST_TEST(cplInterface.isDataGradientRequired(dataID) == true);
      Vector3d valueGradDataA(1.0, 1.0, 1.0);
      cplInterface.writeScalarGradientData(dataID, 0, valueGradDataA.data());
    }

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

// Unidirectional Nearest Neighbor Gradient Read Mapping
BOOST_AUTO_TEST_CASE(NNG_Unidirectional_Read_Block_Scalar)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank))
  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, pathToTests + "nng-unidirectional-scalar.xml", 0, 1);
  if (context.isNamed("A")) {

    int meshOneID = cplInterface.getMeshID("MeshA");
    int dataID    = cplInterface.getDataID("DataA", meshOneID);

    Vector3d posOne = Vector3d::Constant(0.0);
    Vector3d posTwo = Vector3d::Constant(1.0);
    cplInterface.setMeshVertex(meshOneID, posOne.data());
    cplInterface.setMeshVertex(meshOneID, posTwo.data());

    // Initialize, thus sending the mesh.
    double maxDt = cplInterface.initialize();
    BOOST_TEST(cplInterface.isCouplingOngoing(), "Sending participant should have to advance once!");

    double values[2]  = {1.0, 2.0};
    int    indices[2] = {0, 1};
    cplInterface.writeBlockScalarData(dataID, 2, indices, values);

    if (cplInterface.isDataGradientRequired(dataID)) {
      BOOST_TEST(cplInterface.isDataGradientRequired(dataID) == true);

      double gradientValues[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
      cplInterface.writeBlockScalarGradientData(dataID, 2, indices, gradientValues, true);
    }

    // Participant must make move after writing
    maxDt = cplInterface.advance(maxDt);

    BOOST_TEST(!cplInterface.isCouplingOngoing(), "Sending participant should have to advance once!");
    cplInterface.finalize();

  } else {
    BOOST_TEST(context.isNamed("B"));
    int meshTwoID = cplInterface.getMeshID("MeshB");
    int dataID    = cplInterface.getDataID("DataA", meshTwoID);

    Vector3d posOne = Vector3d::Constant(0.1);
    Vector3d posTwo = Vector3d::Constant(1.1);
    cplInterface.setMeshVertex(meshTwoID, posOne.data());
    cplInterface.setMeshVertex(meshTwoID, posTwo.data());

    double maxDt = cplInterface.initialize();
    BOOST_TEST(cplInterface.isCouplingOngoing(), "Receiving participant should have to advance once!");

    double valueData[2];
    int    indices[2] = {0, 1};
    cplInterface.readBlockScalarData(dataID, 2, indices, valueData);
    double expected[2] = {1.9, 3.2};
    // without romMajor : double expected[2] = {1.6, 3.5};
    BOOST_TEST(valueData == expected);

    cplInterface.advance(maxDt);
    BOOST_TEST(!cplInterface.isCouplingOngoing(), "Receiving participant should have to advance once!");

    cplInterface.finalize();
  }
}

// Unidirectional Nearest Neighbor Gradient Read Mapping
BOOST_AUTO_TEST_CASE(NNG_Unidirectional_Read_Block_Vector)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank))
  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, pathToTests + "nng-unidirectional-vector.xml", 0, 1);
  if (context.isNamed("A")) {

    int meshOneID = cplInterface.getMeshID("MeshA");
    int dataID    = cplInterface.getDataID("DataA", meshOneID);

    Vector3d posOne = Vector3d::Constant(0.0);
    Vector3d posTwo = Vector3d::Constant(1.0);
    cplInterface.setMeshVertex(meshOneID, posOne.data());
    cplInterface.setMeshVertex(meshOneID, posTwo.data());

    // Initialize, thus sending the mesh.
    double maxDt = cplInterface.initialize();
    BOOST_TEST(cplInterface.isCouplingOngoing(), "Sending participant should have to advance once!");

    double values[6]  = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    int    indices[2] = {0, 1};
    cplInterface.writeBlockVectorData(dataID, 2, indices, values);

    if (cplInterface.isDataGradientRequired(dataID)) {
      BOOST_TEST(cplInterface.isDataGradientRequired(dataID) == true);

      double gradientValues[18] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
                                   10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0};
      cplInterface.writeBlockVectorGradientData(dataID, 2, indices, gradientValues);
    }

    // Participant must make move after writing
    maxDt = cplInterface.advance(maxDt);

    BOOST_TEST(!cplInterface.isCouplingOngoing(), "Sending participant should have to advance once!");
    cplInterface.finalize();

  } else {
    BOOST_TEST(context.isNamed("B"));
    int meshTwoID = cplInterface.getMeshID("MeshB");
    int dataID    = cplInterface.getDataID("DataA", meshTwoID);

    Vector3d posOne = Vector3d::Constant(0.1);
    Vector3d posTwo = Vector3d::Constant(1.1);
    cplInterface.setMeshVertex(meshTwoID, posOne.data());
    cplInterface.setMeshVertex(meshTwoID, posTwo.data());

    double maxDt = cplInterface.initialize();
    BOOST_TEST(cplInterface.isCouplingOngoing(), "Receiving participant should have to advance once!");

    double valueData[6];
    int    indices[2] = {0, 1};
    cplInterface.readBlockVectorData(dataID, 2, indices, valueData);
    // To test rowMajor = true parameter : double expected[6] = {3.1, 4.4, 5.7, 7.0, 8.3, 9.6};
    double expected[6] = {1.6, 3.5, 5.4, 7.3, 9.2, 11.1};
    BOOST_TEST(valueData == expected);

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
  using Eigen::Vector2d;
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

    double   valueDataB = 1.0;
    Vector3d valueGradDataB(1.0, 1.0, 1.0);
    cplInterface.writeScalarData(dataBID, 0, valueDataB);
    cplInterface.writeScalarGradientData(dataBID, 0, valueGradDataB.data());

    //tell preCICE that data has been written and call initializeData
    cplInterface.markActionFulfilled(precice::constants::actionWriteInitialData());
    cplInterface.initializeData();

    Vector3d valueDataA;
    cplInterface.readVectorData(dataAID, 0, valueDataA.data());
    Vector3d expected(1.0, 1.0, 1.0);
    BOOST_TEST(valueDataA == expected);

    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(dataBID, 0, 1.5);
      Vector3d valueGradDataA(1.0, 1.0, 1.0);
      cplInterface.writeScalarGradientData(dataBID, 0, valueGradDataA.data());

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

    Vector2d                    valueDataB(2.0, 3.0);
    Eigen::Matrix<double, 3, 3> gradient;
    gradient << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;
    cplInterface.writeVectorData(dataBID, 0, valueDataB.data());
    cplInterface.writeVectorGradientData(dataBID, 0, gradient.data());

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
      cplInterface.writeVectorGradientData(dataBID, 0, gradient.data());

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
      Vector3d                    valueDataA(1.0, 1.0, 1.0);
      Eigen::Matrix<double, 3, 3> gradient;
      gradient << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;
      cplInterface.writeVectorData(dataAID, 0, valueDataA.data());
      cplInterface.writeVectorGradientData(dataAID, 0, gradient.data());

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

// Bidirectional test : Read: Vector & NN - Write: Scalar & NNG (Parallel coupling)
BOOST_AUTO_TEST_CASE(NNG_Bidirectional_Read_Scalar)
{

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
    BOOST_TEST(valueDataB == 1.0);

    while (cplInterface.isCouplingOngoing()) {

      cplInterface.writeScalarData(dataAID, 0, 3.0);
      Vector3d valueGradDataA(1.0, 2.0, 3.0);
      cplInterface.writeScalarGradientData(dataAID, 0, valueGradDataA.data());

      maxDt = cplInterface.advance(maxDt);

      cplInterface.readScalarData(dataBID, 0, valueDataB);
      BOOST_TEST(valueDataB == 1.5);
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

    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(dataBID, 0, 1.5);

      maxDt = cplInterface.advance(maxDt);

      double valueDataA;
      cplInterface.readScalarData(dataAID, 0, valueDataA);
      BOOST_TEST(valueDataA == 2.4);
    }
    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
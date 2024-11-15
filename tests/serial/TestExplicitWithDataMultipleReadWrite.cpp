#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
/**
 * @brief Tests the reading and writing of data multiple times within one timestep.
 *
 * The first solver performs multiple consistent readings of data sent by the second solver.
 * The second solver performs multiple sendings, of which the last is expected by the first solver.
 */
BOOST_AUTO_TEST_CASE(TestExplicitWithDataMultipleReadWrite)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  using Eigen::Vector3d;

  precice::Participant cplInterface(context.name, context.config(), 0, 1);
  if (context.isNamed("SolverOne")) {
    auto            meshName = "MeshOne";
    int             size     = 1;
    Eigen::VectorXi vertexIDs(size);
    Eigen::VectorXd readDataA(size * 3);
    Eigen::VectorXd readDataB(size);
    Eigen::VectorXd readPositions(size * 3);
    vertexIDs[0] = cplInterface.setMeshVertex(meshName, readPositions);

    auto dataAID = "DataOne";
    auto dataBID = "DataTwo";
    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    // multiple readBlockScalarData
    cplInterface.readData(meshName, dataBID, vertexIDs, maxDt, readDataB);
    // expected data value received
    BOOST_TEST(3.0 == readDataB[0]);
    // change value at read destination
    readDataB[0] = -1.11;
    BOOST_TEST(-1.11 == readDataB[0]);
    cplInterface.readData(meshName, dataBID, vertexIDs, maxDt, readDataB);
    // expected data value received
    BOOST_TEST(3.0 == readDataB[0]);
    // change value at read destination
    readDataB[0] = -1.12;
    BOOST_TEST(-1.12 == readDataB[0]);

    // multiple readScalarData
    cplInterface.readData(meshName, dataBID, vertexIDs, maxDt, readDataB);
    // expected data value received
    BOOST_TEST(3.0 == readDataB[0]);
    // change value at read destination
    readDataB[0] = -1.21;
    BOOST_TEST(-1.21 == readDataB[0]);
    cplInterface.readData(meshName, dataBID, vertexIDs, maxDt, readDataB);
    // expected data value received
    BOOST_TEST(3.0 == readDataB[0]);

    // multiple readBlockVectorData
    cplInterface.readData(meshName, dataAID, vertexIDs, maxDt, readDataA);
    // expected data value received
    BOOST_TEST(Vector3d(7.0, 7.0, 7.0) == readDataA);
    // change value at read destination
    readDataA[0] = -1.31;
    readDataA[1] = -1.31;
    readDataA[2] = -1.31;
    BOOST_TEST(Vector3d(-1.31, -1.31, -1.31) == readDataA);
    cplInterface.readData(meshName, dataAID, vertexIDs, maxDt, readDataA);
    // expected data value received
    BOOST_TEST(Vector3d(7.0, 7.0, 7.0) == readDataA);
    // change value at read destination
    readDataA[0] = -1.32;
    readDataA[1] = -1.32;
    readDataA[2] = -1.32;
    BOOST_TEST(Vector3d(-1.32, -1.32, -1.32) == readDataA);

    // multiple readVectorData
    cplInterface.readData(meshName, dataAID, vertexIDs, maxDt, readDataA);
    // expected data value received
    BOOST_TEST(Vector3d(7.0, 7.0, 7.0) == readDataA);
    // change value at read destination
    readDataA[0] = -1.41;
    readDataA[1] = -1.41;
    readDataA[2] = -1.41;
    BOOST_TEST(Vector3d(-1.41, -1.41, -1.41) == readDataA);
    // multiple readVectorData
    cplInterface.readData(meshName, dataAID, vertexIDs, maxDt, readDataA);
    // expected data value received
    BOOST_TEST(Vector3d(7.0, 7.0, 7.0) == readDataA);

    while (cplInterface.isCouplingOngoing()) {
      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      // multiple readBlockScalarData
      cplInterface.readData(meshName, dataBID, vertexIDs, maxDt, readDataB);
      // expected data value received
      BOOST_TEST(5.0 == readDataB[0]);
      // change value at read destination
      readDataB[0] = -1.51;
      BOOST_TEST(-1.51 == readDataB[0]);
      cplInterface.readData(meshName, dataBID, vertexIDs, maxDt, readDataB);
      // expected data value received
      BOOST_TEST(5.0 == readDataB[0]);
      // change value at read destination
      readDataB[0] = -1.52;
      BOOST_TEST(-1.52 == readDataB[0]);

      // multiple readScalarData
      cplInterface.readData(meshName, dataBID, vertexIDs, maxDt, readDataB);
      // expected data value received
      BOOST_TEST(5.0 == readDataB[0]);
      // change value at read destination
      readDataB[0] = -1.61;
      BOOST_TEST(-1.61 == readDataB[0]);
      cplInterface.readData(meshName, dataBID, vertexIDs, maxDt, readDataB);
      // expected data value received
      BOOST_TEST(5.0 == readDataB[0]);

      // multiple readBlockVectorData
      cplInterface.readData(meshName, dataAID, vertexIDs, maxDt, readDataA);
      // expected data value received
      BOOST_TEST(Vector3d(9.0, 9.0, 9.0) == readDataA);
      // change value at read destination
      readDataA[0] = -1.71;
      readDataA[1] = -1.71;
      readDataA[2] = -1.71;
      BOOST_TEST(Vector3d(-1.71, -1.71, -1.71) == readDataA);
      cplInterface.readData(meshName, dataAID, vertexIDs, maxDt, readDataA);
      // expected data value received
      BOOST_TEST(Vector3d(9.0, 9.0, 9.0) == readDataA);
      // change value at read destination
      readDataA[0] = -1.72;
      readDataA[1] = -1.72;
      readDataA[2] = -1.72;
      BOOST_TEST(Vector3d(-1.72, -1.72, -1.72) == readDataA);

      // multiple readVectorData
      cplInterface.readData(meshName, dataAID, vertexIDs, maxDt, readDataA);
      // expected data value received
      BOOST_TEST(Vector3d(9.0, 9.0, 9.0) == readDataA);
      // change value at read destination
      readDataA[0] = -1.81;
      readDataA[1] = -1.81;
      readDataA[2] = -1.81;
      BOOST_TEST(Vector3d(-1.81, -1.81, -1.81) == readDataA);
      cplInterface.readData(meshName, dataAID, vertexIDs, maxDt, readDataA);
      // expected data value received
      BOOST_TEST(Vector3d(9.0, 9.0, 9.0) == readDataA);
    }
    cplInterface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    auto            meshName = "MeshTwo";
    int             size     = 1;
    Eigen::VectorXi vertexIDs(size);
    Eigen::VectorXd writeDataA(size * 3);
    Eigen::VectorXd writeDataB(size);
    Eigen::VectorXd writePositions(size * 3);
    vertexIDs[0] = cplInterface.setMeshVertex(meshName, writePositions);

    BOOST_REQUIRE(cplInterface.requiresInitialData());

    auto dataAID = "DataOne";
    auto dataBID = "DataTwo";

    // multiple writeBlockScalarData
    writeDataB[0] = -2.11;
    BOOST_TEST(-2.11 == writeDataB[0]);
    cplInterface.writeData(meshName, dataBID, vertexIDs, writeDataB);
    // changed data value sent, overwriting previous one
    writeDataB[0] = -2.12;
    BOOST_TEST(-2.12 == writeDataB[0]);
    cplInterface.writeData(meshName, dataBID, vertexIDs, writeDataB);
    // changed data value sent, overwriting previous one
    writeDataB[0] = -2.13;
    BOOST_TEST(-2.13 == writeDataB[0]);

    // multiple writeScalarData
    cplInterface.writeData(meshName, dataBID, vertexIDs, writeDataB);
    // expected data value sent, overwriting previous ones
    writeDataB[0] = 3.0;
    BOOST_TEST(3.0 == writeDataB[0]);
    cplInterface.writeData(meshName, dataBID, vertexIDs, writeDataB);

    // multiple writeBlockVectorData
    writeDataA[0] = -2.31;
    writeDataA[1] = -2.31;
    writeDataA[2] = -2.31;
    BOOST_TEST(Vector3d(-2.31, -2.31, -2.31) == writeDataA);
    cplInterface.writeData(meshName, dataAID, vertexIDs, writeDataA);
    // changed data value sent, overwriting previous one
    writeDataA[0] = -2.32;
    writeDataA[1] = -2.32;
    writeDataA[2] = -2.32;
    BOOST_TEST(Vector3d(-2.32, -2.32, -2.32) == writeDataA);
    cplInterface.writeData(meshName, dataAID, vertexIDs, writeDataA);
    // changed data value sent, overwriting previous one
    writeDataA[0] = -2.33;
    writeDataA[1] = -2.33;
    writeDataA[2] = -2.33;
    BOOST_TEST(Vector3d(-2.33, -2.33, -2.33) == writeDataA);

    // multiple writeVectorData
    cplInterface.writeData(meshName, dataAID, vertexIDs, writeDataA);
    // expected data value sent, overwriting previous one
    writeDataA[0] = 7.0;
    writeDataA[1] = 7.0;
    writeDataA[2] = 7.0;
    BOOST_TEST(Vector3d(7.0, 7.0, 7.0) == writeDataA);
    cplInterface.writeData(meshName, dataAID, vertexIDs, writeDataA);

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();
    while (cplInterface.isCouplingOngoing()) {
      // multiple writeBlockScalarData
      writeDataB[0] = -2.51;
      BOOST_TEST(-2.51 == writeDataB[0]);
      cplInterface.writeData(meshName, dataBID, vertexIDs, writeDataB);
      // changed data value sent, overwriting previous one
      writeDataB[0] = -2.52;
      BOOST_TEST(-2.52 == writeDataB[0]);
      cplInterface.writeData(meshName, dataBID, vertexIDs, writeDataB);
      // changed data value sent, overwriting previous one
      writeDataB[0] = -2.53;
      BOOST_TEST(-2.53 == writeDataB[0]);

      // multiple writeScalarData
      cplInterface.writeData(meshName, dataBID, vertexIDs, writeDataB);
      // expected data value sent, overwriting previous one
      writeDataB[0] = 5.0;
      BOOST_TEST(5.0 == writeDataB[0]);
      cplInterface.writeData(meshName, dataBID, vertexIDs, writeDataB);

      // multiple writeBlockVectorData
      writeDataA[0] = -2.71;
      writeDataA[1] = -2.71;
      writeDataA[2] = -2.71;
      BOOST_TEST(Vector3d(-2.71, -2.71, -2.71) == writeDataA);
      cplInterface.writeData(meshName, dataAID, vertexIDs, writeDataA);
      // changed data value sent, overwriting previous one
      writeDataA[0] = -2.72;
      writeDataA[1] = -2.72;
      writeDataA[2] = -2.72;
      BOOST_TEST(Vector3d(-2.72, -2.72, -2.72) == writeDataA);
      cplInterface.writeData(meshName, dataAID, vertexIDs, writeDataA);
      // changed data value sent, overwriting previous one
      writeDataA[0] = -2.73;
      writeDataA[1] = -2.73;
      writeDataA[2] = -2.73;
      BOOST_TEST(Vector3d(-2.73, -2.73, -2.73) == writeDataA);

      // multiple writeVectorData
      cplInterface.writeData(meshName, dataAID, vertexIDs, writeDataA);
      // expected data value sent, overwriting previous one
      writeDataA[0] = 9.0;
      writeDataA[1] = 9.0;
      writeDataA[2] = 9.0;
      BOOST_TEST(Vector3d(9.0, 9.0, 9.0) == writeDataA);
      cplInterface.writeData(meshName, dataAID, vertexIDs, writeDataA);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();
    }
    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial

#endif // PRECICE_NO_MPI

#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
/**
 * @brief Tests the reading of data using available API functions
 */
BOOST_AUTO_TEST_CASE(TestReadAPI)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  using Eigen::Vector3d;

  precice::Participant cplInterface(context.name, context.config(), 0, 1);
  if (context.isNamed("SolverOne")) {
    auto            meshName = "MeshOne";
    int             size     = 1;
    Eigen::VectorXi vertexIDs(size);
    Eigen::VectorXd writeDataA(size * 3);
    Eigen::VectorXd readDataB(size);
    Eigen::VectorXd readPositions(size * 3);
    vertexIDs[0] = cplInterface.setMeshVertex(meshName, readPositions);

    auto dataAID = "DataOne";
    auto dataBID = "DataTwo";

    BOOST_REQUIRE(cplInterface.requiresInitialData());

    // writeVectorData
    writeDataA[0] = 7.0;
    writeDataA[1] = 7.0;
    writeDataA[2] = 7.0;
    cplInterface.writeData(meshName, dataAID, vertexIDs, writeDataA);

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    // readBlockScalarData without waveform
    cplInterface.readData(meshName, dataBID, vertexIDs, maxDt, readDataB);
    // expected data value received
    BOOST_TEST(3.0 == readDataB[0]);
    // reset value at read destination
    readDataB[0] = 0;
    // readBlockScalarData with waveform
    cplInterface.readData(meshName, dataBID, vertexIDs, 0.5, readDataB);
    // expected data value received
    BOOST_TEST(3.0 == readDataB[0]);
    // reset value at read destination
    readDataB[0] = 0;

    // readScalarData without waveform
    cplInterface.readData(meshName, dataBID, vertexIDs, maxDt, readDataB);
    // expected data value received
    BOOST_TEST(3.0 == readDataB[0]);
    // reset value at read destination
    readDataB[0] = 0;
    // readScalarData with waveform
    cplInterface.readData(meshName, dataBID, vertexIDs, 0.5, readDataB);
    // expected data value received
    BOOST_TEST(3.0 == readDataB[0]);
    // reset value at read destination
    readDataB[0] = 0;

    while (cplInterface.isCouplingOngoing()) {
      if (cplInterface.requiresWritingCheckpoint()) {
      }
      if (cplInterface.requiresReadingCheckpoint()) {
      }

      // writeVectorData
      writeDataA[0] = 14.0;
      writeDataA[1] = 14.0;
      writeDataA[2] = 14.0;
      cplInterface.writeData(meshName, dataAID, vertexIDs, writeDataA);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      if (cplInterface.isCouplingOngoing()) {

        // readBlockScalarData without waveform
        cplInterface.readData(meshName, dataBID, vertexIDs, maxDt, readDataB);
        // expected data value received
        BOOST_TEST(6.0 == readDataB[0]);
        // reset value at read destination
        readDataB[0] = 0;
        // readBlockScalarData with waveform
        cplInterface.readData(meshName, dataBID, vertexIDs, 0.5, readDataB);
        // expected data value received
        BOOST_TEST(4.5 == readDataB[0]);
        // reset value at read destination
        readDataB[0] = 0;

        // readScalarData without waveform
        cplInterface.readData(meshName, dataBID, vertexIDs, maxDt, readDataB);
        // expected data value received
        BOOST_TEST(6.0 == readDataB[0]);
        // reset value at read destination
        readDataB[0] = 0;
        // readScalarData with waveform
        cplInterface.readData(meshName, dataBID, vertexIDs, 0.5, readDataB);
        // expected data value received
        BOOST_TEST(4.5 == readDataB[0]);
        // reset value at read destination
        readDataB[0] = 0;
      }
    }
    cplInterface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    auto            meshName = "MeshTwo";
    int             size     = 1;
    Eigen::VectorXi vertexIDs(size);
    Eigen::VectorXd readDataA(size * 3);
    Eigen::VectorXd writeDataB(size);
    Eigen::VectorXd writePositions(size * 3);
    vertexIDs[0] = cplInterface.setMeshVertex(meshName, writePositions);

    auto dataAID = "DataOne";
    auto dataBID = "DataTwo";

    BOOST_REQUIRE(cplInterface.requiresInitialData());

    // writeScalarData
    writeDataB[0] = 3.0;
    cplInterface.writeData(meshName, dataBID, vertexIDs, writeDataB);

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    // readBlockVectorData without waveform
    cplInterface.readData(meshName, dataAID, vertexIDs, maxDt, readDataA);
    // expected data value received
    BOOST_TEST(Vector3d(7.0, 7.0, 7.0) == readDataA);
    // reset value at read destination
    readDataA[0] = 0;
    readDataA[1] = 0;
    readDataA[2] = 0;
    // readBlockVectorData with waveform
    cplInterface.readData(meshName, dataAID, vertexIDs, 0.5, readDataA);
    // expected data value received
    BOOST_TEST(Vector3d(7.0, 7.0, 7.0) == readDataA);
    // reset value at read destination
    readDataA[0] = 0;
    readDataA[1] = 0;
    readDataA[2] = 0;

    // readVectorData without waveform
    cplInterface.readData(meshName, dataAID, vertexIDs, maxDt, readDataA);
    // expected data value received
    BOOST_TEST(Vector3d(7.0, 7.0, 7.0) == readDataA);
    // reset value at read destination
    readDataA[0] = 0;
    readDataA[1] = 0;
    readDataA[2] = 0;
    // readVectorData with waveform
    cplInterface.readData(meshName, dataAID, vertexIDs, 0.5, readDataA);
    // expected data value received
    BOOST_TEST(Vector3d(7.0, 7.0, 7.0) == readDataA);
    // reset value at read destination
    readDataA[0] = 0;
    readDataA[1] = 0;
    readDataA[2] = 0;

    while (cplInterface.isCouplingOngoing()) {
      if (cplInterface.requiresWritingCheckpoint()) {
      }
      if (cplInterface.requiresReadingCheckpoint()) {
      }

      // writeScalarData
      writeDataB[0] = 6.0;
      cplInterface.writeData(meshName, dataBID, vertexIDs, writeDataB);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      if (cplInterface.isCouplingOngoing()) {
        // readBlockVectorData without waveform
        cplInterface.readData(meshName, dataAID, vertexIDs, maxDt, readDataA);
        // expected data value received
        BOOST_TEST(Vector3d(14.0, 14.0, 14.0) == readDataA);
        // reset value at read destination
        readDataA[0] = 0;
        readDataA[1] = 0;
        readDataA[2] = 0;
        // readBlockVectorData with waveform
        cplInterface.readData(meshName, dataAID, vertexIDs, 0.5, readDataA);
        // expected data value received
        BOOST_TEST(Vector3d(10.5, 10.5, 10.5) == readDataA);
        // reset value at read destination
        readDataA[0] = 0;
        readDataA[1] = 0;
        readDataA[2] = 0;

        // readVectorData without waveform
        cplInterface.readData(meshName, dataAID, vertexIDs, maxDt, readDataA);
        // expected data value received
        BOOST_TEST(Vector3d(14.0, 14.0, 14.0) == readDataA);
        // reset value at read destination
        readDataA[0] = 0;
        readDataA[1] = 0;
        readDataA[2] = 0;
        // readVectorData with waveform
        cplInterface.readData(meshName, dataAID, vertexIDs, 0.5, readDataA);
        // expected data value received
        BOOST_TEST(Vector3d(10.5, 10.5, 10.5) == readDataA);
        // reset value at read destination
        readDataA[0] = 0;
        readDataA[1] = 0;
        readDataA[2] = 0;
      }
    }
    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial

#endif // PRECICE_NO_MPI

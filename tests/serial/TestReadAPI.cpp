#include <boost/test/tools/old/interface.hpp>
#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
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

  precice::SolverInterface cplInterface(context.name, context.config(), 0, 1);
  if (context.isNamed("SolverOne")) {
    auto            meshName = "MeshOne";
    int             size     = 1;
    Eigen::VectorXi vertexIDs(size);
    Eigen::VectorXd writeDataA(size * 3);
    Eigen::VectorXd readDataB(size);
    Eigen::VectorXd readPositions(size * 3);
    vertexIDs[0] = cplInterface.setMeshVertex(meshName, readPositions.data());

    auto dataAID = "DataOne"; //  meshOneID
    auto dataBID = "DataTwo"; //  meshOneID

    BOOST_REQUIRE(cplInterface.requiresInitialData());

    // writeVectorData
    writeDataA[0] = 7.0;
    writeDataA[1] = 7.0;
    writeDataA[2] = 7.0;
    cplInterface.writeVectorData(meshName, dataAID, vertexIDs[0], writeDataA.data());

    double maxDt = cplInterface.initialize();

    // readBlockScalarData without waveform
    cplInterface.readBlockScalarData(meshName, dataBID, 1, vertexIDs.data(), readDataB.data());
    // expected data value received
    BOOST_TEST(3.0 == readDataB[0]);
    // reset value at read destination
    readDataB[0] = 0;
    // readBlockScalarData with waveform
    cplInterface.readBlockScalarData(meshName, dataBID, 1, vertexIDs.data(), 0.5, readDataB.data());
    // expected data value received
    BOOST_TEST(3.0 == readDataB[0]);
    // reset value at read destination
    readDataB[0] = 0;

    // readScalarData without waveform
    cplInterface.readScalarData(meshName, dataBID, vertexIDs[0], readDataB[0]);
    // expected data value received
    BOOST_TEST(3.0 == readDataB[0]);
    // reset value at read destination
    readDataB[0] = 0;
    // readScalarData with waveform
    cplInterface.readScalarData(meshName, dataBID, vertexIDs[0], 0.5, readDataB[0]);
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
      cplInterface.writeVectorData(meshName, dataAID, vertexIDs[0], writeDataA.data());

      maxDt = cplInterface.advance(maxDt);

      if (cplInterface.isCouplingOngoing()) {

        // readBlockScalarData without waveform
        cplInterface.readBlockScalarData(meshName, dataBID, 1, vertexIDs.data(), readDataB.data());
        // expected data value received
        BOOST_TEST(6.0 == readDataB[0]);
        // reset value at read destination
        readDataB[0] = 0;
        // readBlockScalarData with waveform
        cplInterface.readBlockScalarData(meshName, dataBID, 1, vertexIDs.data(), 0.5, readDataB.data());
        // expected data value received
        BOOST_TEST(4.5 == readDataB[0]);
        // reset value at read destination
        readDataB[0] = 0;

        // readScalarData without waveform
        cplInterface.readScalarData(meshName, dataBID, vertexIDs[0], readDataB[0]);
        // expected data value received
        BOOST_TEST(6.0 == readDataB[0]);
        // reset value at read destination
        readDataB[0] = 0;
        // readScalarData with waveform
        cplInterface.readScalarData(meshName, dataBID, vertexIDs[0], 0.5, readDataB[0]);
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
    vertexIDs[0] = cplInterface.setMeshVertex(meshName, writePositions.data());

    auto dataAID = "DataOne"; //  meshTwoID
    auto dataBID = "DataTwo"; //  meshTwoID

    BOOST_REQUIRE(cplInterface.requiresInitialData());

    // writeScalarData
    writeDataB[0] = 3.0;
    cplInterface.writeScalarData(meshName, dataBID, vertexIDs[0], writeDataB[0]);

    double maxDt = cplInterface.initialize();

    // readBlockVectorData without waveform
    cplInterface.readBlockVectorData(meshName, dataAID, 1, vertexIDs.data(), readDataA.data());
    // expected data value received
    BOOST_TEST(Vector3d(7.0, 7.0, 7.0) == readDataA);
    // reset value at read destination
    readDataA[0] = 0;
    readDataA[1] = 0;
    readDataA[2] = 0;
    // readBlockVectorData with waveform
    cplInterface.readBlockVectorData(meshName, dataAID, 1, vertexIDs.data(), 0.5, readDataA.data());
    // expected data value received
    BOOST_TEST(Vector3d(7.0, 7.0, 7.0) == readDataA);
    // reset value at read destination
    readDataA[0] = 0;
    readDataA[1] = 0;
    readDataA[2] = 0;

    // readVectorData without waveform
    cplInterface.readVectorData(meshName, dataAID, vertexIDs[0], readDataA.data());
    // expected data value received
    BOOST_TEST(Vector3d(7.0, 7.0, 7.0) == readDataA);
    // reset value at read destination
    readDataA[0] = 0;
    readDataA[1] = 0;
    readDataA[2] = 0;
    // readVectorData with waveform
    cplInterface.readVectorData(meshName, dataAID, vertexIDs[0], 0.5, readDataA.data());
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
      cplInterface.writeScalarData(meshName, dataBID, vertexIDs[0], writeDataB[0]);

      maxDt = cplInterface.advance(maxDt);

      if (cplInterface.isCouplingOngoing()) {
        // readBlockVectorData without waveform
        cplInterface.readBlockVectorData(meshName, dataAID, 1, vertexIDs.data(), readDataA.data());
        // expected data value received
        BOOST_TEST(Vector3d(14.0, 14.0, 14.0) == readDataA);
        // reset value at read destination
        readDataA[0] = 0;
        readDataA[1] = 0;
        readDataA[2] = 0;
        // readBlockVectorData with waveform
        cplInterface.readBlockVectorData(meshName, dataAID, 1, vertexIDs.data(), 0.5, readDataA.data());
        // expected data value received
        BOOST_TEST(Vector3d(10.5, 10.5, 10.5) == readDataA);
        // reset value at read destination
        readDataA[0] = 0;
        readDataA[1] = 0;
        readDataA[2] = 0;

        // readVectorData without waveform
        cplInterface.readVectorData(meshName, dataAID, vertexIDs[0], readDataA.data());
        // expected data value received
        BOOST_TEST(Vector3d(14.0, 14.0, 14.0) == readDataA);
        // reset value at read destination
        readDataA[0] = 0;
        readDataA[1] = 0;
        readDataA[2] = 0;
        // readVectorData with waveform
        cplInterface.readVectorData(meshName, dataAID, vertexIDs[0], 0.5, readDataA.data());
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

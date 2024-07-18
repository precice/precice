#ifndef PRECICE_NO_MPI

#include "helpers.hpp"

#include "precice/precice.hpp"
#include "testing/Testing.hpp"

using namespace precice;

// Simple case of A <==> B
void multiCouplingTwoSolvers(const std::string configFile, const TestContext &context)
{
  Eigen::Vector2d coordOneA{0.0, 0.0};

  double valueA = 1.0;
  double valueB = 2.0;
  double valueC = 3.0;

  if (context.isNamed("SolverA")) {
    Participant cplInterface("SolverA", configFile, 0, 1);
    auto        meshName = "MeshA";
    int         vertexID = cplInterface.setMeshVertex(meshName, coordOneA);
    auto        dataABID = "DataAB";
    auto        dataBAID = "DataBA";

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();
    double valueRead;

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeData(meshName, dataABID, {&vertexID, 1}, {&valueA, 1});
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readData(meshName, dataBAID, {&vertexID, 1}, maxDt, {&valueRead, 1});
    }

    BOOST_TEST(valueRead == valueB);

    cplInterface.finalize();
  } else {
    Participant cplInterface("SolverB", configFile, 0, 1);
    auto        meshName = "MeshB";
    int         vertexID = cplInterface.setMeshVertex(meshName, coordOneA);
    auto        dataABID = "DataAB";
    auto        dataBAID = "DataBA";

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();
    double valueRead;

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeData(meshName, dataBAID, {&vertexID, 1}, {&valueB, 1});
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readData(meshName, dataABID, {&vertexID, 1}, maxDt, {&valueRead, 1});
    }

    BOOST_TEST(valueRead == valueA);

    cplInterface.finalize();
  }
}

// Simple case of A <==> B <==> C
void multiCouplingThreeSolvers(const std::string configFile, const TestContext &context)
{
  Eigen::Vector2d coordOneA{0.0, 0.0};

  double valueA = 1.0;
  double valueB = 2.0;
  double valueC = 3.0;

  if (context.isNamed("SolverA")) {
    Participant cplInterface("SolverA", configFile, 0, 1);
    auto        meshName = "MeshA";
    int         vertexID = cplInterface.setMeshVertex(meshName, coordOneA);
    auto        dataABID = "DataAB";
    auto        dataBAID = "DataBA";

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();
    double valueRead;

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeData(meshName, dataABID, {&vertexID, 1}, {&valueA, 1});
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readData(meshName, dataBAID, {&vertexID, 1}, maxDt, {&valueRead, 1});
    }

    BOOST_TEST(valueRead == valueB);

    cplInterface.finalize();
  } else if (context.isNamed("SolverB")) {
    Participant cplInterface("SolverB", configFile, 0, 1);
    auto        meshName1 = "MeshB1";
    auto        meshName2 = "MeshB2";
    int         vertexID1 = cplInterface.setMeshVertex(meshName1, coordOneA);
    int         vertexID2 = cplInterface.setMeshVertex(meshName2, coordOneA);
    auto        dataABID  = "DataAB";
    auto        dataBAID  = "DataBA";
    auto        dataCBID  = "DataCB";
    auto        dataBCID  = "DataBC";

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();
    double valueReadA, valueReadC;

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeData(meshName1, dataBAID, {&vertexID1, 1}, {&valueB, 1});
      cplInterface.writeData(meshName2, dataBCID, {&vertexID2, 1}, {&valueB, 1});
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readData(meshName1, dataABID, {&vertexID1, 1}, maxDt, {&valueReadA, 1});
      cplInterface.readData(meshName2, dataCBID, {&vertexID2, 1}, maxDt, {&valueReadC, 1});
    }

    BOOST_TEST(valueReadA == 1.0);
    BOOST_TEST(valueReadC == 3.0);

    cplInterface.finalize();

  } else {
    Participant cplInterface("SolverC", configFile, 0, 1);
    auto        meshName = "MeshC";
    int         vertexID = cplInterface.setMeshVertex(meshName, coordOneA);
    auto        dataCBID = "DataCB";
    auto        dataBCID = "DataBC";

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();
    double valueRead;

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {

      cplInterface.writeData(meshName, dataCBID, {&vertexID, 1}, {&valueC, 1});
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readData(meshName, dataBCID, {&vertexID, 1}, maxDt, {&valueRead, 1});
    }

    BOOST_TEST(valueRead == 2.0);

    cplInterface.finalize();
  }
}

void multiCouplingFourSolvers(const std::string configFile, const TestContext &context)
{
  Eigen::Vector2d coordOneA{0.0, 0.0};

  if (context.isNamed("SolverA")) {
    Participant cplInterface("SolverA", configFile, 0, 1);
    auto        meshName = "MeshA";
    int         vertexID = cplInterface.setMeshVertex(meshName, coordOneA);
    auto        dataABID = "DataAB";
    auto        dataBAID = "DataBA";

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();
    double valueRead;
    double valueWrite = 1.0;

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeData(meshName, dataABID, {&vertexID, 1}, {&valueWrite, 1});
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readData(meshName, dataBAID, {&vertexID, 1}, maxDt, {&valueRead, 1});
    }
    cplInterface.finalize();
  } else if (context.isNamed("SolverB")) {
    Participant cplInterface("SolverB", configFile, 0, 1);
    auto        meshName1 = "MeshB1";
    auto        meshName2 = "MeshB2";
    int         vertexID1 = cplInterface.setMeshVertex(meshName1, coordOneA);
    int         vertexID2 = cplInterface.setMeshVertex(meshName2, coordOneA);
    auto        dataABID  = "DataAB";
    auto        dataBAID  = "DataBA";
    auto        dataCBID  = "DataCB";
    auto        dataBCID  = "DataBC";

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();
    double valueReadA, valueReadC;
    double valueWriteA{1.0}, valueWriteC{1.0};

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeData(meshName1, dataBAID, {&vertexID1, 1}, {&valueWriteA, 1});
      cplInterface.writeData(meshName2, dataBCID, {&vertexID2, 1}, {&valueWriteC, 1});
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readData(meshName1, dataABID, {&vertexID1, 1}, maxDt, {&valueReadA, 1});
      cplInterface.readData(meshName2, dataCBID, {&vertexID2, 1}, maxDt, {&valueReadC, 1});
    }
    cplInterface.finalize();

  } else if (context.isNamed("SolverC")) {
    Participant cplInterface("SolverC", configFile, 0, 1);
    auto        meshName1 = "MeshC1";
    auto        meshName2 = "MeshC2";
    int         vertexID1 = cplInterface.setMeshVertex(meshName1, coordOneA);
    int         vertexID2 = cplInterface.setMeshVertex(meshName2, coordOneA);
    auto        dataBCID  = "DataBC";
    auto        dataCBID  = "DataCB";
    auto        dataCDID  = "DataCD";
    auto        dataDCID  = "DataDC";

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();
    double valueReadA, valueReadC;
    double valueWriteA{1.0}, valueWriteC{1.0};

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeData(meshName1, dataCBID, {&vertexID1, 1}, {&valueWriteA, 1});
      cplInterface.writeData(meshName2, dataCDID, {&vertexID2, 1}, {&valueWriteC, 1});
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readData(meshName1, dataBCID, {&vertexID1, 1}, maxDt, {&valueReadA, 1});
      cplInterface.readData(meshName2, dataDCID, {&vertexID2, 1}, maxDt, {&valueReadC, 1});
    }
    cplInterface.finalize();
  } else {
    Participant cplInterface("SolverD", configFile, 0, 1);
    auto        meshName = "MeshD";
    int         vertexID = cplInterface.setMeshVertex(meshName, coordOneA);
    auto        dataCDID = "DataCD";
    auto        dataDCID = "DataDC";

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();
    double valueRead;
    double valueWrite = 1.0;

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeData(meshName, dataDCID, {&vertexID, 1}, {&valueWrite, 1});
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readData(meshName, dataCDID, {&vertexID, 1}, maxDt, {&valueRead, 1});
    }
    cplInterface.finalize();
  }
}

#endif

#ifndef PRECICE_NO_MPI

#include "helpers.hpp"

#include "precice/SolverInterface.hpp"
#include "testing/Testing.hpp"

using namespace precice;

// Simple case of A <==> B <==> C
void multiCouplingThreeSolvers(const std::string configFile, const TestContext &context)
{
  Eigen::Vector2d coordOneA{0.0, 0.0};

  double valueA = 1.0;
  double valueB = 2.0;
  double valueC = 3.0;

  if (context.isNamed("SolverA")) {
    SolverInterface cplInterface("SolverA", configFile, 0, 1);
    const int       meshID   = cplInterface.getMeshID("MeshA");
    int             vertexID = cplInterface.setMeshVertex(meshID, coordOneA.data());
    int             dataABID = cplInterface.getDataID("DataAB", meshID);
    int             dataBAID = cplInterface.getDataID("DataBA", meshID);

    double maxDt = cplInterface.initialize();
    double valueRead;

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(dataABID, vertexID, valueA);
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readScalarData(dataBAID, vertexID, valueRead);
    }

    BOOST_TEST(valueRead == valueB);

    cplInterface.finalize();
  } else if (context.isNamed("SolverB")) {
    SolverInterface cplInterface("SolverB", configFile, 0, 1);
    const int       meshID1   = cplInterface.getMeshID("MeshB1");
    const int       meshID2   = cplInterface.getMeshID("MeshB2");
    int             vertexID1 = cplInterface.setMeshVertex(meshID1, coordOneA.data());
    int             vertexID2 = cplInterface.setMeshVertex(meshID2, coordOneA.data());
    int             dataABID  = cplInterface.getDataID("DataAB", meshID1);
    int             dataBAID  = cplInterface.getDataID("DataBA", meshID1);
    int             dataCBID  = cplInterface.getDataID("DataCB", meshID2);
    int             dataBCID  = cplInterface.getDataID("DataBC", meshID2);

    double maxDt = cplInterface.initialize();
    double valueReadA, valueReadC;

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(dataBAID, vertexID1, valueB);
      cplInterface.writeScalarData(dataBCID, vertexID2, valueB);
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readScalarData(dataABID, vertexID1, valueReadA);
      cplInterface.readScalarData(dataCBID, vertexID2, valueReadC);
    }

    BOOST_TEST(valueReadA == 1.0);
    BOOST_TEST(valueReadC == 3.0);

    cplInterface.finalize();

  } else {
    SolverInterface cplInterface("SolverC", configFile, 0, 1);
    const int       meshID   = cplInterface.getMeshID("MeshC");
    int             vertexID = cplInterface.setMeshVertex(meshID, coordOneA.data());
    int             dataCBID = cplInterface.getDataID("DataCB", meshID);
    int             dataBCID = cplInterface.getDataID("DataBC", meshID);

    double maxDt = cplInterface.initialize();
    double valueRead;

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {

      cplInterface.writeScalarData(dataCBID, vertexID, valueC);
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readScalarData(dataBCID, vertexID, valueRead);
    }

    BOOST_TEST(valueRead == 2.0);

    cplInterface.finalize();
  }
}

void multiCouplingFourSolvers(const std::string configFile, const TestContext &context)
{
  Eigen::Vector2d coordOneA{0.0, 0.0};

  if (context.isNamed("SolverA")) {
    SolverInterface cplInterface("SolverA", configFile, 0, 1);
    const int       meshID   = cplInterface.getMeshID("MeshA");
    int             vertexID = cplInterface.setMeshVertex(meshID, coordOneA.data());
    int             dataABID = cplInterface.getDataID("DataAB", meshID);
    int             dataBAID = cplInterface.getDataID("DataBA", meshID);

    double maxDt = cplInterface.initialize();
    double valueRead;
    double valueWrite = 1.0;

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(dataABID, vertexID, valueWrite);
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readScalarData(dataBAID, vertexID, valueRead);
    }
    cplInterface.finalize();
  } else if (context.isNamed("SolverB")) {
    SolverInterface cplInterface("SolverB", configFile, 0, 1);
    const int       meshID1   = cplInterface.getMeshID("MeshB1");
    const int       meshID2   = cplInterface.getMeshID("MeshB2");
    int             vertexID1 = cplInterface.setMeshVertex(meshID1, coordOneA.data());
    int             vertexID2 = cplInterface.setMeshVertex(meshID2, coordOneA.data());
    int             dataABID  = cplInterface.getDataID("DataAB", meshID1);
    int             dataBAID  = cplInterface.getDataID("DataBA", meshID1);
    int             dataCBID  = cplInterface.getDataID("DataCB", meshID2);
    int             dataBCID  = cplInterface.getDataID("DataBC", meshID2);

    double maxDt = cplInterface.initialize();
    double valueReadA, valueReadC;
    double valueWriteA{1.0}, valueWriteC{1.0};

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(dataBAID, vertexID1, valueWriteA);
      cplInterface.writeScalarData(dataBCID, vertexID2, valueWriteC);
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readScalarData(dataABID, vertexID1, valueReadA);
      cplInterface.readScalarData(dataCBID, vertexID2, valueReadC);
    }
    cplInterface.finalize();

  } else if (context.isNamed("SolverC")) {
    SolverInterface cplInterface("SolverC", configFile, 0, 1);
    const int       meshID1   = cplInterface.getMeshID("MeshC1");
    const int       meshID2   = cplInterface.getMeshID("MeshC2");
    int             vertexID1 = cplInterface.setMeshVertex(meshID1, coordOneA.data());
    int             vertexID2 = cplInterface.setMeshVertex(meshID2, coordOneA.data());
    int             dataBCID  = cplInterface.getDataID("DataBC", meshID1);
    int             dataCBID  = cplInterface.getDataID("DataCB", meshID1);
    int             dataCDID  = cplInterface.getDataID("DataCD", meshID2);
    int             dataDCID  = cplInterface.getDataID("DataDC", meshID2);

    double maxDt = cplInterface.initialize();
    double valueReadA, valueReadC;
    double valueWriteA{1.0}, valueWriteC{1.0};

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(dataCBID, vertexID1, valueWriteA);
      cplInterface.writeScalarData(dataCDID, vertexID2, valueWriteC);
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readScalarData(dataBCID, vertexID1, valueReadA);
      cplInterface.readScalarData(dataDCID, vertexID2, valueReadC);
    }
    cplInterface.finalize();
  } else {
    SolverInterface cplInterface("SolverD", configFile, 0, 1);
    const int       meshID   = cplInterface.getMeshID("MeshD");
    int             vertexID = cplInterface.setMeshVertex(meshID, coordOneA.data());
    int             dataCDID = cplInterface.getDataID("DataCD", meshID);
    int             dataDCID = cplInterface.getDataID("DataDC", meshID);

    double maxDt = cplInterface.initialize();
    double valueRead;
    double valueWrite = 1.0;

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(dataDCID, vertexID, valueWrite);
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readScalarData(dataCDID, vertexID, valueRead);
    }
    cplInterface.finalize();
  }
}

#endif

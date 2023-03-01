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
    auto            meshID   = "MeshA";
    int             vertexID = cplInterface.setMeshVertex(meshID, coordOneA.data());
    auto            dataABID = "DataAB"; //  meshID
    auto            dataBAID = "DataBA"; //  meshID

    double maxDt = cplInterface.initialize();
    double valueRead;

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(meshID, dataABID, vertexID, valueA);
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readScalarData(meshID, dataBAID, vertexID, valueRead);
    }

    BOOST_TEST(valueRead == valueB);

    cplInterface.finalize();
  } else if (context.isNamed("SolverB")) {
    SolverInterface cplInterface("SolverB", configFile, 0, 1);
    auto            meshID1   = "MeshB1";
    auto            meshID2   = "MeshB2";
    int             vertexID1 = cplInterface.setMeshVertex(meshID1, coordOneA.data());
    int             vertexID2 = cplInterface.setMeshVertex(meshID2, coordOneA.data());
    auto            dataABID  = "DataAB"; //  meshID1
    auto            dataBAID  = "DataBA"; //  meshID1
    auto            dataCBID  = "DataCB"; //  meshID2
    auto            dataBCID  = "DataBC"; //  meshID2

    double maxDt = cplInterface.initialize();
    double valueReadA, valueReadC;

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(meshID, dataBAID, vertexID1, valueB);
      cplInterface.writeScalarData(meshID, dataBCID, vertexID2, valueB);
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readScalarData(meshID, dataABID, vertexID1, valueReadA);
      cplInterface.readScalarData(meshID, dataCBID, vertexID2, valueReadC);
    }

    BOOST_TEST(valueReadA == 1.0);
    BOOST_TEST(valueReadC == 3.0);

    cplInterface.finalize();

  } else {
    SolverInterface cplInterface("SolverC", configFile, 0, 1);
    auto            meshID   = "MeshC";
    int             vertexID = cplInterface.setMeshVertex(meshID, coordOneA.data());
    auto            dataCBID = "DataCB"; //  meshID
    auto            dataBCID = "DataBC"; //  meshID

    double maxDt = cplInterface.initialize();
    double valueRead;

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {

      cplInterface.writeScalarData(meshID, dataCBID, vertexID, valueC);
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readScalarData(meshID, dataBCID, vertexID, valueRead);
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
    auto            meshID   = "MeshA";
    int             vertexID = cplInterface.setMeshVertex(meshID, coordOneA.data());
    auto            dataABID = "DataAB"; //  meshID
    auto            dataBAID = "DataBA"; //  meshID

    double maxDt = cplInterface.initialize();
    double valueRead;
    double valueWrite = 1.0;

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(meshID, dataABID, vertexID, valueWrite);
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readScalarData(meshID, dataBAID, vertexID, valueRead);
    }
    cplInterface.finalize();
  } else if (context.isNamed("SolverB")) {
    SolverInterface cplInterface("SolverB", configFile, 0, 1);
    auto            meshID1   = "MeshB1";
    auto            meshID2   = "MeshB2";
    int             vertexID1 = cplInterface.setMeshVertex(meshID1, coordOneA.data());
    int             vertexID2 = cplInterface.setMeshVertex(meshID2, coordOneA.data());
    auto            dataABID  = "DataAB"; //  meshID1
    auto            dataBAID  = "DataBA"; //  meshID1
    auto            dataCBID  = "DataCB"; //  meshID2
    auto            dataBCID  = "DataBC"; //  meshID2

    double maxDt = cplInterface.initialize();
    double valueReadA, valueReadC;
    double valueWriteA{1.0}, valueWriteC{1.0};

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(meshID, dataBAID, vertexID1, valueWriteA);
      cplInterface.writeScalarData(meshID, dataBCID, vertexID2, valueWriteC);
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readScalarData(meshID, dataABID, vertexID1, valueReadA);
      cplInterface.readScalarData(meshID, dataCBID, vertexID2, valueReadC);
    }
    cplInterface.finalize();

  } else if (context.isNamed("SolverC")) {
    SolverInterface cplInterface("SolverC", configFile, 0, 1);
    auto            meshID1   = "MeshC1";
    auto            meshID2   = "MeshC2";
    int             vertexID1 = cplInterface.setMeshVertex(meshID1, coordOneA.data());
    int             vertexID2 = cplInterface.setMeshVertex(meshID2, coordOneA.data());
    auto            dataBCID  = "DataBC"; //  meshID1
    auto            dataCBID  = "DataCB"; //  meshID1
    auto            dataCDID  = "DataCD"; //  meshID2
    auto            dataDCID  = "DataDC"; //  meshID2

    double maxDt = cplInterface.initialize();
    double valueReadA, valueReadC;
    double valueWriteA{1.0}, valueWriteC{1.0};

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(meshID, dataCBID, vertexID1, valueWriteA);
      cplInterface.writeScalarData(meshID, dataCDID, vertexID2, valueWriteC);
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readScalarData(meshID, dataBCID, vertexID1, valueReadA);
      cplInterface.readScalarData(meshID, dataDCID, vertexID2, valueReadC);
    }
    cplInterface.finalize();
  } else {
    SolverInterface cplInterface("SolverD", configFile, 0, 1);
    auto            meshID   = "MeshD";
    int             vertexID = cplInterface.setMeshVertex(meshID, coordOneA.data());
    auto            dataCDID = "DataCD"; //  meshID
    auto            dataDCID = "DataDC"; //  meshID

    double maxDt = cplInterface.initialize();
    double valueRead;
    double valueWrite = 1.0;

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(meshID, dataDCID, vertexID, valueWrite);
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readScalarData(meshID, dataCDID, vertexID, valueRead);
    }
    cplInterface.finalize();
  }
}

#endif

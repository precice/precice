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
    auto            meshName = "MeshA";
    int             vertexID = cplInterface.setMeshVertex(meshName, coordOneA.data());
    auto            dataABID = "DataAB"; //  meshName
    auto            dataBAID = "DataBA"; //  meshName

    double maxDt = cplInterface.initialize();
    double valueRead;

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(meshName, dataABID, vertexID, valueA);
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readScalarData(meshName, dataBAID, vertexID, valueRead);
    }

    BOOST_TEST(valueRead == valueB);

    cplInterface.finalize();
  } else if (context.isNamed("SolverB")) {
    SolverInterface cplInterface("SolverB", configFile, 0, 1);
    auto            meshName1 = "MeshB1";
    auto            meshName2 = "MeshB2";
    int             vertexID1 = cplInterface.setMeshVertex(meshName1, coordOneA.data());
    int             vertexID2 = cplInterface.setMeshVertex(meshName2, coordOneA.data());
    auto            dataABID  = "DataAB"; //  meshName1
    auto            dataBAID  = "DataBA"; //  meshName1
    auto            dataCBID  = "DataCB"; //  meshName2
    auto            dataBCID  = "DataBC"; //  meshName2

    double maxDt = cplInterface.initialize();
    double valueReadA, valueReadC;

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(meshName1, dataBAID, vertexID1, valueB);
      cplInterface.writeScalarData(meshName2, dataBCID, vertexID2, valueB);
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readScalarData(meshName1, dataABID, vertexID1, valueReadA);
      cplInterface.readScalarData(meshName2, dataCBID, vertexID2, valueReadC);
    }

    BOOST_TEST(valueReadA == 1.0);
    BOOST_TEST(valueReadC == 3.0);

    cplInterface.finalize();

  } else {
    SolverInterface cplInterface("SolverC", configFile, 0, 1);
    auto            meshName = "MeshC";
    int             vertexID = cplInterface.setMeshVertex(meshName, coordOneA.data());
    auto            dataCBID = "DataCB"; //  meshName
    auto            dataBCID = "DataBC"; //  meshName

    double maxDt = cplInterface.initialize();
    double valueRead;

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {

      cplInterface.writeScalarData(meshName, dataCBID, vertexID, valueC);
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readScalarData(meshName, dataBCID, vertexID, valueRead);
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
    auto            meshName = "MeshA";
    int             vertexID = cplInterface.setMeshVertex(meshName, coordOneA.data());
    auto            dataABID = "DataAB"; //  meshName
    auto            dataBAID = "DataBA"; //  meshName

    double maxDt = cplInterface.initialize();
    double valueRead;
    double valueWrite = 1.0;

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(meshName, dataABID, vertexID, valueWrite);
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readScalarData(meshName, dataBAID, vertexID, valueRead);
    }
    cplInterface.finalize();
  } else if (context.isNamed("SolverB")) {
    SolverInterface cplInterface("SolverB", configFile, 0, 1);
    auto            meshName1 = "MeshB1";
    auto            meshName2 = "MeshB2";
    int             vertexID1 = cplInterface.setMeshVertex(meshName1, coordOneA.data());
    int             vertexID2 = cplInterface.setMeshVertex(meshName2, coordOneA.data());
    auto            dataABID  = "DataAB"; //  meshName1
    auto            dataBAID  = "DataBA"; //  meshName1
    auto            dataCBID  = "DataCB"; //  meshName2
    auto            dataBCID  = "DataBC"; //  meshName2

    double maxDt = cplInterface.initialize();
    double valueReadA, valueReadC;
    double valueWriteA{1.0}, valueWriteC{1.0};

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(meshName1, dataBAID, vertexID1, valueWriteA);
      cplInterface.writeScalarData(meshName2, dataBCID, vertexID2, valueWriteC);
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readScalarData(meshName1, dataABID, vertexID1, valueReadA);
      cplInterface.readScalarData(meshName2, dataCBID, vertexID2, valueReadC);
    }
    cplInterface.finalize();

  } else if (context.isNamed("SolverC")) {
    SolverInterface cplInterface("SolverC", configFile, 0, 1);
    auto            meshName1 = "MeshC1";
    auto            meshName2 = "MeshC2";
    int             vertexID1 = cplInterface.setMeshVertex(meshName1, coordOneA.data());
    int             vertexID2 = cplInterface.setMeshVertex(meshName2, coordOneA.data());
    auto            dataBCID  = "DataBC"; //  meshName1
    auto            dataCBID  = "DataCB"; //  meshName1
    auto            dataCDID  = "DataCD"; //  meshName2
    auto            dataDCID  = "DataDC"; //  meshName2

    double maxDt = cplInterface.initialize();
    double valueReadA, valueReadC;
    double valueWriteA{1.0}, valueWriteC{1.0};

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(meshName1, dataCBID, vertexID1, valueWriteA);
      cplInterface.writeScalarData(meshName2, dataCDID, vertexID2, valueWriteC);
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readScalarData(meshName1, dataBCID, vertexID1, valueReadA);
      cplInterface.readScalarData(meshName2, dataDCID, vertexID2, valueReadC);
    }
    cplInterface.finalize();
  } else {
    SolverInterface cplInterface("SolverD", configFile, 0, 1);
    auto            meshName = "MeshD";
    int             vertexID = cplInterface.setMeshVertex(meshName, coordOneA.data());
    auto            dataCDID = "DataCD"; //  meshName
    auto            dataDCID = "DataDC"; //  meshName

    double maxDt = cplInterface.initialize();
    double valueRead;
    double valueWrite = 1.0;

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(meshName, dataDCID, vertexID, valueWrite);
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);

      if (cplInterface.requiresReadingCheckpoint()) {
      }
      cplInterface.readScalarData(meshName, dataCDID, vertexID, valueRead);
    }
    cplInterface.finalize();
  }
}

#endif

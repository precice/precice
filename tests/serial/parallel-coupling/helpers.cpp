#ifndef PRECICE_NO_MPI

#include "helpers.hpp"

#include "precice/precice.hpp"
#include "testing/Testing.hpp"

using namespace precice;

void parallelCoupling(const std::string configFile, const TestContext &context)
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

#endif

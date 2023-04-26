#ifndef PRECICE_NO_MPI

#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include <mpi.h>
#include <string>
#include <vector>

#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "math/constants.hpp"
#include "math/geometry.hpp"
#include "mesh/Mesh.hpp"
#include "precice/SolverInterface.hpp"
#include "precice/config/Configuration.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"
#include "precice/types.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/IntraComm.hpp"
#include "utils/Parallel.hpp"

using namespace precice;
using testing::TestContext;

struct ParallelTestFixture : testing::WhiteboxAccessor {

  std::string _pathToTests;

  ParallelTestFixture()
  {
    _pathToTests = testing::getPathToSources() + "/precice/tests/";
  }
};

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_FIXTURE_TEST_SUITE(Parallel, ParallelTestFixture)

// Simple case of A <==> B <==> C
void multiCouplingThreeSolversParallelControl(const std::string configFile, const TestContext &context)
{
  Eigen::Vector2d coordOneA{0.0, 0.0};
  Eigen::Vector2d coordOneB{1.0, 0.0};

  double valueA1 = 1.0;
  double valueA2 = 1.5;
  double valueB1 = 2.0;
  double valueB2 = 2.5;
  double valueC1 = 3.0;
  double valueC2 = 3.5;

  if (context.isNamed("SolverA")) {
    SolverInterface cplInterface("SolverA", configFile, context.rank, context.size);
    auto            meshName = "MeshA";
    auto            dataABID = "DataAB";
    auto            dataBAID = "DataBA";

    if (context.isPrimary()) {
      int vertex1 = cplInterface.setMeshVertex(meshName, coordOneA.data());

      cplInterface.initialize();
      double maxDt = cplInterface.getMaxTimeStepSize();
      double valueRead;

      BOOST_TEST(cplInterface.isCouplingOngoing());
      while (cplInterface.isCouplingOngoing()) {
        cplInterface.writeScalarData(meshName, dataABID, vertex1, valueA1);
        if (cplInterface.requiresWritingCheckpoint()) {
        }

        cplInterface.advance(maxDt);

        if (cplInterface.requiresReadingCheckpoint()) {
        }

        maxDt = cplInterface.getMaxTimeStepSize();
        cplInterface.readScalarData(meshName, dataBAID, vertex1, maxDt, valueRead);
      }

      BOOST_TEST(valueRead == valueB1);

      cplInterface.finalize();

    } else {
      int vertex2 = cplInterface.setMeshVertex(meshName, coordOneB.data());

      cplInterface.initialize();
      double maxDt = cplInterface.getMaxTimeStepSize();
      double valueRead;

      BOOST_TEST(cplInterface.isCouplingOngoing());
      while (cplInterface.isCouplingOngoing()) {
        cplInterface.writeScalarData(meshName, dataABID, vertex2, valueA2);
        if (cplInterface.requiresWritingCheckpoint()) {
        }

        cplInterface.advance(maxDt);

        if (cplInterface.requiresReadingCheckpoint()) {
        }

        maxDt = cplInterface.getMaxTimeStepSize();
        cplInterface.readScalarData(meshName, dataBAID, vertex2, maxDt, valueRead);
      }

      BOOST_TEST(valueRead == valueB2);

      cplInterface.finalize();
    }

  } else if (context.isNamed("SolverB")) {
    SolverInterface cplInterface("SolverB", configFile, 0, 1);
    auto            meshName1 = "MeshB1";
    auto            meshName2 = "MeshB2";
    int             vertex1   = cplInterface.setMeshVertex(meshName1, coordOneA.data());
    int             vertex2   = cplInterface.setMeshVertex(meshName1, coordOneB.data());
    int             vertex3   = cplInterface.setMeshVertex(meshName2, coordOneA.data());
    int             vertex4   = cplInterface.setMeshVertex(meshName2, coordOneB.data());

    auto dataABID = "DataAB"; // meshName1
    auto dataBAID = "DataBA"; // meshName1
    auto dataCBID = "DataCB"; // meshName2;
    auto dataBCID = "DataBC"; // meshName2;

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();
    double valueReadA1, valueReadA2, valueReadC1, valueReadC2;

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeScalarData(meshName1, dataBAID, vertex1, valueB1);
      cplInterface.writeScalarData(meshName1, dataBAID, vertex2, valueB2);
      cplInterface.writeScalarData(meshName2, dataBCID, vertex3, valueB1);
      cplInterface.writeScalarData(meshName2, dataBCID, vertex4, valueB2);
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);

      if (cplInterface.requiresReadingCheckpoint()) {
      }

      maxDt = cplInterface.getMaxTimeStepSize();
      cplInterface.readScalarData(meshName1, dataABID, vertex1, maxDt, valueReadA1);
      cplInterface.readScalarData(meshName1, dataABID, vertex2, maxDt, valueReadA2);
      cplInterface.readScalarData(meshName2, dataCBID, vertex1, maxDt, valueReadC1);
      cplInterface.readScalarData(meshName2, dataCBID, vertex2, maxDt, valueReadC2);
    }

    BOOST_TEST(valueReadA1 == valueA1);
    BOOST_TEST(valueReadA2 == valueA2);
    BOOST_TEST(valueReadC1 == valueC1);
    BOOST_TEST(valueReadC2 == valueC2);

    cplInterface.finalize();

  } else {
    SolverInterface cplInterface("SolverC", configFile, 0, 1);
    auto            meshName = "MeshC";
    int             vertex1  = cplInterface.setMeshVertex(meshName, coordOneA.data());
    int             vertex2  = cplInterface.setMeshVertex(meshName, coordOneB.data());
    auto            dataCBID = "DataCB";
    auto            dataBCID = "DataBC";

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();
    double valueRead1, valueRead2;

    BOOST_TEST(cplInterface.isCouplingOngoing());
    while (cplInterface.isCouplingOngoing()) {

      cplInterface.writeScalarData(meshName, dataCBID, vertex1, valueC1);
      cplInterface.writeScalarData(meshName, dataCBID, vertex2, valueC2);
      if (cplInterface.requiresWritingCheckpoint()) {
      }

      cplInterface.advance(maxDt);

      if (cplInterface.requiresReadingCheckpoint()) {
      }

      maxDt = cplInterface.getMaxTimeStepSize();
      cplInterface.readScalarData(meshName, dataBCID, vertex1, maxDt, valueRead1);
      cplInterface.readScalarData(meshName, dataBCID, vertex2, maxDt, valueRead2);
    }

    BOOST_TEST(valueRead1 == valueB1);
    BOOST_TEST(valueRead2 == valueB2);

    cplInterface.finalize();
  }
}

// BOOST_AUTO_TEST_CASE(MultiCouplingThreeSolversParallelCentral1)
// {
//   PRECICE_TEST("SolverA"_on(2_ranks), "SolverB"_on(1_rank), "SolverC"_on(1_rank));
//   const std::string configFile = _pathToTests + "multi-coupling-three-solver-1.xml";
//   multiCouplingThreeSolversParallelControl(configFile, context);
// }

// BOOST_AUTO_TEST_CASE(MultiCouplingThreeSolversParallelCentral2)
// {
//   PRECICE_TEST("SolverA"_on(2_ranks), "SolverB"_on(1_rank), "SolverC"_on(1_rank));
//   const std::string configFile = _pathToTests + "multi-coupling-three-solver-2.xml";
//   multiCouplingThreeSolversParallelControl(configFile, context);
// }

// BOOST_AUTO_TEST_CASE(MultiCouplingThreeSolversParallelCentral3)
// {
//   PRECICE_TEST("SolverA"_on(2_ranks), "SolverB"_on(1_rank), "SolverC"_on(1_rank));
//   const std::string configFile = _pathToTests + "multi-coupling-three-solver-3.xml";
//   multiCouplingThreeSolversParallelControl(configFile, context);
// }

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
#endif // PRECICE_NO_MPI

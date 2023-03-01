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

//std::string pathToTests = testing::getPathToSources() + "/tests/serial/mapping-nearest-neighbor-gradient/";

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MappingNearestNeighborGradient)

// Bidirectional test : Read: Vector & NN - Write: Scalar & NNG (Serial coupling)
BOOST_AUTO_TEST_CASE(GradientTestBidirectionalWriteScalar)
{

  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  using Eigen::Vector2d;
  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, context.config(), 0, 1);
  if (context.isNamed("SolverOne")) {
    auto     meshOneID = "MeshOne";
    Vector3d vec1      = Vector3d::Constant(0.1);
    cplInterface.setMeshVertex(meshOneID, vec1.data());
    auto dataAID = "DataOne"; //  meshOneID
    auto dataBID = "DataTwo"; //  meshOneID

    double valueDataB = 0.0;
    double maxDt      = cplInterface.initialize();
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
    auto     meshTwoID = "MeshTwo";
    Vector3d vec2      = Vector3d::Constant(0.0);
    cplInterface.setMeshVertex(meshTwoID, vec2.data());

    auto dataAID = "DataOne"; //  meshTwoID
    auto dataBID = "DataTwo"; //  meshTwoID
    BOOST_REQUIRE(cplInterface.requiresInitialData());

    double   valueDataB = 1.0;
    Vector3d valueGradDataB(1.0, 1.0, 1.0);
    cplInterface.writeScalarData(dataBID, 0, valueDataB);
    cplInterface.writeScalarGradientData(dataBID, 0, valueGradDataB.data());

    //tell preCICE that data has been written and call initialize
    double maxDt = cplInterface.initialize();

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

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI

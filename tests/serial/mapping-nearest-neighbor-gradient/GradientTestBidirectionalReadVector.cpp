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

// Bidirectional test : Read: Vector & NNG - Write: Vector & NN (Serial coupling)
BOOST_AUTO_TEST_CASE(GradientTestBidirectionalReadVector)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))

  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, context.config(), 0, 1);
  if (context.isNamed("SolverOne")) {
    auto     meshName = "MeshOne";
    Vector3d posOne   = Vector3d::Constant(0.0);
    cplInterface.setMeshVertex(meshName, posOne.data());
    auto dataAID = "DataOne";
    auto dataBID = "DataTwo";

    Vector3d valueDataB;

    double maxDt = cplInterface.initialize();
    cplInterface.readVectorData(meshName, dataBID, 0, valueDataB.data());
    Vector3d expected(2.0, 3.0, 4.0);
    BOOST_TEST(valueDataB == expected);

    while (cplInterface.isCouplingOngoing()) {
      Vector3d                    valueDataA(1.0, 1.0, 1.0);
      Eigen::Matrix<double, 3, 3> gradient;
      gradient << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;
      cplInterface.writeVectorData(meshName, dataAID, 0, valueDataA.data());
      cplInterface.writeVectorGradientData(meshName, dataAID, 0, gradient.data());

      maxDt = cplInterface.advance(maxDt);

      cplInterface.readVectorData(meshName, dataBID, 0, valueDataB.data());
      expected << 2.5, 3.5, 4.5;
      BOOST_TEST(valueDataB == expected);
    }
    cplInterface.finalize();

  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    auto     meshName = "MeshTwo";
    Vector3d pos      = Vector3d::Constant(1.0);
    cplInterface.setMeshVertex(meshName, pos.data());

    auto dataAID = "DataOne";
    auto dataBID = "DataTwo";

    BOOST_REQUIRE(cplInterface.requiresInitialData());
    Vector3d valueDataB(2.0, 3.0, 4.0);
    cplInterface.writeVectorData(meshName, dataBID, 0, valueDataB.data());

    //tell preCICE that data has been written and call initialize
    double maxDt = cplInterface.initialize();

    Vector3d valueDataA;
    cplInterface.readVectorData(meshName, dataAID, 0, valueDataA.data());
    Vector3d expected(4.0, 4.0, 4.0);
    BOOST_TEST(valueDataA == expected);

    while (cplInterface.isCouplingOngoing()) {

      valueDataB << 2.5, 3.5, 4.5;
      cplInterface.writeVectorData(meshName, dataBID, 0, valueDataB.data());

      maxDt = cplInterface.advance(maxDt);
      cplInterface.readVectorData(meshName, dataAID, 0, valueDataA.data());
      BOOST_TEST(valueDataA == expected);
    }
    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI

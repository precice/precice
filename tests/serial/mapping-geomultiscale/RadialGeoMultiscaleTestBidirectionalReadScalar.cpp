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
#include "precice/Participant.hpp"
#include "precice/impl/MeshContext.hpp"
#include "precice/impl/ParticipantImpl.hpp"
#include "precice/impl/ParticipantState.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "precice/types.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using precice::testing::TestContext;

//std::string pathToTests = testing::getPathToSources() + "/tests/serial/mapping-nearest-neighbor-gradient/";

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(RadialGeoMultiscaleMapping)

// Bidirectional test : Read: Vector & NNG - Write: Vector & NN (Serial coupling)
BOOST_AUTO_TEST_CASE(RadialGeoTestBidirectionalReadVector)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  using Eigen::Vector3d;

  Participant cplInterface(context.name, context.config(), 0, 1);
  if (context.isNamed("SolverOne")) {
    auto     meshName = "MeshOne";
    Vector3d vec1     = Vector3d::Constant(0.0);
    auto     vid      = cplInterface.setMeshVertex(meshName, vec1);
    auto     dataAID  = "DataOne";
    auto     dataBID  = "DataTwo";

    double valueDataB = 1.0;
    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();
    cplInterface.readData(meshName, dataBID, {&vid, 1}, maxDt, {&valueDataB, 1});
    BOOST_TEST(valueDataB == 1.0);

    while (cplInterface.isCouplingOngoing()) {

      double data[] = {1.0};
      cplInterface.writeData(meshName, dataAID, {&vid, 1}, data);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      cplInterface.readData(meshName, dataBID, {&vid, 1}, maxDt, {&valueDataB, 1});
      BOOST_TEST(valueDataB == 1.0);
    }
    cplInterface.finalize();

  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    auto     meshName = "MeshTwo";
    Vector3d vec2     = Vector3d::Constant(1.0);
    auto     vid      = cplInterface.setMeshVertex(meshName, vec2);

    auto dataAID = "DataOne";
    auto dataBID = "DataTwo";
    BOOST_REQUIRE(cplInterface.requiresInitialData());
    BOOST_TEST(cplInterface.requiresGradientDataFor(meshName, dataBID) == false);

    double valueDataB = 1.0;
    cplInterface.writeData(meshName, dataBID, {&vid, 1}, {&valueDataB, 1});

    //tell preCICE that data has been written and call initialize
    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    while (cplInterface.isCouplingOngoing()) {
      double data[] = {1.0};
      cplInterface.writeData(meshName, dataBID, {&vid, 1}, data);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      double valueDataA;
      cplInterface.readData(meshName, dataAID, {&vid, 1}, maxDt, {&valueDataA, 1});
      BOOST_TEST(valueDataA == 1.0);
    }
    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI

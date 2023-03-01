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

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MappingNearestNeighborGradient)

// Unidirectional Nearest Neighbor Gradient Read Mapping
// Also to test writeBlockScalarGradientData method
BOOST_AUTO_TEST_CASE(GradientTestUnidirectionalReadScalar)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank))
  using Eigen::Vector3d;

  SolverInterface cplInterface(context.name, context.config(), 0, 1);
  if (context.isNamed("A")) {

    auto meshOneID = "MeshA";
    auto dataID    = "DataA"; //  meshOneID

    Vector3d posOne = Vector3d::Constant(0.0);
    Vector3d posTwo = Vector3d::Constant(1.0);
    cplInterface.setMeshVertex(meshOneID, posOne.data());
    cplInterface.setMeshVertex(meshOneID, posTwo.data());

    // Initialize, thus sending the mesh.
    double maxDt = cplInterface.initialize();
    BOOST_TEST(cplInterface.isCouplingOngoing(), "Sending participant should have to advance once!");

    double values[2]  = {1.0, 2.0};
    int    indices[2] = {0, 1};
    cplInterface.writeBlockScalarData(meshID, dataID, 2, indices, values);

    BOOST_TEST(cplInterface.requiresGradientDataFor(meshID, dataID) == true);

    if (cplInterface.requiresGradientDataFor(meshID, dataID)) {
      double gradientValues[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
      cplInterface.writeBlockScalarGradientData(meshID, dataID, 2, indices, gradientValues);
    }

    // Participant must make move after writing
    maxDt = cplInterface.advance(maxDt);

    BOOST_TEST(!cplInterface.isCouplingOngoing(), "Sending participant should have to advance once!");
    cplInterface.finalize();

  } else {
    BOOST_TEST(context.isNamed("B"));
    auto meshTwoID = "MeshB";
    auto dataID    = "DataA"; //  meshTwoID

    Vector3d posOne = Vector3d::Constant(0.1);
    Vector3d posTwo = Vector3d::Constant(1.1);
    cplInterface.setMeshVertex(meshTwoID, posOne.data());
    cplInterface.setMeshVertex(meshTwoID, posTwo.data());

    double maxDt = cplInterface.initialize();
    BOOST_TEST(cplInterface.isCouplingOngoing(), "Receiving participant should have to advance once!");

    double valueData[2];
    int    indices[2] = {0, 1};
    cplInterface.readBlockScalarData(meshID, dataID, 2, indices, valueData);
    double expected[2] = {1.6, 3.5};
    BOOST_TEST(valueData == expected);

    cplInterface.advance(maxDt);
    BOOST_TEST(!cplInterface.isCouplingOngoing(), "Receiving participant should have to advance once!");

    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI

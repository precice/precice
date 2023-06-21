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
#include "precice/impl/MeshContext.hpp"
#include "precice/impl/ParticipantImpl.hpp"
#include "precice/impl/ParticipantState.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "precice/precice.hpp"
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

  Participant cplInterface(context.name, context.config(), 0, 1);
  if (context.isNamed("A")) {

    auto meshName = "MeshA";
    auto dataName = "DataA";

    Vector3d posOne = Vector3d::Constant(0.0);
    Vector3d posTwo = Vector3d::Constant(1.0);
    cplInterface.setMeshVertex(meshName, posOne);
    cplInterface.setMeshVertex(meshName, posTwo);

    // Initialize, thus sending the mesh.
    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();
    BOOST_TEST(cplInterface.isCouplingOngoing(), "Sending participant should have to advance once!");

    double values[2]  = {1.0, 2.0};
    int    indices[2] = {0, 1};
    cplInterface.writeData(meshName, dataName, indices, values);

    BOOST_TEST(cplInterface.requiresGradientDataFor(meshName, dataName) == true);

    if (cplInterface.requiresGradientDataFor(meshName, dataName)) {
      double gradients[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
      cplInterface.writeGradientData(meshName, dataName, indices, gradients);
    }

    // Participant must make move after writing
    cplInterface.advance(maxDt);
    maxDt = cplInterface.getMaxTimeStepSize();

    BOOST_TEST(!cplInterface.isCouplingOngoing(), "Sending participant should have to advance once!");
    cplInterface.finalize();

  } else {
    BOOST_TEST(context.isNamed("B"));
    auto meshName = "MeshB";
    auto dataName = "DataA";

    Vector3d posOne = Vector3d::Constant(0.1);
    Vector3d posTwo = Vector3d::Constant(1.1);
    cplInterface.setMeshVertex(meshName, posOne);
    cplInterface.setMeshVertex(meshName, posTwo);

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();
    BOOST_TEST(cplInterface.isCouplingOngoing(), "Receiving participant should have to advance once!");

    double valueData[2];
    int    indices[2] = {0, 1};
    cplInterface.readData(meshName, dataName, indices, maxDt, valueData);
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

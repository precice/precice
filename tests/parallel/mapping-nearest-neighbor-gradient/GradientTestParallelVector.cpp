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

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(MappingNearestNeighborGradient)

// Bidirectional test : Read: Vector & NNG - Write: Scalar & NN (Parallel Coupling)
BOOST_AUTO_TEST_CASE(GradientTestParallelVector)
{

  PRECICE_TEST("SolverOne"_on(3_ranks), "SolverTwo"_on(1_rank));

  if (context.isNamed("SolverOne")) {
    Participant interface(context.name, context.config(), context.rank, context.size);
    auto        meshName  = "MeshOne";
    auto        dataName1 = "Data1";
    auto        dataName2 = "Data2";

    int    vertexIDs[2];
    double xCoord       = context.rank * 0.4 + 0.05;
    double positions[4] = {xCoord, 0.0, xCoord + 0.2, 0.0};
    interface.setMeshVertices(meshName, positions, vertexIDs);
    BOOST_TEST(interface.requiresGradientDataFor(meshName, dataName1) == false);
    BOOST_TEST(interface.requiresGradientDataFor(meshName, dataName2) == false);
    interface.initialize();
    BOOST_TEST(interface.requiresGradientDataFor(meshName, dataName1) == false);
    BOOST_TEST(interface.requiresGradientDataFor(meshName, dataName2) == false);
    Eigen::Vector4d values;
    interface.advance(1.0);
    interface.readData(meshName, dataName2, vertexIDs, interface.getMaxTimeStepSize(), values);
    Eigen::Vector4d expected(context.rank * 2.0 + 1.0 + 0.05, context.rank * 2.0 + 1.0 + 0.05,
                             2.0 * (context.rank + 1) + 0.05, 2.0 * (context.rank + 1) + 0.05);
    BOOST_TEST(values == expected);
    interface.finalize();
  } else {
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    Participant interface(context.name, context.config(), context.rank, context.size);
    auto        meshName = "MeshTwo";
    int         vertexIDs[6];
    double      positions[12] = {0.0, 0.0, 0.2, 0.0, 0.4, 0.0, 0.6, 0.0, 0.8, 0.0, 1.0, 0.0};
    interface.setMeshVertices(meshName, positions, vertexIDs);

    auto dataName1 = "Data1";
    auto dataName2 = "Data2";
    BOOST_TEST(interface.requiresGradientDataFor(meshName, dataName1) == false);
    BOOST_TEST(interface.requiresGradientDataFor(meshName, dataName2) == true);

    interface.initialize();
    double values[12] = {1.0, 1.0,
                         2.0, 2.0,
                         3.0, 3.0,
                         4.0, 4.0,
                         5.0, 5.0,
                         6.0, 6.0};

    interface.writeData(meshName, dataName2, vertexIDs, values);

    BOOST_TEST(interface.requiresGradientDataFor(meshName, dataName1) == false);
    BOOST_TEST(interface.requiresGradientDataFor(meshName, dataName2) == true);

    if (interface.requiresGradientDataFor(meshName, dataName2)) {
      std::vector<double> gradientValues(24, 1.0);
      interface.writeGradientData(meshName, dataName2, vertexIDs, gradientValues);
    }
    interface.advance(1.0);
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI

#ifndef PRECICE_NO_MPI

#include <Eigen/Core>
#include <algorithm>
#include <deque>
#include <fstream>
#include <istream>
#include <vector>

#include "logging/LogMacros.hpp"
#include "math/constants.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Utils.hpp"
#include "precice/impl/ParticipantImpl.hpp"
#include "precice/impl/ParticipantState.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "precice/impl/Types.hpp"
#include "precice/precice.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using precice::testing::TestContext;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(MappingNearestNeighborGradient)

// Bidirectional test : Read: Vector & NNG - Write: Scalar & NN
PRECICE_TEST_SETUP("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks))
BOOST_AUTO_TEST_CASE(GradientTestParallelWriteVector)
{
  PRECICE_TEST();

  if (context.isNamed("SolverOne")) {
    Participant interface(context.name, context.config(), context.rank, context.size);
    auto        meshName = "MeshOne";
    auto        dataName = "Data2";

    if (context.isPrimary()) {
      std::vector<int>    vertexIDs(2);
      std::vector<double> positions = {1.0, 1.0, 2.0, 2., 2., 3.0};
      interface.setMeshVertices(meshName, positions, vertexIDs);
      interface.initialize();
      Eigen::Vector3d values;
      interface.advance(1.0);
      interface.readData(meshName, dataName, {&vertexIDs[0], 1}, interface.getMaxTimeStepSize(), values);
      Eigen::Vector3d expected(21.1, 24.8, 28.5);
      BOOST_TEST(interface.requiresGradientDataFor(meshName, dataName) == false);
      BOOST_TEST(testing::equals(values, expected));
      interface.readData(meshName, dataName, {&vertexIDs[1], 1}, interface.getMaxTimeStepSize(), values);
      Eigen::Vector3d expected2(2.3, 4.2, 6.1);
      BOOST_TEST(testing::equals(values, expected2));
    } else {
      std::vector<int>    vertexIDs(1);
      std::vector<double> positions = {4.0, 4.0, 4.0};
      interface.setMeshVertices(meshName, positions, vertexIDs);
      interface.initialize();
      Eigen::Vector3d values;
      interface.advance(1.0);
      interface.readData(meshName, dataName, vertexIDs, interface.getMaxTimeStepSize(), values);
      Eigen::Vector3d expected(1., 2., 3.);
      BOOST_TEST(interface.requiresGradientDataFor(meshName, dataName) == false);
      BOOST_TEST(testing::equals(values, expected));
    }
    interface.finalize();
  } else {
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    Participant interface(context.name, context.config(), context.rank, context.size);
    auto        meshName = "MeshTwo";
    if (context.isPrimary()) {
      std::vector<int>    vertexIDs(4);
      std::vector<double> positions = {4.0, 4.0, 4.0, 0.0, 0.4, 0.0, 0.7, 0.7, 1.7, 0.0, 1.0, 0.0};
      interface.setMeshVertices(meshName, positions, vertexIDs);
      interface.initialize();
      auto                dataName = "Data2";
      std::vector<double> values   = {1.0, 2.0, 3.0,
                                    -1.0, -1.0, -1.0,
                                    4.0, 5.0, 6.0,
                                    0.0, 0.0, 0.0};

      interface.writeData(meshName, dataName, vertexIDs, values);

      BOOST_TEST(interface.requiresGradientDataFor(meshName, dataName) == true);

      if (interface.requiresGradientDataFor(meshName, dataName)) {
        std::vector<double> gradients;
        for (unsigned int i = 0; i < 36; ++i) {
          gradients.emplace_back(i);
        }
        interface.writeGradientData(meshName, dataName, vertexIDs, gradients);
      }
    } else {
      // Assigned to the first rank
      std::vector<int>    vertexIDs(1);
      std::vector<double> positions = {2.1, 2.1, 3.1};
      interface.setMeshVertices(meshName, positions, vertexIDs);
      interface.initialize();
      auto                dataName = "Data2";
      std::vector<double> values   = {2.0, 3.0, 4.0};

      interface.writeData(meshName, dataName, vertexIDs, values);

      BOOST_TEST(interface.requiresGradientDataFor(meshName, dataName) == true);

      if (interface.requiresGradientDataFor(meshName, dataName)) {
        std::vector<double> gradients;
        for (int i = 0; i < 9; ++i) {
          gradients.emplace_back(-i);
        }
        interface.writeGradientData(meshName, dataName, {&vertexIDs[0], 1}, gradients);
      }
    }
    interface.advance(1.0);
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI

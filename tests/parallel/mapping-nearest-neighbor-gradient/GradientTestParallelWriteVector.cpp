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
#include "precice/SolverInterface.hpp"
#include "precice/impl/Participant.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"
#include "precice/types.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using precice::testing::TestContext;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(MappingNearestNeighborGradient)

// Bidirectional test : Read: Vector & NNG - Write: Scalar & NN
BOOST_AUTO_TEST_CASE(GradientTestParallelWriteVector)
{

  PRECICE_TEST("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks));

  if (context.isNamed("SolverOne")) {
    SolverInterface interface(context.name, context.config(), context.rank, context.size);
    auto            meshName = "MeshOne";
    auto            dataName = "Data2";

    std::vector<int> vertexIDs(2);
    if (context.isPrimary()) {
      std::vector<double> positions = {1.0, 1.0, 2.0, 2., 2., 3.0};
      interface.setMeshVertices(meshName, 2, positions.data(), vertexIDs.data());
      interface.initialize();
      Eigen::Vector3d values;
      double          preciceDt = interface.advance(1.0);
      interface.readBlockVectorData(meshName, dataName, 1, &vertexIDs[0], preciceDt, values.data());
      Eigen::Vector3d expected(21.1, 24.8, 28.5);
      BOOST_TEST(interface.requiresGradientDataFor(meshName, dataName) == false);
      BOOST_TEST(testing::equals(values, expected));
      interface.readBlockVectorData(meshName, dataName, 1, &vertexIDs[1], preciceDt, values.data());
      Eigen::Vector3d expected2(2.3, 4.2, 6.1);
      BOOST_TEST(testing::equals(values, expected2));
    } else {
      std::vector<double> positions = {4.0, 4.0, 4.0};
      interface.setMeshVertices(meshName, 1, positions.data(), vertexIDs.data());
      interface.initialize();
      Eigen::Vector3d values;
      double          preciceDt = interface.advance(1.0);
      interface.readBlockVectorData(meshName, dataName, 1, vertexIDs.data(), preciceDt, values.data());
      Eigen::Vector3d expected(1., 2., 3.);
      BOOST_TEST(interface.requiresGradientDataFor(meshName, dataName) == false);
      BOOST_TEST(testing::equals(values, expected));
    }
    interface.finalize();
  } else {
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    SolverInterface  interface(context.name, context.config(), context.rank, context.size);
    auto             meshName = "MeshTwo";
    std::vector<int> vertexIDs(4);
    if (context.isPrimary()) {
      std::vector<double> positions = {4.0, 4.0, 4.0, 0.0, 0.4, 0.0, 0.7, 0.7, 1.7, 0.0, 1.0, 0.0};
      interface.setMeshVertices(meshName, 4, positions.data(), vertexIDs.data());
      interface.initialize();
      auto                dataName = "Data2";
      std::vector<double> values   = {1.0, 2.0, 3.0,
                                    -1.0, -1.0, -1.0,
                                    4.0, 5.0, 6.0,
                                    0.0, 0.0, 0.0};

      interface.writeBlockVectorData(meshName, dataName, 4, vertexIDs.data(), values.data());

      BOOST_TEST(interface.requiresGradientDataFor(meshName, dataName) == true);

      if (interface.requiresGradientDataFor(meshName, dataName)) {
        std::vector<double> gradientValues;
        for (unsigned int i = 0; i < 36; ++i) {
          gradientValues.emplace_back(i);
        }
        interface.writeBlockVectorGradientData(meshName, dataName, 4, vertexIDs.data(), gradientValues.data());
      }
    } else {
      // Assigned to the first rank
      std::vector<double> positions = {2.1, 2.1, 3.1};
      interface.setMeshVertices(meshName, 1, positions.data(), vertexIDs.data());
      interface.initialize();
      auto                dataName = "Data2";
      std::vector<double> values   = {2.0, 3.0, 4.0};

      interface.writeBlockVectorData(meshName, dataName, 1, vertexIDs.data(), values.data());

      BOOST_TEST(interface.requiresGradientDataFor(meshName, dataName) == true);

      if (interface.requiresGradientDataFor(meshName, dataName)) {
        std::vector<double> gradientValues;
        for (int i = 0; i < 9; ++i) {
          gradientValues.emplace_back(-i);
        }
        interface.writeVectorGradientData(meshName, dataName, vertexIDs[0], gradientValues.data());
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

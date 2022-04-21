#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>
#include "helpers.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(DirectMeshAccess)
BOOST_AUTO_TEST_CASE(AccessReceivedMeshNoOverlapTwoLevelInit)
{
  // Test case for parallel mesh partitioning without any mapping. Each solver
  // runs on two ranks. SolverTwo defines 5(2 and 3) vertices which need to be
  // repartitioned on SolverOne according to the defined boundingBoxes
  // (resulting in 3 and 2 vertices per rank. The boundingBoxes don't have any
  // overlap.
  // Config using the two-level-initialization
  PRECICE_TEST("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks));
  const std::vector<double> boundingBoxSecondary      = std::vector<double>{0.0, 1.0, 4.0, 7};
  const std::vector<double> expectedPositionSecondary = std::vector<double>{0.0, 4.0, 0.0, 5.0};
  const std::vector<double> writeDataSecondary        = std::vector<double>({4, 5});
  const std::vector<double> expectedReadDataSecondary = std::vector<double>({3, 4, 5});
  runTestAccessReceivedMesh(context, boundingBoxSecondary, writeDataSecondary, expectedPositionSecondary, expectedReadDataSecondary, 0);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // DirectMeshAccess

#endif // PRECICE_NO_MPI

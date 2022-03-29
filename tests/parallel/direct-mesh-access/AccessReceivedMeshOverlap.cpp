#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>
#include "helpers.hpp"

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(DirectMeshAccess)
BOOST_AUTO_TEST_CASE(AccessReceivedMeshOverlap)
{

  // Test case for parallel mesh partitioning without any mapping. Each solver
  // runs on two ranks. SolverTwo defines 5(2 and 3) vertices which need to be
  // repartitioned on SolverOne according to the defined boundingBoxes
  // (resulting in 3 and 3 vertices per rank). The boundingBoxes of the other
  // participant have an overlap including the common vertex with position 3.
  // The vertex is 'written' by both ranks and summed up on the receiver side.  BOOST_TEST(false);

  PRECICE_TEST("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks));
  const std::vector<double> boundingBoxSlave      = std::vector<double>{0.0, 1.0, 3.0, 7};
  const std::vector<double> expectedPositionSlave = std::vector<double>{0.0, 3.0, 0.0, 4.0, 0.0, 5.0};
  const std::vector<double> writeDataSlave        = std::vector<double>({4, 5, 6});
  const std::vector<double> expectedReadDataSlave = std::vector<double>({7, 5, 6});
  runTestAccessReceivedMesh(context, boundingBoxSlave, writeDataSlave, expectedPositionSlave, expectedReadDataSlave, 0);
}

BOOST_AUTO_TEST_SUITE_END() // PreciceTests
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // DirectMeshAccess

#endif // PRECICE_NO_MPI

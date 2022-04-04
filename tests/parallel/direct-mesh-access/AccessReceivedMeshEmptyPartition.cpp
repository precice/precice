#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>
#include "helpers.hpp"

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(DirectMeshAccess)
BOOST_AUTO_TEST_CASE(AccessReceivedMeshEmptyPartition)
{
  // Same as above, but one rank of SolverOne does not receive any vertices due
  // to the defined bounding box. Test case for parallel mesh partitioning without
  // any mapping. Each solver runs on two ranks. SolverTwo defines 5(2 and 3)
  // vertices which need to be repartitioned on SolverOne according to the defined
  // boundingBoxes (resulting in 3  vertices on one rank and 2 completely filtered
  // vertices). Filtered vertices are filled with zero data values

  PRECICE_TEST("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks));
  const std::vector<double> boundingBoxSlave      = std::vector<double>{10.0, 10.0, 13.0, 17};
  const std::vector<double> expectedPositionSlave = std::vector<double>{};
  const std::vector<double> writeDataSlave        = std::vector<double>({});
  const std::vector<double> expectedReadDataSlave = std::vector<double>({3., 0., 0.});
  runTestAccessReceivedMesh(context, boundingBoxSlave, writeDataSlave, expectedPositionSlave, expectedReadDataSlave, 0);
}

BOOST_AUTO_TEST_SUITE_END() // PreciceTests
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // DirectMeshAccess

#endif // PRECICE_NO_MPI

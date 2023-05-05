#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include "helpers.hpp"

using namespace precice;
using precice::testing::TestContext;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MappingNearestNeighborGradient)

// Unidirectional Nearest Neighbor Gradient Read Mapping
BOOST_AUTO_TEST_CASE(GradientTestUnidirectionalReadVector)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank))
  testVectorGradientFunctions(context);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI

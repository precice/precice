#ifndef PRECICE_NO_MPI

#include "helper.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(MapInitialData)
BOOST_AUTO_TEST_CASE(NonZeroData)
{
  PRECICE_TEST("One"_on(1_rank), "Two"_on(2_ranks));

  testMapInitialDataP(context, 1.0, 1.0, 2, 0);
}

BOOST_AUTO_TEST_SUITE_END() // MapInitialData
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

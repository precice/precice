#ifndef PRECICE_NO_MPI

#include "helper.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(MapInitialData)
PRECICE_TEST_SETUP("One"_on(1_rank), "Two"_on(2_ranks))
BOOST_AUTO_TEST_CASE(ZeroData)
{
  PRECICE_TEST();

  testMapInitialDataP(context, 0.0, 0.0, 1, 0); // End of time window is not zero
}

BOOST_AUTO_TEST_SUITE_END() // MapInitialData
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

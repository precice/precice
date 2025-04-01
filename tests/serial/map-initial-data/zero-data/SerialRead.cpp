#ifndef PRECICE_NO_MPI

#include "../helper.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MapInitialData)
BOOST_AUTO_TEST_SUITE(ZeroData)
PRECICE_TEST_SETUP("One"_on(1_rank), "Two"_on(1_rank))
BOOST_AUTO_TEST_CASE(SerialRead)
{
  PRECICE_TEST();

  testMapInitialData(context, 0.0, 0.0, 1, 0); // End of time window is not zero
}

BOOST_AUTO_TEST_SUITE_END() // NonZeroData
BOOST_AUTO_TEST_SUITE_END() // MapInitialData
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

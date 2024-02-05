#ifndef PRECICE_NO_MPI

#include "../helper.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MapInitialData)
BOOST_AUTO_TEST_SUITE(NonZeroData)
BOOST_AUTO_TEST_CASE(SerialRead)
{
  PRECICE_TEST("One"_on(1_rank), "Two"_on(1_rank));

  testMapInitialData(context, 1.0, 1.0, 2, 0);
}

BOOST_AUTO_TEST_SUITE_END() // NonZeroData
BOOST_AUTO_TEST_SUITE_END() // MapInitialData
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

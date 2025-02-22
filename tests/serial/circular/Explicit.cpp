#ifndef PRECICE_NO_MPI

#include "helper.hpp"
#include "testing/Testing.hpp"

#include <precice/precice.hpp>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Circular)
PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank), "C"_on(1_rank))
BOOST_AUTO_TEST_CASE(Explicit)
{
  PRECICE_TEST();
  tests::cyclicExplicit(context);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Circular

#endif // PRECICE_NO_MPI

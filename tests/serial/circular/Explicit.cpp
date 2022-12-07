#ifndef PRECICE_NO_MPI

#include "helper.hpp"
#include "testing/Testing.hpp"

#include <precice/precice.hpp>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Circular)
BOOST_AUTO_TEST_CASE(Explicit)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), "C"_on(1_rank));
  tests::cyclicExplicit(context);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Circular

#endif // PRECICE_NO_MPI

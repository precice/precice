#ifndef PRECICE_NO_MPI

#include "../helper.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Remeshing)
BOOST_AUTO_TEST_SUITE(ParallelImplicit)
BOOST_AUTO_TEST_SUITE(Acceleration)
BOOST_AUTO_TEST_SUITE(Aitken)
PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank))
BOOST_AUTO_TEST_CASE(RemeshFirstSerial)
{
  PRECICE_TEST();
  precice::tests::remesh::parallelImplicit::acceleration::runResetAAitken(context);
}

BOOST_AUTO_TEST_SUITE_END() // Convergence
BOOST_AUTO_TEST_SUITE_END() // Constant
BOOST_AUTO_TEST_SUITE_END() // ParallelImplicit
BOOST_AUTO_TEST_SUITE_END() // Remeshing
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

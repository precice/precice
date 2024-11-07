#ifndef PRECICE_NO_MPI

#include "../helper.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Remeshing)
BOOST_AUTO_TEST_SUITE(ParallelImplicit)
BOOST_AUTO_TEST_SUITE(Acceleration)
BOOST_AUTO_TEST_SUITE(IQNIMVJ)
BOOST_AUTO_TEST_CASE(RemeshSecondSerial)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank));
  precice::tests::remesh::parallelImplicit::acceleration::runResetAIQNIMVJ(context);
}

BOOST_AUTO_TEST_SUITE_END() // Convergence
BOOST_AUTO_TEST_SUITE_END() // IQNIMVJ
BOOST_AUTO_TEST_SUITE_END() // ParallelImplicit
BOOST_AUTO_TEST_SUITE_END() // Remeshing
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

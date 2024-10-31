#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include "../helper.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Remeshing)
BOOST_AUTO_TEST_SUITE(ParallelImplicit)
BOOST_AUTO_TEST_SUITE(Noop)
BOOST_AUTO_TEST_CASE(RemeshSecond2LI)
{
  PRECICE_TEST("A"_on(2_ranks), "B"_on(2_ranks));
  precice::tests::remesh::parallelImplicit::noop::runResetA(context);
}

BOOST_AUTO_TEST_SUITE_END() // Noop
BOOST_AUTO_TEST_SUITE_END() // ParallelImplicit
BOOST_AUTO_TEST_SUITE_END() // Remeshing
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

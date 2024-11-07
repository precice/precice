#ifndef PRECICE_NO_MPI

#include "../helper.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Remeshing)
BOOST_AUTO_TEST_SUITE(ParallelExplicit)
BOOST_AUTO_TEST_SUITE(Noop)
BOOST_AUTO_TEST_CASE(RemeshBoth2LI)
{
  PRECICE_TEST("A"_on(2_ranks), "B"_on(2_ranks));
  precice::tests::remesh::parallelExplicit::noop::runResetBoth(context);
}

BOOST_AUTO_TEST_SUITE_END() // Noop
BOOST_AUTO_TEST_SUITE_END() // ParallelExplicit
BOOST_AUTO_TEST_SUITE_END() // Remeshing
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

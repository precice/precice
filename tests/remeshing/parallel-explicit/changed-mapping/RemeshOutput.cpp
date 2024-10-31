#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include "../helper.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Remeshing)
BOOST_AUTO_TEST_SUITE(ParallelExplicit)
BOOST_AUTO_TEST_SUITE(ChangedMapping)
BOOST_AUTO_TEST_CASE(RemeshOutput)
{
  PRECICE_TEST("A"_on(2_ranks), "B"_on(2_ranks));
  precice::tests::remesh::parallelExplicit::changemapping::runResetOutput(context);
}

BOOST_AUTO_TEST_SUITE_END() // ChangedMapping
BOOST_AUTO_TEST_SUITE_END() // ParallelExplicit
BOOST_AUTO_TEST_SUITE_END() // Remeshing
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

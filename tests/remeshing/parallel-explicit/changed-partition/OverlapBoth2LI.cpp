#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include "../helper.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Remeshing)
BOOST_AUTO_TEST_SUITE(ParallelExplicit)
BOOST_AUTO_TEST_SUITE(ChangedPartition)
PRECICE_TEST_SETUP("A"_on(2_ranks), "B"_on(2_ranks))
BOOST_AUTO_TEST_CASE(OverlapBoth2LI)
{
  PRECICE_TEST();
  precice::tests::remesh::parallelExplicit::changepartition::runOverlapBoth(context);
}

BOOST_AUTO_TEST_SUITE_END() // ChangedPartition
BOOST_AUTO_TEST_SUITE_END() // ParallelExplicit
BOOST_AUTO_TEST_SUITE_END() // Remeshing
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

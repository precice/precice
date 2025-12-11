#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

// Test representing the minimal lifecylce, which consists out of construction only.
// The destructor has to cleanup correctly.
BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(Lifecycle)
PRECICE_TEST_SETUP("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks))
BOOST_AUTO_TEST_CASE(ConstructOnly)
{
  PRECICE_TEST();

  precice::Participant interface(context.name, context.config(), context.rank, context.size);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // Lifecycle

#endif // PRECICE_NO_MPI

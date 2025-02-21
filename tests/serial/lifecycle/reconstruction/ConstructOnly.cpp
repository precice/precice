#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Lifecycle)
BOOST_AUTO_TEST_SUITE(Reconstruction)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(ConstructOnly)
{
  PRECICE_TEST();
  for (auto n : {1, 2, 3})
    BOOST_TEST_CONTEXT("construction #" << n)
    {
      precice::Participant interface(context.name, context.config(), context.rank, context.size, context.comm());
    }
}

BOOST_AUTO_TEST_SUITE_END() // Reconstruction
BOOST_AUTO_TEST_SUITE_END() // Lifecycle
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

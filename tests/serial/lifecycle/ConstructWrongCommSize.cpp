#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Lifecycle)
BOOST_AUTO_TEST_CASE(ConstructWrongCommSize)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  // tests non-matching communicator sizes
  BOOST_CHECK_THROW((precice::Participant{context.name, context.config(), context.rank, 3}), ::precice::Error);

  BOOST_CHECK_THROW((precice::Participant{context.name, context.config(), context.rank, 5, context.comm()}), ::precice::Error);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Lifecycle

#endif // PRECICE_NO_MPI

#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Lifecycle)
BOOST_AUTO_TEST_CASE(ConstructNullComm)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  // tests passing null as communicator
  BOOST_CHECK_THROW((precice::Participant{context.name, context.config(), context.rank, 3, nullptr}), ::precice::Error);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Lifecycle

#endif // PRECICE_NO_MPI

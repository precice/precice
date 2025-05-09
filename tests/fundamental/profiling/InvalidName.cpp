#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Fundamental)
BOOST_AUTO_TEST_SUITE(Profiling)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank));
BOOST_AUTO_TEST_CASE(InvalidName)
{
  PRECICE_TEST();

  precice::Participant p(context.name, context.config(), context.rank, context.size);

  BOOST_CHECK_THROW(p.startProfilingSection("not/allowed"), ::precice::Error);
}

BOOST_AUTO_TEST_SUITE_END() // Profiling
BOOST_AUTO_TEST_SUITE_END() // Fundamental
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

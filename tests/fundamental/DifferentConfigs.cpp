#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Fundamental)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(DifferentConfigs)
{
  PRECICE_TEST();

  auto config = context.prefix(precice::testing::getTestName());
  if (context.isNamed("SolverOne")) {
    config.append("-A.xml");
  } else {
    config.append("-B.xml");
  }

  precice::Participant p(context.name, config, context.rank, context.size);

  BOOST_CHECK_THROW(p.initialize(), ::precice::Error);
}

BOOST_AUTO_TEST_SUITE_END() // Fundamental
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(TimeHandling)
BOOST_AUTO_TEST_CASE(SimpleMaxTime, *boost::unit_test::tolerance(1e-18))
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  precice::Participant p(context.name, context.config(), context.rank, context.size);
  double               coord[] = {0, 0};
  p.setMeshVertex(context.name + "-Mesh", coord);
  p.initialize();

  const double expected = 0.01;

  BOOST_TEST(p.getMaxTimeStepSize() == expected);
  BOOST_TEST(p.isCouplingOngoing());
  p.advance(p.getMaxTimeStepSize());

  BOOST_TEST(p.getMaxTimeStepSize() == expected);
  BOOST_TEST(p.isCouplingOngoing());
  p.advance(p.getMaxTimeStepSize());

  BOOST_REQUIRE(p.isCouplingOngoing());
  // BUG was: return value 0.009999999999999998
  BOOST_TEST(p.getMaxTimeStepSize() == expected);
  p.advance(p.getMaxTimeStepSize());

  BOOST_TEST(p.isCouplingOngoing());
  // BUG was: return value 0.010000000000000002
  BOOST_TEST(p.getMaxTimeStepSize() == expected);
  p.advance(p.getMaxTimeStepSize());

  BOOST_TEST(p.getMaxTimeStepSize() == 0.);
  BOOST_REQUIRE(!p.isCouplingOngoing());
}

BOOST_AUTO_TEST_SUITE_END() // TimeHandling
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI

#include "precice/SolverInterface.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using precice::testing::TestContext;

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(ExceptionTests)

BOOST_AUTO_TEST_CASE(ConfigurationFileNotFound)
{
  PRECICE_TEST(1_rank);

  BOOST_REQUIRE_THROW(SolverInterface s("a", "b", 0, 1), ::precice::Error);
}

BOOST_AUTO_TEST_CASE(WrongParticipant)
{
  PRECICE_TEST(1_rank);

  BOOST_REQUIRE_THROW(precice::SolverInterface s("a", testing::getPathToSources() + "/precice/tests/config-checker.xml", 0, 1), ::precice::Error);
}

BOOST_AUTO_TEST_CASE(FinalizeTwice)
{
  PRECICE_TEST(1_rank);
  precice::SolverInterface s("SolverOne", testing::getPathToSources() + "/precice/tests/config-checker.xml", 0, 1);
  s.finalize();
  BOOST_REQUIRE_THROW(s.finalize(), ::precice::Error); // this should be a preCICE exception type!
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

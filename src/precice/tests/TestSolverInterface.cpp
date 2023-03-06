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
  //BOOST_REQUIRE_THROW does not work here. Expanding it manually.
  bool allFine = false;
  try {
    SolverInterface s("a", "b", 0, 1);
  } catch (::precice::Error) { // this should be a preCICE exception type!
    allFine = true;
  }
  BOOST_TEST(allFine);
}

BOOST_AUTO_TEST_CASE(WrongParticipant)
{
  PRECICE_TEST(1_rank);
  //BOOST_REQUIRE_THROW does not work here. Expanding it manually.
  bool allFine = false;
  try {
    precice::SolverInterface s("a", testing::getPathToSources() + "/precice/tests/config-checker.xml", 0, 1);
  } catch (::precice::Error) { // this should be a preCICE exception type!
    allFine = true;
  }
  BOOST_TEST(allFine);
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

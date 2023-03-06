#include "precice/SolverInterface.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using precice::testing::TestContext;

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(ExceptionTests)

BOOST_AUTO_TEST_CASE(ConfigurationFileNotFound)
{
  //BOOST_REQUIRE_THROW does not work here. Expanding it manually.
  bool allFine = false;
  try {
    SolverInterface s("a", "b", 0, 1);
  } catch (std::exception) { // this should be a preCICE exception type!
    allFine = true;
  }
  BOOST_TEST(allFine);
}

BOOST_AUTO_TEST_CASE(WrongParticipant)
{
  //BOOST_REQUIRE_THROW does not work here. Expanding it manually.
  bool allFine = false;
  try {
    precice::SolverInterface s("a", testing::getPathToSources() + "/precice/tests/aitken-acceleration.xml", 0, 1);
  } catch (std::exception) { // this should be a preCICE exception type!
    allFine = true;
  }
  BOOST_TEST(allFine);
}

BOOST_AUTO_TEST_CASE(FinalizeTwice)
{
  precice::SolverInterface s("A", testing::getPathToSources() + "/precice/tests/aitken-acceleration.xml", 0, 1);
  s.finalize();
  BOOST_REQUIRE_THROW(s.finalize(), std::exception); // this should be a preCICE exception type!
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

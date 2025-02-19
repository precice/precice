#include <sstream>
#include <string>
#include "precice/Tooling.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(Tooling)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(MarkdownReference)
{
  PRECICE_TEST();
  const std::string ref = [] {
    std::ostringstream oss;
    precice::tooling::printConfigReference(oss, precice::tooling::ConfigReferenceType::MD);
    return oss.str();
  }();

  BOOST_TEST(ref.size() > 0);
  for (const auto &keyword : {"# precice-configuration", "<precice-configuration", "</precice-configuration>", "Example", "Valid Subtags:", "Attribute"}) {
    BOOST_TEST(ref.find(keyword) != std::string::npos, "The output should include \"" << keyword << '"');
  }
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(XMLReference)
{
  PRECICE_TEST();
  const std::string ref = [] {
    std::ostringstream oss;
    precice::tooling::printConfigReference(oss, precice::tooling::ConfigReferenceType::XML);
    return oss.str();
  }();

  BOOST_TEST(ref.size() > 0);
  for (const auto &keyword : {"<!-- TAG precice-configuration", "<!-- TAG mesh", "ATTR name:"}) {
    BOOST_TEST(ref.find(keyword) != std::string::npos, "The output should include \"" << keyword << '"');
  }
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(DTDReference)
{
  PRECICE_TEST();
  const std::string ref = [] {
    std::ostringstream oss;
    precice::tooling::printConfigReference(oss, precice::tooling::ConfigReferenceType::DTD);
    return oss.str();
  }();

  BOOST_TEST(ref.size() > 0);
  for (const auto &keyword : {"<!ELEMENT precice-configuration", "<!ATTLIST mesh name"}) {
    BOOST_TEST(ref.find(keyword) != std::string::npos, "The output should include \"" << keyword << '"');
  }
}

BOOST_AUTO_TEST_SUITE(ConfigCheck)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(Serial)
{
  PRECICE_TEST();
  BOOST_REQUIRE_NO_THROW(
      precice::tooling::checkConfiguration(
          precice::testing::getPathToSources() + "/precice/tests/config-checker.xml",
          "SolverTwo",
          1));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(Parallel)
{
  PRECICE_TEST();
  BOOST_REQUIRE_NO_THROW(
      precice::tooling::checkConfiguration(
          precice::testing::getPathToSources() + "/precice/tests/config-checker.xml",
          "SolverTwo",
          4));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

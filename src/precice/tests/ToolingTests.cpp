#include <sstream>
#include <string>
#include "precice/Tooling.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(Tooling)

BOOST_AUTO_TEST_CASE(MarkdownReference)
{
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

BOOST_AUTO_TEST_CASE(XMLReference)
{
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

BOOST_AUTO_TEST_CASE(DTDReference)
{
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

BOOST_AUTO_TEST_CASE(Serial)
{
  BOOST_REQUIRE_NO_THROW(
      precice::tooling::checkConfiguration(
          precice::testing::getPathToSources() + "/precice/tests/config-checker.xml",
          "SolverTwo",
          1));
}

BOOST_AUTO_TEST_CASE(Parallel)
{
  BOOST_REQUIRE_NO_THROW(
      precice::tooling::checkConfiguration(
          precice::testing::getPathToSources() + "/precice/tests/config-checker.xml",
          "SolverTwo",
          4));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

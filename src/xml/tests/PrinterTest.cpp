#include <sstream>
#include <string>
#include "precice/config/Configuration.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "xml/Printer.hpp"

BOOST_AUTO_TEST_SUITE(XML)
BOOST_AUTO_TEST_SUITE(Printer)

BOOST_AUTO_TEST_CASE(Documentation)
{
  PRECICE_TEST(1_rank);
  std::ostringstream             oss;
  precice::config::Configuration config;
  precice::xml::toDocumentation(oss, config.getXMLTag());
  BOOST_TEST(!oss.str().empty());
}

BOOST_AUTO_TEST_CASE(DTD)
{
  PRECICE_TEST(1_rank);
  std::ostringstream             oss;
  precice::config::Configuration config;
  precice::xml::toDTD(oss, config.getXMLTag());
  BOOST_TEST(!oss.str().empty());
}

BOOST_AUTO_TEST_CASE(Markdown)
{
  PRECICE_TEST(1_rank);
  std::ostringstream             oss;
  precice::config::Configuration config;
  precice::xml::toMarkdown(oss, config.getXMLTag());
  BOOST_TEST(!oss.str().empty());
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

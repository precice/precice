#include "testing/Testing.hpp"
#include "xml/Printer.hpp"
#include "precice/config/Configuration.hpp"

#include <sstream>

BOOST_AUTO_TEST_SUITE(XML)
BOOST_AUTO_TEST_SUITE(Printer, *precice::testing::OnMaster())

BOOST_AUTO_TEST_CASE(Documentation)
{
  std::ostringstream             oss;
  precice::config::Configuration config;
  precice::xml::toDocumentation(oss, config.getXMLTag());
  BOOST_TEST(!oss.str().empty());
}

BOOST_AUTO_TEST_CASE(DTD)
{
  std::ostringstream             oss;
  precice::config::Configuration config;
  precice::xml::toDTD(oss, config.getXMLTag());
  BOOST_TEST(!oss.str().empty());
}

BOOST_AUTO_TEST_CASE(Markdown)
{
  std::ostringstream             oss;
  precice::config::Configuration config;
  precice::xml::toMarkdown(oss, config.getXMLTag());
  BOOST_TEST(!oss.str().empty());
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

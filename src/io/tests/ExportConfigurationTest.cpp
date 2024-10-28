#include <list>
#include "io/ExportContext.hpp"
#include "io/SharedPointer.hpp"
#include "io/config/ExportConfiguration.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "xml/XMLTag.hpp"

BOOST_AUTO_TEST_SUITE(IOTests)
BOOST_AUTO_TEST_SUITE(Configuration)

using namespace precice;

BOOST_AUTO_TEST_CASE(VTKEvery10)
{
  PRECICE_TEST(1_rank);
  using xml::XMLTag;
  XMLTag                  tag = xml::getRootTag();
  io::ExportConfiguration config(tag);
  xml::configure(tag, xml::ConfigurationContext{}, testing::getPathToSources() + "/io/tests/config1.xml");
  BOOST_TEST(config.exportContexts().size() == 1);
  const io::ExportContext &econtext = config.exportContexts().front();
  BOOST_TEST(econtext.type == "vtk");
  BOOST_TEST(econtext.everyNTimeWindows == 10);
}

BOOST_AUTO_TEST_CASE(VTKLocation)
{
  PRECICE_TEST(1_rank);
  using xml::XMLTag;
  XMLTag                  tag = xml::getRootTag();
  io::ExportConfiguration config(tag);
  xml::configure(tag, xml::ConfigurationContext{}, testing::getPathToSources() + "/io/tests/config2.xml");
  BOOST_TEST(config.exportContexts().size() == 1);
  const io::ExportContext &econtext = config.exportContexts().front();
  BOOST_TEST(econtext.type == "vtk");
  BOOST_TEST(econtext.everyNTimeWindows == 1);
  BOOST_TEST(econtext.location == "somepath");
}

BOOST_AUTO_TEST_SUITE_END() // Configuration
BOOST_AUTO_TEST_SUITE_END() // IOTests

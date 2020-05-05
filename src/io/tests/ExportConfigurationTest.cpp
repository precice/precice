#include "io/config/ExportConfiguration.hpp"
#include "testing/Testing.hpp"
#include "xml/XMLTag.hpp"

BOOST_AUTO_TEST_SUITE(IOTests)

using namespace precice;

BOOST_AUTO_TEST_CASE(Configuration)
{
  PRECICE_TEST(1_rank);
  using xml::XMLTag;
  XMLTag tag = xml::getRootTag();
  {
    io::ExportConfiguration config(tag);
    xml::configure(tag, xml::ConfigurationContext{}, testing::getPathToSources() + "/io/tests/config1.xml");
    BOOST_TEST(config.exportContexts().size() == 1);
    const io::ExportContext &context = config.exportContexts().front();
    BOOST_TEST(context.type == "vtk");
    BOOST_TEST(context.everyNTimeWindows == 10);
    BOOST_TEST(context.triggerSolverPlot);
  }
  {
    tag.clear();
    io::ExportConfiguration config(tag);
    xml::configure(tag, xml::ConfigurationContext{}, testing::getPathToSources() + "/io/tests/config2.xml");
    BOOST_TEST(config.exportContexts().size() == 1);
    const io::ExportContext &context = config.exportContexts().front();
    BOOST_TEST(context.type == "vtk");
    BOOST_TEST(context.everyNTimeWindows == 1);
    BOOST_TEST(context.location == "somepath");
    BOOST_TEST(not context.triggerSolverPlot);
  }
}

BOOST_AUTO_TEST_SUITE_END() // IOTests

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

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(VTKEvery10)
{
  PRECICE_TEST();
  using xml::XMLTag;
  XMLTag                  tag = xml::getRootTag();
  io::ExportConfiguration config(tag);
  xml::configure(tag, xml::ConfigurationContext{}, testing::getPathToSources() + "/io/tests/config1.xml");
  BOOST_TEST(config.exportContexts().size() == 1);
  const io::ExportContext &econtext = config.exportContexts().front();
  BOOST_TEST(econtext.type == "vtk");
  BOOST_TEST(econtext.everyNTimeWindows == 10);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(VTKLocation)
{
  PRECICE_TEST();
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

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(VTUSingleMesh)
{
  PRECICE_TEST();
  using xml::XMLTag;
  XMLTag                  tag = xml::getRootTag();
  io::ExportConfiguration config(tag);

  xml::configure(tag, xml::ConfigurationContext{},
                 testing::getPathToSources() + "/io/tests/config3.xml");

  BOOST_TEST(config.exportContexts().size() == 1);
  const io::ExportContext &econtext = config.exportContexts().front();

  BOOST_TEST(econtext.type == "vtu");
  BOOST_TEST(econtext.selectedMeshes.size() == 1);
  BOOST_TEST(econtext.selectedMeshes.front() == "SolidMesh");
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(VTUMultipleMeshes)
{
  PRECICE_TEST();
  using xml::XMLTag;
  XMLTag                  tag = xml::getRootTag();
  io::ExportConfiguration config(tag);

  xml::configure(tag, xml::ConfigurationContext{},
                 testing::getPathToSources() + "/io/tests/config4.xml");

  const io::ExportContext &econtext = config.exportContexts().front();

  BOOST_TEST(econtext.selectedMeshes.size() == 2);
  BOOST_TEST(econtext.selectedMeshes[0] == "SolidMesh");
  BOOST_TEST(econtext.selectedMeshes[1] == "FluidMesh");
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(VTUBothMeshAttributes)
{
  PRECICE_TEST();
  using xml::XMLTag;
  XMLTag                  tag = xml::getRootTag();
  io::ExportConfiguration config(tag);

  BOOST_CHECK_THROW(
      xml::configure(tag, xml::ConfigurationContext{},
                     testing::getPathToSources() + "/io/tests/config5.xml"),
      std::exception);
}

BOOST_AUTO_TEST_SUITE_END() // Configuration
BOOST_AUTO_TEST_SUITE_END() // IOTests

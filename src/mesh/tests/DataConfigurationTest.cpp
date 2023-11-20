#include <memory>
#include <string>
#include <vector>
#include "logging/Logger.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "xml/XMLTag.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(MeshTests)

BOOST_AUTO_TEST_CASE(DataConfig)
{
  PRECICE_TEST(1_rank);
  std::string filename(testing::getPathToSources() + "/mesh/tests/data-config.xml");
  using xml::XMLTag;
  XMLTag                  tag = xml::getRootTag();
  mesh::DataConfiguration dataConfig(tag);
  xml::configure(tag, xml::ConfigurationContext{}, filename);
  BOOST_TEST(dataConfig.data().size() == 3);
  BOOST_TEST(dataConfig.data().at(0).name == "vector-data");
  BOOST_TEST(dataConfig.data().at(1).name == "floating-data");
  BOOST_TEST(dataConfig.data().at(2).name == "second-vector-data");
}

BOOST_AUTO_TEST_SUITE_END()

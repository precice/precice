#include "mesh/config/DataConfiguration.hpp"
#include "testing/Testing.hpp"
#include "utils/xml/XMLTag.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(Mesh)

BOOST_AUTO_TEST_CASE(DataConfig)
{
  std::string filename(utils::Globals::getPathToSources() + "/mesh/boosttests/data-config.xml");
  int dim = 3;
  using utils::XMLTag;
  XMLTag tag = utils::getRootTag();
  mesh::DataConfiguration dataConfig(tag);
  dataConfig.setDimensions(dim);
  utils::configure(tag, filename);
  BOOST_TEST(dataConfig.data().size() == 3);
  BOOST_TEST(dataConfig.data()[0].name == "vector-data");
  BOOST_TEST(dataConfig.data()[0].dimensions == 3);
  BOOST_TEST(dataConfig.data()[1].name == "floating-data");
  BOOST_TEST(dataConfig.data()[1].dimensions == 1);
  BOOST_TEST(dataConfig.data()[2].name == "second-vector-data");
  BOOST_TEST(dataConfig.data()[2].dimensions == 3);
}

BOOST_AUTO_TEST_SUITE_END()

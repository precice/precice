#include <memory>
#include <string>
#include <vector>
#include "logging/Logger.hpp"
#include "mesh/config/GradientConfiguration.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "xml/XMLTag.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(MeshTests)

BOOST_AUTO_TEST_CASE(GradientConfig)
{
  PRECICE_TEST(1_rank);
  std::string filename(testing::getPathToSources() + "/mesh/tests/gradient-config.xml");
  int         dim = 3;
  using xml::XMLTag;
  XMLTag                  tag = xml::getRootTag();
  mesh::GradientConfiguration gradientConfig(tag);
  gradientConfig.setDimensions(dim);
  xml::configure(tag, xml::ConfigurationContext{}, filename);
  BOOST_TEST(gradientConfig.gradients().size() == 3);
  BOOST_TEST(gradientConfig.gradients()[0].name == "vector-data");
  BOOST_TEST(gradientConfig.gradients()[0].dimensions == 3);
  BOOST_TEST(gradientConfig.gradients()[1].name == "floating-data");
  BOOST_TEST(gradientConfig.gradients()[1].dimensions == 1);
  BOOST_TEST(gradientConfig.gradients()[2].name == "second-vector-data");
  BOOST_TEST(gradientConfig.gradients()[2].dimensions == 3);
}

BOOST_AUTO_TEST_SUITE_END()

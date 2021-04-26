#include <memory>
#include <string>
#include <vector>
#include "logging/Logger.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/config/MappingConfiguration.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "xml/XMLTag.hpp"

using namespace precice;
using namespace precice::mapping;

BOOST_AUTO_TEST_SUITE(MappingTests)
BOOST_AUTO_TEST_SUITE(Configuration)

BOOST_AUTO_TEST_CASE(Configuration)
{
  PRECICE_TEST(1_rank);

  std::string pathToTests = testing::getPathToSources() + "/mapping/tests/";
  std::string file(pathToTests + "mapping-config.xml");
  using xml::XMLTag;
  XMLTag                     tag = xml::getRootTag();
  mesh::PtrDataConfiguration dataConfig(new mesh::DataConfiguration(tag));
  dataConfig->setDimensions(3);
  mesh::PtrMeshConfiguration meshConfig(new mesh::MeshConfiguration(tag, dataConfig));
  meshConfig->setDimensions(3);
  mapping::MappingConfiguration mappingConfig(tag, meshConfig);
  xml::configure(tag, xml::ConfigurationContext{}, file);

  BOOST_TEST(meshConfig->meshes().size() == 3);
  BOOST_TEST(mappingConfig.mappings().size() == 3);
  BOOST_TEST(mappingConfig.mappings().at(0).timing == MappingConfiguration::ON_DEMAND);
  BOOST_TEST(mappingConfig.mappings().at(0).fromMesh == meshConfig->meshes().at(0));
  BOOST_TEST(mappingConfig.mappings().at(0).toMesh == meshConfig->meshes().at(2));
  BOOST_TEST(mappingConfig.mappings().at(0).direction == MappingConfiguration::WRITE);

  BOOST_TEST(mappingConfig.mappings().at(1).timing == MappingConfiguration::INITIAL);
  BOOST_TEST(mappingConfig.mappings().at(1).fromMesh == meshConfig->meshes().at(2));
  BOOST_TEST(mappingConfig.mappings().at(1).toMesh == meshConfig->meshes().at(1));
  BOOST_TEST(mappingConfig.mappings().at(1).direction == MappingConfiguration::READ);

  BOOST_TEST(mappingConfig.mappings().at(2).timing == MappingConfiguration::ON_ADVANCE);
  BOOST_TEST(mappingConfig.mappings().at(2).fromMesh == meshConfig->meshes().at(1));
  BOOST_TEST(mappingConfig.mappings().at(2).toMesh == meshConfig->meshes().at(0));
  BOOST_TEST(mappingConfig.mappings().at(2).direction == MappingConfiguration::WRITE);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

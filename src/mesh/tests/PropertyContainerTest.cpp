#include "mesh/PropertyContainer.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::mesh;

BOOST_AUTO_TEST_SUITE(MeshTests)
BOOST_AUTO_TEST_SUITE(PropertyContainerTest)

BOOST_AUTO_TEST_CASE(SinglePropertyContainer)
{
  PropertyContainer propertyContainer;
  int               propertyIndex = 0;
  int               integerValue  = 1;
  double            doubleValue   = 2.0;
  Eigen::Vector3d   vectorValue(0, 1, 2);

  BOOST_CHECK(not propertyContainer.hasProperty(propertyIndex));

  propertyContainer.setProperty(propertyIndex, integerValue);
  BOOST_CHECK(propertyContainer.hasProperty(propertyIndex));
  integerValue = 0;
  BOOST_TEST(propertyContainer.getProperty<int>(propertyIndex) == 1);

  propertyContainer.setProperty(propertyIndex, doubleValue);
  BOOST_TEST(propertyContainer.getProperty<double>(propertyIndex) == 2);

  propertyContainer.setProperty(propertyIndex, vectorValue);
  for (int dim = 0; dim < 3; dim++) {
    BOOST_TEST(
        propertyContainer.getProperty<Eigen::Vector3d>(propertyIndex)(dim) ==
        static_cast<double>(dim));
  }

  BOOST_CHECK(propertyContainer.deleteProperty(propertyIndex));
  BOOST_CHECK(not propertyContainer.hasProperty(propertyIndex));
}

BOOST_AUTO_TEST_CASE(HierarchicalPropertyContainers)
{
  PropertyContainer parent, child;
  int               index = 0;

  child.addParent(parent);
  BOOST_CHECK(not child.hasProperty(index));
  parent.setProperty(index, 1);
  std::vector<int> properties;
  child.getProperties(index, properties);
  BOOST_TEST(properties.size() == 1);
  BOOST_TEST(properties[0] == 1);
  child.setProperty(index, 2);
  BOOST_TEST(child.getProperty<int>(index) == 2);
  child.getProperties(index, properties);
  BOOST_TEST(properties.size() == 2);
  BOOST_TEST(properties[1] == 2);
}

BOOST_AUTO_TEST_CASE(MultipleParents)
{
  PropertyContainer parent0, parent1, child;

  child.addParent(parent0);
  child.addParent(parent1);
  int id = 0;
  BOOST_CHECK(!child.hasProperty(id));
  parent0.setProperty(id, 0);
  std::vector<int> properties;
  child.getProperties(id, properties);
  BOOST_TEST(properties.size() == 1);
  BOOST_TEST(properties[0] == 0);

  parent1.setProperty(id, 1);
  properties.clear();
  child.getProperties(id, properties);
  BOOST_TEST(properties.size() == 2);
  BOOST_TEST(properties[0] == 0);
  BOOST_TEST(properties[1] == 1);
}

BOOST_AUTO_TEST_SUITE_END() // PropertyContainerTEst
BOOST_AUTO_TEST_SUITE_END() // Mesh

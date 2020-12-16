#include <Eigen/Core>
#include <algorithm>
#include <list>
#include <memory>
#include <string>
#include <vector>
#include "action/Action.hpp"
#include "action/SharedPointer.hpp"
#include "action/SummationAction.hpp"
#include "action/config/ActionConfiguration.hpp"
#include "logging/Logger.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "xml/XMLTag.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(ActionTests)
BOOST_AUTO_TEST_SUITE(Summation)

BOOST_AUTO_TEST_CASE(SummationOneDimensional)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;
  PtrMesh          mesh(new Mesh("Mesh", 3, true, testing::nextMeshID()));
  int              dimension   = 1;
  PtrData          sourceData1 = mesh->createData("SourceData1", dimension);
  PtrData          sourceData2 = mesh->createData("SourceData2", dimension);
  PtrData          sourceData3 = mesh->createData("SourceData3", dimension);
  PtrData          targetData  = mesh->createData("TargetData", dimension);
  std::vector<int> sourceDataIDs{sourceData1->getID(), sourceData2->getID(), sourceData3->getID()};
  int              targetDataID = targetData->getID();
  mesh->createVertex(Eigen::Vector3d::Constant(0.0));
  mesh->createVertex(Eigen::Vector3d::Constant(1.0));
  mesh->createVertex(Eigen::Vector3d::Constant(2.0));

  mesh->allocateDataValues();
  auto &sourceValues1 = sourceData1->values();
  auto &sourceValues2 = sourceData2->values();
  auto &sourceValues3 = sourceData3->values();
  auto &targetValues  = targetData->values();
  sourceValues1 << 2.0, 3.0, 4.0;
  sourceValues2 << 1.0, 2.0, 3.0;
  sourceValues3 << 2.0, 3.0, 4.0;
  targetValues = Eigen::VectorXd::Zero(targetValues.size());

  action::SummationAction sum(
      action::SummationAction::WRITE_MAPPING_PRIOR, sourceDataIDs, targetDataID, mesh);

  sum.performAction(0.0, 0.25, 0.0, 0.25);
  BOOST_TEST(sourceValues1(0) == 2.0);
  BOOST_TEST(sourceValues1(1) == 3.0);
  BOOST_TEST(sourceValues1(2) == 4.0);
  BOOST_TEST(sourceValues2(0) == 1.0);
  BOOST_TEST(sourceValues2(1) == 2.0);
  BOOST_TEST(sourceValues2(2) == 3.0);
  BOOST_TEST(sourceValues3(0) == 2.0);
  BOOST_TEST(sourceValues3(1) == 3.0);
  BOOST_TEST(sourceValues3(2) == 4.0);
  BOOST_TEST(targetValues(0) == 5.0);
  BOOST_TEST(targetValues(1) == 8.0);
  BOOST_TEST(targetValues(2) == 11.0);

  sum.performAction(0.0, 0.25, 0.25, 0.25);
  BOOST_TEST(sourceValues1(0) == 2.0);
  BOOST_TEST(sourceValues1(1) == 3.0);
  BOOST_TEST(sourceValues1(2) == 4.0);
  BOOST_TEST(sourceValues2(0) == 1.0);
  BOOST_TEST(sourceValues2(1) == 2.0);
  BOOST_TEST(sourceValues2(2) == 3.0);
  BOOST_TEST(sourceValues3(0) == 2.0);
  BOOST_TEST(sourceValues3(1) == 3.0);
  BOOST_TEST(sourceValues3(2) == 4.0);
  BOOST_TEST(targetValues(0) == 5.0);
  BOOST_TEST(targetValues(1) == 8.0);
  BOOST_TEST(targetValues(2) == 11.0);
}

BOOST_AUTO_TEST_CASE(SummationThreeDimensional)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;
  int              dimension = 3;
  PtrMesh          mesh(new Mesh("Mesh", dimension, true, testing::nextMeshID()));
  PtrData          sourceData1 = mesh->createData("SourceData1", dimension);
  PtrData          sourceData2 = mesh->createData("SourceData2", dimension);
  PtrData          targetData  = mesh->createData("TargetData", dimension);
  std::vector<int> sourceDataIDs{sourceData1->getID(), sourceData2->getID()};
  int              targetDataID = targetData->getID();
  mesh->createVertex(Eigen::Vector3d::Constant(0.0));
  mesh->createVertex(Eigen::Vector3d::Constant(1.0));
  mesh->createVertex(Eigen::Vector3d::Constant(2.0));
  mesh->allocateDataValues();
  auto &sourceValues1 = sourceData1->values();
  auto &sourceValues2 = sourceData2->values();
  auto &targetValues  = targetData->values();
  sourceValues1 << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
  sourceValues2 << 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0;
  targetValues = Eigen::VectorXd::Zero(targetValues.size());

  action::SummationAction sum(
      action::SummationAction::WRITE_MAPPING_PRIOR, sourceDataIDs, targetDataID, mesh);

  sum.performAction(0.0, 0.25, 0.0, 0.25);
  BOOST_TEST(sourceValues1(0) == 1.0);
  BOOST_TEST(sourceValues2(0) == 2.0);
  BOOST_TEST(targetValues(0) == 3.0);

  BOOST_TEST(sourceValues1(1) == 2.0);
  BOOST_TEST(sourceValues2(1) == 3.0);
  BOOST_TEST(targetValues(1) == 5.0);

  BOOST_TEST(sourceValues1(2) == 3.0);
  BOOST_TEST(sourceValues2(2) == 4.0);
  BOOST_TEST(targetValues(2) == 7.0);

  BOOST_TEST(sourceValues1(3) == 4.0);
  BOOST_TEST(sourceValues2(3) == 5.0);
  BOOST_TEST(targetValues(3) == 9.0);

  BOOST_TEST(sourceValues1(4) == 5.0);
  BOOST_TEST(sourceValues2(4) == 6.0);
  BOOST_TEST(targetValues(4) == 11.0);

  BOOST_TEST(sourceValues1(5) == 6.0);
  BOOST_TEST(sourceValues2(5) == 7.0);
  BOOST_TEST(targetValues(5) == 13.0);

  BOOST_TEST(sourceValues1(6) == 7.0);
  BOOST_TEST(sourceValues2(6) == 8.0);
  BOOST_TEST(targetValues(6) == 15.0);

  BOOST_TEST(sourceValues1(7) == 8.0);
  BOOST_TEST(sourceValues2(7) == 9.0);
  BOOST_TEST(targetValues(7) == 17.0);

  BOOST_TEST(sourceValues1(8) == 9.0);
  BOOST_TEST(sourceValues2(8) == 10.0);
  BOOST_TEST(targetValues(8) == 19.0);
}

BOOST_AUTO_TEST_CASE(Configuration)
{
  PRECICE_TEST(1_rank);
  std::string                filename = testing::getPathToSources() + "/action/tests/SummationActionTest-testConfiguration-1.xml";
  xml::XMLTag                tag      = xml::getRootTag();
  mesh::PtrDataConfiguration dataConfig(new mesh::DataConfiguration(tag));
  dataConfig->setDimensions(3);
  mesh::PtrMeshConfiguration meshConfig(new mesh::MeshConfiguration(tag, dataConfig));
  meshConfig->setDimensions(3);
  action::ActionConfiguration config(tag, meshConfig);
  xml::configure(tag, xml::ConfigurationContext{}, filename);
  BOOST_TEST(config.actions().size() == 1);
  action::PtrAction action = config.actions().front();
  BOOST_TEST(action);
}

BOOST_AUTO_TEST_SUITE_END() // Summation
BOOST_AUTO_TEST_SUITE_END() // ActionTest

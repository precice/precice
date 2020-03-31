#include "action/SummationAction.hpp"
#include "action/config/ActionConfiguration.hpp"
#include "mesh/Data.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "testing/Testing.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(ActionTests)
BOOST_AUTO_TEST_SUITE(Summation, *testing::OnMaster())

BOOST_AUTO_TEST_CASE(SummationTwoSources)
{
  using namespace mesh;
  PtrMesh mesh(new Mesh("Mesh", 3, true, testing::nextMeshID()));
  PtrData sourceData1   = mesh->createData("SourceData1", 1);
  PtrData sourceData2	= mesh->createData("SourceData2", 1);
  PtrData targetData   = mesh->createData("TargetData", 1);
  std::vector<int> sourceDataID{sourceData1->getID(), sourceData2->getID()};	
  int     		   targetDataID = targetData->getID();
  mesh->createVertex(Eigen::Vector3d::Constant(0.0));
  mesh->createVertex(Eigen::Vector3d::Constant(1.0));
  mesh->createVertex(Eigen::Vector3d::Constant(2.0));

  mesh->allocateDataValues();
  auto &sourceValues1 = sourceData1->values();
  auto &sourceValues2 = sourceData2->values();
  auto &targetValues = targetData->values();
  sourceValues1 << 2.0, 3.0, 4.0;
  sourceValues2 << 4.0, 5.0, 1.0;
  targetValues = Eigen::VectorXd::Zero(targetValues.size());

  action::SummationAction sum(
      action::SummationAction::ALWAYS_PRIOR, sourceDataID, targetDataID, mesh);

  sum.performAction(0.0, 0.0, 0.0, 1.0);
  BOOST_TEST(sourceValues1[0] == 2.0);
  BOOST_TEST(sourceValues1[1] == 3.0);
  BOOST_TEST(sourceValues1[2] == 4.0);
  BOOST_TEST(sourceValues2[0] == 4.0);
  BOOST_TEST(sourceValues2[1] == 5.0);
  BOOST_TEST(sourceValues2[2] == 1.0);
  BOOST_TEST(targetValues[0] == 6.0);
  BOOST_TEST(targetValues[1] == 8.0);
  BOOST_TEST(targetValues[2] == 5.0);

  sum.performAction(0.0, 0.5, 0.5, 1.0);
  BOOST_TEST(sourceValues1[0] == 2.0);
  BOOST_TEST(sourceValues1[1] == 3.0);
  BOOST_TEST(sourceValues1[2] == 4.0);
  BOOST_TEST(sourceValues2[0] == 4.0);
  BOOST_TEST(sourceValues2[1] == 5.0);
  BOOST_TEST(sourceValues2[2] == 1.0);
  BOOST_TEST(targetValues[0] == 6.0);
  BOOST_TEST(targetValues[1] == 8.0);
  BOOST_TEST(targetValues[2] == 5.0);

  sum.performAction(0.0, 0.25, 0.75, 1.0);
  BOOST_TEST(sourceValues1[0] == 2.0);
  BOOST_TEST(sourceValues1[1] == 3.0);
  BOOST_TEST(sourceValues1[2] == 4.0);
  BOOST_TEST(sourceValues2[0] == 4.0);
  BOOST_TEST(sourceValues2[1] == 5.0);
  BOOST_TEST(sourceValues2[2] == 1.0);
  BOOST_TEST(targetValues[0] == 6.0);
  BOOST_TEST(targetValues[1] == 8.0);
  BOOST_TEST(targetValues[2] == 5.0);

  sum.performAction(0.0, 0.25, 1.0, 1.0);
  BOOST_TEST(sourceValues1[0] == 2.0);
  BOOST_TEST(sourceValues1[1] == 3.0);
  BOOST_TEST(sourceValues1[2] == 4.0);
  BOOST_TEST(sourceValues2[0] == 4.0);
  BOOST_TEST(sourceValues2[1] == 5.0);
  BOOST_TEST(sourceValues2[2] == 1.0);
  BOOST_TEST(targetValues[0] == 6.0);
  BOOST_TEST(targetValues[1] == 8.0);
  BOOST_TEST(targetValues[2] == 5.0);
}

BOOST_AUTO_TEST_CASE(SummationTripleSource)
{
  using namespace mesh;
  PtrMesh mesh(new Mesh("Mesh", 3, true, testing::nextMeshID()));
  PtrData sourceData1   = mesh->createData("SourceData1", 1);
  PtrData sourceData2   = mesh->createData("SourceData2", 1);
  PtrData sourceData3   = mesh->createData("SourceData3", 1);
  PtrData targetData   = mesh->createData("TargetData", 1);
  std::vector<int> sourceDataID{sourceData1->getID(), sourceData2->getID(), sourceData3->getID()};
  int     targetDataID = targetData->getID();
  mesh->createVertex(Eigen::Vector3d::Constant(0.0));
  mesh->createVertex(Eigen::Vector3d::Constant(1.0));
  mesh->createVertex(Eigen::Vector3d::Constant(2.0));

  mesh->allocateDataValues();
  auto &sourceValues1 = sourceData1->values();
  auto &sourceValues2 = sourceData2->values();
  auto &sourceValues3 = sourceData3->values();
  auto &targetValues = targetData->values();
  sourceValues1 << 2.0, 3.0, 4.0;
  sourceValues2 << 1.0, 2.0, 3.0;
  sourceValues3 << 2.0, 3.0, 4.0;
  targetValues = Eigen::VectorXd::Zero(targetValues.size());

  action::SummationAction sum(
      action::SummationAction::ALWAYS_PRIOR, sourceDataID, targetDataID, mesh);

  sum.performAction(0.0, 0.0, 0.0, 1.0);
  BOOST_TEST(sourceValues1[0] == 2.0);
  BOOST_TEST(sourceValues1[1] == 3.0);
  BOOST_TEST(sourceValues1[2] == 4.0);
  BOOST_TEST(sourceValues2[0] == 1.0);
  BOOST_TEST(sourceValues2[1] == 2.0);
  BOOST_TEST(sourceValues2[2] == 3.0);
  BOOST_TEST(sourceValues3[0] == 2.0);
  BOOST_TEST(sourceValues3[1] == 3.0);
  BOOST_TEST(sourceValues3[2] == 4.0);
  BOOST_TEST(targetValues[0] == 5.0);
  BOOST_TEST(targetValues[1] == 8.0);
  BOOST_TEST(targetValues[2] == 11.0);

  sum.performAction(0.0, 0.5, 0.5, 1.0);
  BOOST_TEST(sourceValues1[0] == 2.0);
  BOOST_TEST(sourceValues1[1] == 3.0);
  BOOST_TEST(sourceValues1[2] == 4.0);
  BOOST_TEST(sourceValues2[0] == 1.0);
  BOOST_TEST(sourceValues2[1] == 2.0);
  BOOST_TEST(sourceValues2[2] == 3.0);
  BOOST_TEST(sourceValues3[0] == 2.0);
  BOOST_TEST(sourceValues3[1] == 3.0);
  BOOST_TEST(sourceValues3[2] == 4.0);
  BOOST_TEST(targetValues[0] == 5.0);
  BOOST_TEST(targetValues[1] == 8.0);
  BOOST_TEST(targetValues[2] == 11.0);

  sum.performAction(0.0, 0.25, 0.75, 1.0);
  BOOST_TEST(sourceValues1[0] == 2.0);
  BOOST_TEST(sourceValues1[1] == 3.0);
  BOOST_TEST(sourceValues1[2] == 4.0);
  BOOST_TEST(sourceValues2[0] == 1.0);
  BOOST_TEST(sourceValues2[1] == 2.0);
  BOOST_TEST(sourceValues2[2] == 3.0);
  BOOST_TEST(sourceValues3[0] == 2.0);
  BOOST_TEST(sourceValues3[1] == 3.0);
  BOOST_TEST(sourceValues3[2] == 4.0);
  BOOST_TEST(targetValues[0] == 5.0);
  BOOST_TEST(targetValues[1] == 8.0);
  BOOST_TEST(targetValues[2] == 11.0);

  sum.performAction(0.0, 0.25, 1.0, 1.0);
  BOOST_TEST(sourceValues1[0] == 2.0);
  BOOST_TEST(sourceValues1[1] == 3.0);
  BOOST_TEST(sourceValues1[2] == 4.0);
  BOOST_TEST(sourceValues2[0] == 1.0);
  BOOST_TEST(sourceValues2[1] == 2.0);
  BOOST_TEST(sourceValues2[2] == 3.0);
  BOOST_TEST(sourceValues3[0] == 2.0);
  BOOST_TEST(sourceValues3[1] == 3.0);
  BOOST_TEST(sourceValues3[2] == 4.0);
  BOOST_TEST(targetValues[0] == 5.0);
  BOOST_TEST(targetValues[1] == 8.0);
  BOOST_TEST(targetValues[2] == 11.0);
}

BOOST_AUTO_TEST_CASE(Configuration)
{
  {
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
  {
    std::string                filename = testing::getPathToSources() + "/action/tests/SummationActionTest-testConfiguration-2.xml";
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
}

BOOST_AUTO_TEST_SUITE_END() // Summation
BOOST_AUTO_TEST_SUITE_END() // ActionTest

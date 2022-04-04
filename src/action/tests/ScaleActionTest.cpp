#include <Eigen/Core>
#include <algorithm>
#include <list>
#include <memory>
#include <string>
#include "action/Action.hpp"
#include "action/ScaleByAreaAction.hpp"
#include "action/ScaleByDtAction.hpp"
#include "action/SharedPointer.hpp"
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

namespace precice {
namespace mesh {
class Vertex;
} // namespace mesh
} // namespace precice

using namespace precice;

BOOST_AUTO_TEST_SUITE(ActionTests)
BOOST_AUTO_TEST_SUITE(Scale)

BOOST_AUTO_TEST_CASE(DivideByArea2D)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;
  PtrMesh mesh(new Mesh("Mesh", 2, testing::nextMeshID()));
  PtrData data   = mesh->createData("test-data", 1, 0_dataID);
  int     dataID = data->getID();
  Vertex &v0     = mesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  Vertex &v1     = mesh->createVertex(Eigen::Vector2d(1.0, 0.0));
  Vertex &v2     = mesh->createVertex(Eigen::Vector2d(1.0, 1.0));
  mesh->createEdge(v0, v1);
  mesh->createEdge(v1, v2);
  mesh->allocateDataValues();
  auto &values = data->values();
  values << 2.0, 3.0, 4.0;

  BOOST_TEST(values(0) == 2.0);
  BOOST_TEST(values(1) == 3.0);
  BOOST_TEST(values(2) == 4.0);

  // Scale properties
  action::ScaleByAreaAction scale(
      action::ScaleByAreaAction::WRITE_MAPPING_POST, dataID, mesh,
      action::ScaleByAreaAction::SCALING_DIVIDE_BY_AREA);

  scale.performAction(0.0, 0.0, 0.0, 0.0);

  BOOST_TEST(values(0) == 4.0);
  BOOST_TEST(values(1) == 3.0);
  BOOST_TEST(values(2) == 8.0);
}

BOOST_AUTO_TEST_CASE(DivideByArea3D)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;
  PtrMesh mesh(new Mesh("Mesh", 3, testing::nextMeshID()));
  PtrData data   = mesh->createData("test-data", 1, 0_dataID);
  int     dataID = data->getID();
  Vertex &v0     = mesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  Vertex &v1     = mesh->createVertex(Eigen::Vector3d(6.0, 2.0, 0.0));
  Vertex &v2     = mesh->createVertex(Eigen::Vector3d(0.0, 2.0, 0.0));
  Vertex &v3     = mesh->createVertex(Eigen::Vector3d(0.0, 0.0, 3.0));
  Vertex &v4     = mesh->createVertex(Eigen::Vector3d(2.0, 0.0, 3.0));
  Edge &  e0     = mesh->createEdge(v0, v1);
  Edge &  e1     = mesh->createEdge(v1, v2);
  Edge &  e2     = mesh->createEdge(v0, v2);
  Edge &  e3     = mesh->createEdge(v2, v3);
  Edge &  e4     = mesh->createEdge(v0, v3);
  Edge &  e5     = mesh->createEdge(v0, v4);
  Edge &  e6     = mesh->createEdge(v3, v4);
  mesh->createTriangle(e0, e1, e2);
  mesh->createTriangle(e2, e3, e4);
  mesh->createTriangle(e4, e5, e6);
  mesh->allocateDataValues();
  auto &values = data->values();
  values << 2.0, 3.0, 6.0, 5.0, 6.0;

  BOOST_TEST(values(0) == 2.0);
  BOOST_TEST(values(1) == 3.0);
  BOOST_TEST(values(2) == 6.0);
  BOOST_TEST(values(3) == 5.0);
  BOOST_TEST(values(4) == 6.0);

  // Scale properties
  action::ScaleByAreaAction scale(
      action::ScaleByAreaAction::WRITE_MAPPING_POST, dataID, mesh,
      action::ScaleByAreaAction::SCALING_DIVIDE_BY_AREA);

  scale.performAction(0.0, 0.0, 0.0, 0.0);

  BOOST_TEST(values(0) == 0.5);
  BOOST_TEST(values(1) == 1.5);
  BOOST_TEST(values(2) == 2.0);
  BOOST_TEST(values(3) == 2.5);
  BOOST_TEST(values(4) == 6.0);
}

BOOST_AUTO_TEST_CASE(ScaleByTimeStepSizeToTimeWindowSize)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;
  PtrMesh mesh(new Mesh("Mesh", 3, testing::nextMeshID()));
  PtrData sourceData   = mesh->createData("SourceData", 1, 0_dataID);
  PtrData targetData   = mesh->createData("TargetData", 1, 1_dataID);
  int     sourceDataID = sourceData->getID();
  int     targetDataID = targetData->getID();
  mesh->createVertex(Eigen::Vector3d::Constant(0.0));
  mesh->createVertex(Eigen::Vector3d::Constant(1.0));
  mesh->createVertex(Eigen::Vector3d::Constant(2.0));

  mesh->allocateDataValues();
  auto &sourceValues = sourceData->values();
  auto &targetValues = targetData->values();
  sourceValues << 2.0, 3.0, 4.0;
  targetValues = Eigen::VectorXd::Zero(targetValues.size());

  action::ScaleByDtAction scale(
      action::ScaleByDtAction::WRITE_MAPPING_POST, sourceDataID, targetDataID, mesh,
      action::ScaleByDtAction::SCALING_BY_TIME_STEP_TO_TIME_WINDOW_RATIO);

  scale.performAction(0.0, 0.0, 0.0, 1.0);
  BOOST_TEST(sourceValues(0) == 2.0);
  BOOST_TEST(sourceValues(1) == 3.0);
  BOOST_TEST(sourceValues(2) == 4.0);
  BOOST_TEST(targetValues(0) == 0.0);
  BOOST_TEST(targetValues(1) == 0.0);
  BOOST_TEST(targetValues(2) == 0.0);

  scale.performAction(0.0, 0.5, 0.5, 1.0);
  BOOST_TEST(sourceValues(0) == 2.0);
  BOOST_TEST(sourceValues(1) == 3.0);
  BOOST_TEST(sourceValues(2) == 4.0);
  BOOST_TEST(targetValues(0) == 1.0);
  BOOST_TEST(targetValues(1) == 1.5);
  BOOST_TEST(targetValues(2) == 2.0);

  scale.performAction(0.0, 0.25, 0.75, 1.0);
  BOOST_TEST(sourceValues(0) == 2.0);
  BOOST_TEST(sourceValues(1) == 3.0);
  BOOST_TEST(sourceValues(2) == 4.0);
  BOOST_TEST(targetValues(0) == 0.5);
  BOOST_TEST(targetValues(1) == 0.75);
  BOOST_TEST(targetValues(2) == 1.0);

  scale.performAction(0.0, 0.25, 1.0, 1.0);
  BOOST_TEST(sourceValues(0) == 2.0);
  BOOST_TEST(sourceValues(1) == 3.0);
  BOOST_TEST(sourceValues(2) == 4.0);
  BOOST_TEST(targetValues(0) == 0.5);
  BOOST_TEST(targetValues(1) == 0.75);
  BOOST_TEST(targetValues(2) == 1.0);
}

BOOST_AUTO_TEST_CASE(ScaleByComputedTimeWindowPart)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;
  PtrMesh mesh(new Mesh("Mesh", 3, testing::nextMeshID()));
  PtrData sourceData   = mesh->createData("SourceData", 1, 0_dataID);
  PtrData targetData   = mesh->createData("TargetData", 1, 1_dataID);
  int     sourceDataID = sourceData->getID();
  int     targetDataID = targetData->getID();
  mesh->createVertex(Eigen::Vector3d::Constant(0.0));
  mesh->createVertex(Eigen::Vector3d::Constant(1.0));
  mesh->createVertex(Eigen::Vector3d::Constant(2.0));
  mesh->allocateDataValues();
  auto &sourceValues = sourceData->values();
  auto &targetValues = targetData->values();
  sourceValues << 2.0, 3.0, 4.0;
  targetValues = Eigen::VectorXd::Zero(targetValues.size());

  action::ScaleByDtAction scale(
      action::ScaleByDtAction::WRITE_MAPPING_POST, sourceDataID, targetDataID, mesh,
      action::ScaleByDtAction::SCALING_BY_COMPUTED_TIME_WINDOW_PART_RATIO);

  scale.performAction(0.0, 0.0, 0.0, 1.0);
  BOOST_TEST(sourceValues(0) == 2.0);
  BOOST_TEST(sourceValues(1) == 3.0);
  BOOST_TEST(sourceValues(2) == 4.0);
  BOOST_TEST(targetValues(0) == 0.0);
  BOOST_TEST(targetValues(1) == 0.0);
  BOOST_TEST(targetValues(2) == 0.0);

  scale.performAction(0.0, 0.5, 0.5, 1.0);
  BOOST_TEST(sourceValues(0) == 2.0);
  BOOST_TEST(sourceValues(1) == 3.0);
  BOOST_TEST(sourceValues(2) == 4.0);
  BOOST_TEST(targetValues(0) == 1.0);
  BOOST_TEST(targetValues(1) == 1.5);
  BOOST_TEST(targetValues(2) == 2.0);

  scale.performAction(0.0, 0.5, 1.0, 1.0);
  BOOST_TEST(sourceValues(0) == 2.0);
  BOOST_TEST(sourceValues(1) == 3.0);
  BOOST_TEST(sourceValues(2) == 4.0);
  BOOST_TEST(targetValues(0) == 2.0);
  BOOST_TEST(targetValues(1) == 3.0);
  BOOST_TEST(targetValues(2) == 4.0);
}

BOOST_AUTO_TEST_CASE(Configuration)
{
  PRECICE_TEST(1_rank);
  {
    std::string                filename = testing::getPathToSources() + "/action/tests/ScaleActionTest-testConfiguration-1.xml";
    xml::XMLTag                tag      = xml::getRootTag();
    mesh::PtrDataConfiguration dataConfig(new mesh::DataConfiguration(tag));
    dataConfig->setDimensions(2);
    mesh::PtrMeshConfiguration meshConfig(new mesh::MeshConfiguration(tag, dataConfig));
    meshConfig->setDimensions(2);
    action::ActionConfiguration config(tag, meshConfig);
    xml::configure(tag, xml::ConfigurationContext{}, filename);
    BOOST_TEST(config.actions().size() == 1);
    auto &action = config.actions().front();
    BOOST_TEST(static_cast<bool>(action));
  }
  {
    std::string                filename = testing::getPathToSources() + "/action/tests/ScaleActionTest-testConfiguration-2.xml";
    xml::XMLTag                tag      = xml::getRootTag();
    mesh::PtrDataConfiguration dataConfig(new mesh::DataConfiguration(tag));
    dataConfig->setDimensions(2);
    mesh::PtrMeshConfiguration meshConfig(new mesh::MeshConfiguration(tag, dataConfig));
    meshConfig->setDimensions(2);
    action::ActionConfiguration config(tag, meshConfig);
    xml::configure(tag, xml::ConfigurationContext{}, filename);
    BOOST_TEST(config.actions().size() == 1);
    auto &action = config.actions().front();
    BOOST_TEST(static_cast<bool>(action));
  }
  {
    std::string                filename = testing::getPathToSources() + "/action/tests/ScaleActionTest-testConfiguration-3.xml";
    xml::XMLTag                tag      = xml::getRootTag();
    mesh::PtrDataConfiguration dataConfig(new mesh::DataConfiguration(tag));
    dataConfig->setDimensions(3);
    mesh::PtrMeshConfiguration meshConfig(new mesh::MeshConfiguration(tag, dataConfig));
    meshConfig->setDimensions(3);
    action::ActionConfiguration config(tag, meshConfig);
    xml::configure(tag, xml::ConfigurationContext{}, filename);
    BOOST_TEST(config.actions().size() == 1);
    auto &action = config.actions().front();
    BOOST_TEST(static_cast<bool>(action));
  }
  {
    std::string                filename = testing::getPathToSources() + "/action/tests/ScaleActionTest-testConfiguration-4.xml";
    xml::XMLTag                tag      = xml::getRootTag();
    mesh::PtrDataConfiguration dataConfig(new mesh::DataConfiguration(tag));
    dataConfig->setDimensions(3);
    mesh::PtrMeshConfiguration meshConfig(new mesh::MeshConfiguration(tag, dataConfig));
    meshConfig->setDimensions(3);
    action::ActionConfiguration config(tag, meshConfig);
    xml::configure(tag, xml::ConfigurationContext{}, filename);
    BOOST_TEST(config.actions().size() == 1);
    auto &action = config.actions().front();
    BOOST_TEST(static_cast<bool>(action));
  }
}

BOOST_AUTO_TEST_SUITE_END() // Scale
BOOST_AUTO_TEST_SUITE_END() // ActionTest

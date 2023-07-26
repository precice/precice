#include <Eigen/Core>
#include <algorithm>
#include <list>
#include <memory>
#include <string>
#include "action/Action.hpp"
#include "action/ScaleByAreaAction.hpp"
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

namespace precice::mesh {
class Vertex;
} // namespace precice::mesh

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
  Eigen::VectorXd v(3);
  v << 2.0, 3.0, 4.0;
  data->setSampleAtTime(1, time::Sample{1, v});

  const auto &values = data->values();
  BOOST_TEST(values(0) == 2.0);
  BOOST_TEST(values(1) == 3.0);
  BOOST_TEST(values(2) == 4.0);

  // Scale properties
  action::ScaleByAreaAction scale(
      action::ScaleByAreaAction::WRITE_MAPPING_POST, dataID, mesh,
      action::ScaleByAreaAction::SCALING_DIVIDE_BY_AREA);

  scale.performAction(0.0);

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
  Eigen::VectorXd v(5);
  v << 2.0, 3.0, 6.0, 5.0, 6.0;
  data->setSampleAtTime(1, time::Sample{1, v});

  const auto &values = data->values();
  BOOST_TEST(values(0) == 2.0);
  BOOST_TEST(values(1) == 3.0);
  BOOST_TEST(values(2) == 6.0);
  BOOST_TEST(values(3) == 5.0);
  BOOST_TEST(values(4) == 6.0);

  // Scale properties
  action::ScaleByAreaAction scale(
      action::ScaleByAreaAction::WRITE_MAPPING_POST, dataID, mesh,
      action::ScaleByAreaAction::SCALING_DIVIDE_BY_AREA);

  scale.performAction(0.0);

  BOOST_TEST(values(0) == 0.5);
  BOOST_TEST(values(1) == 1.5);
  BOOST_TEST(values(2) == 2.0);
  BOOST_TEST(values(3) == 2.5);
  BOOST_TEST(values(4) == 6.0);
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
}

BOOST_AUTO_TEST_SUITE_END() // Scale
BOOST_AUTO_TEST_SUITE_END() // ActionTest

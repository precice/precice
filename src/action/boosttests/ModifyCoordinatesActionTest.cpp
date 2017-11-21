#include "action/ModifyCoordinatesAction.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Vertex.hpp"
#include "testing/Testing.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(ActionTests)
BOOST_AUTO_TEST_SUITE(ModifyCoordinates)

BOOST_AUTO_TEST_CASE(AddToCoordinates)
{
  using namespace mesh;
  PtrMesh mesh(new Mesh("Mesh", 2, false));
  PtrData data   = mesh->createData("test-data", 2);
  int     dataID = data->getID();
  Vertex &v0     = mesh->createVertex(Eigen::Vector2d::Zero());
  Vertex &v1     = mesh->createVertex(Eigen::Vector2d::Constant(1.0));
  Edge &  edge   = mesh->createEdge(v0, v1);
  mesh->computeState();
  mesh->allocateDataValues();
  auto &values = data->values();

  // Create ApplyDisplacementsMeshAction
  action::ModifyCoordinatesAction modifyCoordinates(
      action::ModifyCoordinatesAction::ALWAYS_PRIOR, dataID, mesh,
      action::ModifyCoordinatesAction::ADD_TO_COORDINATES_MODE);

  // Validate coordinates of mesh before modifying it
  BOOST_TEST(testing::equals(v0.getCoords(), Eigen::Vector2d::Constant(0.0)));
  BOOST_TEST(testing::equals(v1.getCoords(), Eigen::Vector2d::Constant(1.0)));
  Eigen::Vector2d normalizedNormal(0.5, -0.5);
  normalizedNormal.normalize();
  BOOST_TEST(testing::equals(edge.getNormal(), normalizedNormal));

  // Set displacements
  values.segment(v0.getID() * 2, 2) = Eigen::VectorXd::Constant(2, 2.0);
  values.segment(v1.getID() * 2, 2) = Eigen::VectorXd::Constant(2, -2.0);

  // Apply displacements to  node coordinates
  modifyCoordinates.performAction(0.0, 0.0, 0.0, 0.0);

  // Validate coordinates of mesh after modifying it
  BOOST_TEST(testing::equals(v0.getCoords(), Eigen::Vector2d::Constant(2.0)));
  BOOST_TEST(testing::equals(v1.getCoords(), Eigen::Vector2d::Constant(-1.0)));
  normalizedNormal << -0.5, 0.5;
  normalizedNormal.normalize();
  BOOST_TEST(testing::equals(edge.getNormal(), normalizedNormal));
}

BOOST_AUTO_TEST_CASE(SubtractFromCoordinates)
{
  using namespace mesh;
  // Create geometryContext by faking a geometry but not using it to create
  // the mesh. The mesh is created by hand, such that references to the vertices
  // to be displaced are obtained.
  PtrMesh mesh(new Mesh("Mesh", 2, false));
  PtrData data   = mesh->createData("test-data", 2);
  int     dataID = data->getID();
  Vertex &v0     = mesh->createVertex(Eigen::Vector2d::Constant(0.0));
  Vertex &v1     = mesh->createVertex(Eigen::Vector2d::Constant(1.0));
  Edge &  edge   = mesh->createEdge(v0, v1);
  mesh->computeState();
  mesh->allocateDataValues();
  auto &values = data->values();

  // Create ApplyDisplacementsMeshAction
  action::ModifyCoordinatesAction modifyCoordinates(
      action::ModifyCoordinatesAction::ALWAYS_PRIOR, dataID, mesh,
      action::ModifyCoordinatesAction::SUBTRACT_FROM_COORDINATES_MODE);
  //   modifyCoordinates.loadMeshContext ( meshContext );

  // Validate coordinates of mesh before modifying it
  BOOST_TEST(testing::equals(v0.getCoords(), Eigen::Vector2d::Constant(0.0)));
  BOOST_TEST(testing::equals(v1.getCoords(), Eigen::Vector2d::Constant(1.0)));
  Eigen::Vector2d normalizedNormal(0.5, -0.5);
  normalizedNormal.normalize();
  BOOST_TEST(testing::equals(edge.getNormal(), normalizedNormal));

  // Set displacements
  values.segment(v0.getID() * 2, 2) = Eigen::VectorXd::Constant(2, -2.0);
  values.segment(v1.getID() * 2, 2) = Eigen::VectorXd::Constant(2, 2.0);

  // Apply displacements to  node coordinates
  modifyCoordinates.performAction(0.0, 0.0, 0.0, 0.0);

  // Validate coordinates of mesh after modifying it
  BOOST_TEST(testing::equals(v0.getCoords(), Eigen::Vector2d::Constant(2.0)));
  BOOST_TEST(testing::equals(v1.getCoords(), Eigen::Vector2d::Constant(-1.0)));
  normalizedNormal << -0.5, 0.5;
  normalizedNormal.normalize();
  BOOST_TEST(testing::equals(edge.getNormal(), normalizedNormal));
}

BOOST_AUTO_TEST_SUITE_END() // ModifyCoordinates
BOOST_AUTO_TEST_SUITE_END() // ActionTest

#include "../impl/WatchPoint.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Vertex.hpp"
#include "testing/Testing.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(PreciceTests)

BOOST_AUTO_TEST_CASE(WatchPoint)
{
  using namespace mesh;
  int dim = 2;
  using Eigen::VectorXd;
  // Setup geometry
  std::string name("rectangle");
  bool        flipNormals = false;
  PtrMesh     mesh(new Mesh(name, dim, flipNormals, testing::nextMeshID()));

  mesh::Vertex &v1 = mesh->createVertex(Eigen::Vector2d(1.0, 1.0));
  mesh::Vertex &v2 = mesh->createVertex(Eigen::Vector2d(2.0, 1.0));
  mesh::Vertex &v3 = mesh->createVertex(Eigen::Vector2d(1.0, 2.0));
  mesh::Vertex &v4 = mesh->createVertex(Eigen::Vector2d(2.0, 2.0));
  mesh->createEdge(v1, v2);
  mesh->createEdge(v1, v3);
  mesh->createEdge(v2, v4);
  mesh->createEdge(v3, v4);

  PtrData doubleData   = mesh->createData("DoubleData", 1);
  PtrData vectorData   = mesh->createData("VectorData", dim);
  auto &  doubleValues = doubleData->values();
  auto &  vectorValues = vectorData->values();
  mesh->computeState();
  mesh->allocateDataValues();
  doubleValues.setConstant(1.0);
  vectorValues.setConstant(1.0);

  // Create watchpoints
  Eigen::Vector2d  pointToWatch0(1.0, 1.0);
  std::string      filename0("precice-WatchPointTest-output0.log");
  impl::WatchPoint watchpoint0(pointToWatch0, mesh, filename0);
  Eigen::Vector2d  pointToWatch1(1.0, 1.5);
  std::string      filename1("precice-WatchPointTest-output1.log");
  impl::WatchPoint watchpoint1(pointToWatch1, mesh, filename1);

  // Initialize
  watchpoint0.initialize();
  watchpoint1.initialize();

  // Write output
  watchpoint0.exportPointData(0.0);
  watchpoint1.exportPointData(0.0);

  // Change data and write output again
  doubleValues.setConstant(3.0);
  vectorValues.setConstant(3.0);
  watchpoint0.exportPointData(1.0);
  watchpoint1.exportPointData(1.0);

  // Change Mesh and data
  mesh::Vertex& v5 = mesh->createVertex(Eigen::Vector2d(3.0, 2.0));
  mesh->createEdge(v4,v5);
  mesh->allocateDataValues();
  mesh->computeState();
  doubleValues.setConstant(5.0);
  vectorValues.setConstant(5.0);

  // Re-Initialize
  watchpoint0.initialize();
  watchpoint1.initialize();

  // Write output
  watchpoint0.exportPointData(2.0);
  watchpoint1.exportPointData(2.0);
}

BOOST_AUTO_TEST_SUITE_END() // Precice

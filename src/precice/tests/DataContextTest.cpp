#include <Eigen/Core>
#include <string>
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "precice/impl/DataContext.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::impl;

BOOST_AUTO_TEST_SUITE(PreciceTests)

BOOST_AUTO_TEST_SUITE(DataContextTests)

BOOST_AUTO_TEST_CASE(testDataContextWriteMapping)
{
  PRECICE_TEST(1_rank);

  // Create mesh object
  int           dimensions = 3;
  mesh::PtrMesh ptrMesh    = std::make_shared<mesh::Mesh>("MyMesh", dimensions, testing::nextMeshID());
  mesh::PtrData ptrData    = ptrMesh->createData("MyData", dimensions);

  auto &v1 = ptrMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  auto &v2 = ptrMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  auto &v3 = ptrMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 1.0));

  DataContext dataContext(ptrData, ptrMesh);
}

BOOST_AUTO_TEST_CASE(testDataContextReadMapping)
{
  PRECICE_TEST(1_rank);

  // Create mesh object
  int           dimensions = 3;
  mesh::PtrMesh ptrMesh    = std::make_shared<mesh::Mesh>("MyMesh", dimensions, testing::nextMeshID());
  mesh::PtrData ptrData    = ptrMesh->createData("MyData", dimensions);

  auto &v1 = ptrMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  auto &v2 = ptrMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  auto &v3 = ptrMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 1.0));

  DataContext dataContext(ptrData, ptrMesh);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

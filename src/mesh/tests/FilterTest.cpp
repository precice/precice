#include <Eigen/Core>
#include <iosfwd>
#include <string>
#include "logging/Logger.hpp"
#include "mesh/Filter.hpp"
#include "mesh/Mesh.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace Eigen;

BOOST_AUTO_TEST_SUITE(MeshTests)
BOOST_AUTO_TEST_SUITE(FilterTests)

BOOST_AUTO_TEST_CASE(Vertices2D)
{
  PRECICE_TEST(1_rank);

  mesh::Mesh dest("2D dest", 2, testing::nextMeshID());
  mesh::Mesh src("2D src", 2, testing::nextMeshID());

  auto &v0 = dest.createVertex(Vector2d::Constant(4.0)); // Add dummy data to check additivity
  auto &v1 = src.createVertex(Vector2d::Constant(3.0));
  src.createVertex(Vector2d::Constant(2.0));

  v1.tag();

  auto p = [](const mesh::Vertex &v) { return v.isTagged(); };

  mesh::filterMesh(dest, src, p);

  // dest should contain Constante(4) and Constant(3), but not Constant(2)
  BOOST_TEST(dest.nVertices() == 2);
  BOOST_TEST(dest.vertices()[0] == v0);
  BOOST_TEST(dest.vertices()[1] == v1);
}

BOOST_AUTO_TEST_CASE(Vertices3D)
{
  PRECICE_TEST(1_rank);

  mesh::Mesh dest("3D dest", 3, testing::nextMeshID());
  mesh::Mesh src("3D src", 3, testing::nextMeshID());

  auto &v0 = dest.createVertex(Vector3d::Constant(4.0)); // Add dummy data to check additivity
  auto &v1 = src.createVertex(Vector3d::Constant(3.0));
  src.createVertex(Vector3d::Constant(2.0));

  v1.tag();

  auto p = [](const mesh::Vertex &v) { return v.isTagged(); };

  mesh::filterMesh(dest, src, p);

  // dest should contain Constante(4) and Constant(3), but not Constant(2)
  BOOST_TEST(dest.nVertices() == 2);
  BOOST_TEST(dest.vertices()[0] == v0);
  BOOST_TEST(dest.vertices()[1] == v1);
}

BOOST_AUTO_TEST_CASE(Edges)
{
  PRECICE_TEST(1_rank);

  mesh::Mesh dest("3D dest", 3, testing::nextMeshID());
  mesh::Mesh src("3D src", 3, testing::nextMeshID());

  auto &v0 = dest.createVertex(Vector3d::Constant(4.0)); // Add dummy data to check additivity
  auto &v1 = src.createVertex(Vector3d::Constant(0.0));
  auto &v2 = src.createVertex(Vector3d::Constant(1.0));
  auto &v3 = src.createVertex(Vector3d::Constant(2.0));

  auto &e0 = src.createEdge(v1, v2);
  src.createEdge(v2, v3);

  v1.tag();
  v2.tag();

  auto p = [](const mesh::Vertex &v) { return v.isTagged(); };

  mesh::filterMesh(dest, src, p);

  BOOST_TEST(dest.nVertices() == 3);
  BOOST_TEST(dest.vertices()[0] == v0);
  BOOST_TEST(dest.vertices()[1] == v1);
  BOOST_TEST(dest.vertices()[2] == v2);

  // Only e0 should survive
  BOOST_TEST(dest.edges().size() == 1);
  BOOST_TEST(dest.edges()[0] == e0);
}

BOOST_AUTO_TEST_CASE(Triangles)
{
  PRECICE_TEST(1_rank);

  mesh::Mesh dest("3D dest", 3, testing::nextMeshID());
  mesh::Mesh src("3D src", 3, testing::nextMeshID());

  dest.createVertex(Vector3d::Constant(4.0)); // Add dummy data to check additivity
  auto &v1 = src.createVertex(Vector3d::Constant(0.0));
  auto &v2 = src.createVertex(Vector3d{1.0, 0.0, 0.0});
  auto &v3 = src.createVertex(Vector3d{0.0, 1.0, 0.0});
  auto &v4 = src.createVertex(Vector3d{0.0, 0.0, 1.0});

  auto &t1 = src.createTriangle(v1, v2, v3);
  src.createTriangle(v2, v3, v4);

  v1.tag();
  v2.tag();
  v3.tag();

  auto p = [](const mesh::Vertex &v) { return v.isTagged(); };

  mesh::filterMesh(dest, src, p);

  BOOST_TEST(dest.nVertices() == 4);

  // Only t1 should survive (because v4 not passed)
  BOOST_TEST(dest.triangles().size() == 1);
  BOOST_TEST(dest.triangles()[0] == t1);
}

BOOST_AUTO_TEST_CASE(Tetrahedra)
{
  PRECICE_TEST(1_rank);

  mesh::Mesh dest("3D dest", 3, testing::nextMeshID());
  mesh::Mesh src("3D src", 3, testing::nextMeshID());

  dest.createVertex(Vector3d::Constant(4.0)); // Add dummy data to check additivity
  auto &v1 = src.createVertex(Vector3d::Constant(0.0));
  auto &v2 = src.createVertex(Vector3d{1.0, 0.0, 0.0});
  auto &v3 = src.createVertex(Vector3d{0.0, 1.0, 0.0});
  auto &v4 = src.createVertex(Vector3d{0.0, 0.0, 1.0});
  auto &v5 = src.createVertex(Vector3d{0.0, 2.0, 1.0});

  auto &t1 = src.createTetrahedron(v1, v2, v3, v4);
  src.createTetrahedron(v2, v3, v4, v5);

  v1.tag();
  v2.tag();
  v3.tag();
  v4.tag();

  auto p = [](const mesh::Vertex &v) { return v.isTagged(); };

  mesh::filterMesh(dest, src, p);

  BOOST_TEST(dest.nVertices() == 5);

  // Only t1 should survive (because v5 not passed)
  BOOST_TEST(dest.tetrahedra().size() == 1);
  BOOST_TEST(dest.tetrahedra()[0] == t1);
}

BOOST_AUTO_TEST_SUITE_END() // Filter
BOOST_AUTO_TEST_SUITE_END() // Mesh

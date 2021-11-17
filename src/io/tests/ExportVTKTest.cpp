#include <Eigen/Core>
#include <algorithm>
#include <string>
#include "io/Export.hpp"
#include "io/ExportVTK.hpp"
#include "mesh/Mesh.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

namespace precice {
namespace mesh {
class Edge;
class Vertex;
} // namespace mesh
} // namespace precice

BOOST_AUTO_TEST_SUITE(IOTests)

BOOST_AUTO_TEST_SUITE(VTKExport)

using namespace precice;

BOOST_AUTO_TEST_CASE(ExportPolygonalMesh)
{
  PRECICE_TEST(1_rank);
  int           dim = 2;
  mesh::Mesh    mesh("MyMesh", dim, testing::nextMeshID());
  mesh::Vertex &v1 = mesh.createVertex(Eigen::Vector2d::Constant(0.0));
  mesh::Vertex &v2 = mesh.createVertex(Eigen::Vector2d::Constant(1.0));
  mesh::Vertex &v3 = mesh.createVertex(Eigen::Vector2d{1.0, 0.0});

  mesh.createEdge(v1, v2);
  mesh.createEdge(v2, v3);
  mesh.createEdge(v3, v1);

  io::ExportVTK exportVTK;
  std::string   filename = "io-VTKExport-ExportPolygonalMesh";
  std::string   location = "";
  exportVTK.doExport(filename, location, mesh);
}

BOOST_AUTO_TEST_CASE(ExportTriangulatedMesh)
{
  PRECICE_TEST(1_rank);
  int           dim = 3;
  mesh::Mesh    mesh("MyMesh", dim, testing::nextMeshID());
  mesh::Vertex &v1 = mesh.createVertex(Eigen::Vector3d::Constant(0.0));
  mesh::Vertex &v2 = mesh.createVertex(Eigen::Vector3d::Constant(1.0));
  mesh::Vertex &v3 = mesh.createVertex(Eigen::Vector3d{1.0, 0.0, 0.0});

  mesh::Edge &e1 = mesh.createEdge(v1, v2);
  mesh::Edge &e2 = mesh.createEdge(v2, v3);
  mesh::Edge &e3 = mesh.createEdge(v3, v1);
  mesh.createTriangle(e1, e2, e3);

  io::ExportVTK exportVTK;
  std::string   filename = "io-VTKExport-ExportTriangulatedMesh";
  std::string   location = "";
  exportVTK.doExport(filename, location, mesh);
}

BOOST_AUTO_TEST_SUITE_END() // ExportVTK
BOOST_AUTO_TEST_SUITE_END() // IOTests

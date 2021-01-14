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
  int             dim           = 2;
  bool            invertNormals = false;
  mesh::Mesh      mesh("MyMesh", dim, invertNormals, testing::nextMeshID());
  mesh::Vertex &  v1      = mesh.createVertex(Eigen::VectorXd::Constant(dim, 0.0));
  mesh::Vertex &  v2      = mesh.createVertex(Eigen::VectorXd::Constant(dim, 1.0));
  Eigen::VectorXd coords3 = Eigen::VectorXd::Constant(dim, 0.0);
  coords3(0)              = 1.0;
  mesh::Vertex &v3        = mesh.createVertex(coords3);

  mesh.createEdge(v1, v2);
  mesh.createEdge(v2, v3);
  mesh.createEdge(v3, v1);

  mesh.computeState();

  bool          exportNormals = true;
  io::ExportVTK exportVTK(exportNormals);
  std::string   filename = "io-VTKExport-ExportPolygonalMesh";
  std::string   location = "";
  exportVTK.doExport(filename, location, mesh);
}

BOOST_AUTO_TEST_CASE(ExportTriangulatedMesh)
{
  PRECICE_TEST(1_rank);
  int             dim           = 3;
  bool            invertNormals = false;
  mesh::Mesh      mesh("MyMesh", dim, invertNormals, testing::nextMeshID());
  mesh::Vertex &  v1      = mesh.createVertex(Eigen::VectorXd::Constant(dim, 0.0));
  mesh::Vertex &  v2      = mesh.createVertex(Eigen::VectorXd::Constant(dim, 1.0));
  Eigen::VectorXd coords3 = Eigen::VectorXd::Zero(dim);
  coords3(0)              = 1.0;
  mesh::Vertex &v3        = mesh.createVertex(coords3);

  mesh::Edge &e1 = mesh.createEdge(v1, v2);
  mesh::Edge &e2 = mesh.createEdge(v2, v3);
  mesh::Edge &e3 = mesh.createEdge(v3, v1);
  mesh.createTriangle(e1, e2, e3);
  mesh.computeState();

  bool          exportNormals = true;
  io::ExportVTK exportVTK(exportNormals);
  std::string   filename = "io-VTKExport-ExportTriangulatedMesh";
  std::string   location = "";
  exportVTK.doExport(filename, location, mesh);
}

BOOST_AUTO_TEST_SUITE_END() // ExportVTK
BOOST_AUTO_TEST_SUITE_END() // IOTests

#include "io/ExportVTK.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"
#include "testing/Testing.hpp"

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
  coords3[0]              = 1.0;
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
  coords3[0]              = 1.0;
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

BOOST_AUTO_TEST_CASE(ExportQuadMesh)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;
  int        dim           = 3;
  bool       invertNormals = false;
  mesh::Mesh mesh("QuadMesh", dim, invertNormals, testing::nextMeshID());
  // z=0 plane
  Vertex &v0 = mesh.createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  Vertex &v1 = mesh.createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  Vertex &v2 = mesh.createVertex(Eigen::Vector3d(1.0, 1.0, 0.0));
  Vertex &v3 = mesh.createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  // z=1 plane
  Vertex &v4 = mesh.createVertex(Eigen::Vector3d(0.0, 0.0, 1.0));
  Vertex &v5 = mesh.createVertex(Eigen::Vector3d(1.0, 0.0, 1.0));
  Vertex &v6 = mesh.createVertex(Eigen::Vector3d(1.0, 1.0, 1.0));
  Vertex &v7 = mesh.createVertex(Eigen::Vector3d(0.0, 1.0, 1.0));

  // z=0 plane
  Edge &e0 = mesh.createEdge(v0, v1);
  Edge &e1 = mesh.createEdge(v1, v2);
  Edge &e2 = mesh.createEdge(v2, v3);
  Edge &e3 = mesh.createEdge(v3, v0);
  // z=1 plane
  Edge &e4 = mesh.createEdge(v4, v5);
  Edge &e5 = mesh.createEdge(v5, v6);
  Edge &e6 = mesh.createEdge(v6, v7);
  Edge &e7 = mesh.createEdge(v7, v4);
  // inbetween edges
  Edge &e8  = mesh.createEdge(v0, v4);
  Edge &e9  = mesh.createEdge(v1, v5);
  Edge &e10 = mesh.createEdge(v2, v6);
  Edge &e11 = mesh.createEdge(v3, v7);

  // x-y plane
  mesh.createQuad(e3, e2, e1, e0);
  mesh.createQuad(e4, e5, e6, e7);
  // x-z plane
  mesh.createQuad(e0, e9, e4, e8);
  mesh.createQuad(e11, e6, e10, e2);
  // y-z plane
  mesh.createQuad(e8, e7, e11, e3);
  mesh.createQuad(e9, e1, e10, e5);

  mesh.computeState();

  bool          exportNormals = true;
  io::ExportVTK exportVTK(exportNormals);
  std::string   filename = "io-VTKExport-ExportQuadMesh";
  std::string   location = "";
  exportVTK.doExport(filename, location, mesh);
}

BOOST_AUTO_TEST_SUITE_END() // ExportVTK
BOOST_AUTO_TEST_SUITE_END() // IOTests

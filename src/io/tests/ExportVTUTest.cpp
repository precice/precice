#ifndef PRECICE_NO_MPI

#include <Eigen/Core>
#include <algorithm>
#include <map>
#include <string>
#include "com/SharedPointer.hpp"
#include "io/Export.hpp"
#include "io/ExportVTU.hpp"
#include "mesh/Mesh.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/Parallel.hpp"

namespace precice {
namespace mesh {
class Edge;
class Vertex;
} // namespace mesh
} // namespace precice

BOOST_AUTO_TEST_SUITE(IOTests)

using namespace precice;

BOOST_AUTO_TEST_SUITE(VTUExport)

BOOST_AUTO_TEST_CASE(ExportPolygonalMeshSerial)
{
  PRECICE_TEST(""_on(1_rank).setupMasterSlaves());
  int           dim = 2;
  mesh::Mesh    mesh("MyMesh", dim, testing::nextMeshID());
  mesh::Vertex &v1 = mesh.createVertex(Eigen::Vector2d::Zero());
  mesh::Vertex &v2 = mesh.createVertex(Eigen::Vector2d::Constant(1));
  mesh::Vertex &v3 = mesh.createVertex(Eigen::Vector2d{1.0, 0.0});

  mesh.createEdge(v1, v2);
  mesh.createEdge(v2, v3);
  mesh.createEdge(v3, v1);

  io::ExportVTU exportVTU;
  std::string   filename = "io-VTUExport-ExportPolygonalMesh";
  std::string   location = "";
  exportVTU.doExport(filename, location, mesh);
}

BOOST_AUTO_TEST_CASE(ExportPolygonalMesh)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  int        dim = 2;
  mesh::Mesh mesh("MyMesh", dim, testing::nextMeshID());

  if (context.isRank(0)) {
    mesh::Vertex &v1 = mesh.createVertex(Eigen::Vector2d::Zero());
    mesh::Vertex &v2 = mesh.createVertex(Eigen::Vector2d::Constant(1));
    mesh::Vertex &v3 = mesh.createVertex(Eigen::Vector2d{1.0, 0});

    mesh.createEdge(v1, v2);
    mesh.createEdge(v2, v3);
    mesh.createEdge(v3, v1);
    mesh.getVertexDistribution()[0] = {0, 1, 2};
    mesh.getVertexDistribution()[1] = {};
    mesh.getVertexDistribution()[2] = {3, 4, 5};
    mesh.getVertexDistribution()[3] = {6};
  } else if (context.isRank(1)) {
    // nothing
  } else if (context.isRank(2)) {
    mesh::Vertex &v1 = mesh.createVertex(Eigen::Vector2d::Constant(1));
    mesh::Vertex &v2 = mesh.createVertex(Eigen::Vector2d::Constant(2));
    mesh::Vertex &v3 = mesh.createVertex(Eigen::Vector2d{1.0, 0.0});

    mesh.createEdge(v1, v2);
    mesh.createEdge(v2, v3);
    mesh.createEdge(v3, v1);
  } else if (context.isRank(3)) {
    mesh.createVertex(Eigen::Vector2d::Constant(3.0));
  }

  io::ExportVTU exportVTU;
  std::string   filename = "io-ExportVTUTest-testExportPolygonalMesh";
  std::string   location = "";
  exportVTU.doExport(filename, location, mesh);
}

BOOST_AUTO_TEST_CASE(ExportTriangulatedMesh)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  int        dim = 3;
  mesh::Mesh mesh("MyMesh", dim, testing::nextMeshID());

  if (context.isRank(0)) {
    mesh::Vertex &v1 = mesh.createVertex(Eigen::Vector3d::Zero());
    mesh::Vertex &v2 = mesh.createVertex(Eigen::Vector3d::Constant(1));
    mesh::Vertex &v3 = mesh.createVertex(Eigen::Vector3d{1.0, 0.0, 0.0});

    mesh::Edge &e1 = mesh.createEdge(v1, v2);
    mesh::Edge &e2 = mesh.createEdge(v2, v3);
    mesh::Edge &e3 = mesh.createEdge(v3, v1);
    mesh.createTriangle(e1, e2, e3);

    mesh.getVertexDistribution()[0] = {0, 1, 2};
    mesh.getVertexDistribution()[1] = {};
    mesh.getVertexDistribution()[2] = {3, 4, 5};
    mesh.getVertexDistribution()[3] = {6};
  } else if (context.isRank(1)) {
    // nothing
  } else if (context.isRank(2)) {
    mesh::Vertex &v1 = mesh.createVertex(Eigen::Vector3d::Constant(1));
    mesh::Vertex &v2 = mesh.createVertex(Eigen::Vector3d::Constant(2));
    mesh::Vertex &v3 = mesh.createVertex(Eigen::Vector3d{0.0, 1.0, 0.0});

    mesh::Edge &e1 = mesh.createEdge(v1, v2);
    mesh::Edge &e2 = mesh.createEdge(v2, v3);
    mesh::Edge &e3 = mesh.createEdge(v3, v1);
    mesh.createTriangle(e1, e2, e3);
  } else if (context.isRank(3)) {
    mesh.createVertex(Eigen::Vector3d::Constant(3.0));
  }

  io::ExportVTU exportVTU;
  std::string   filename = "io-ExportVTUTest-testExportTriangulatedMesh";
  std::string   location = "";
  exportVTU.doExport(filename, location, mesh);
}

BOOST_AUTO_TEST_SUITE_END() // IOTests
BOOST_AUTO_TEST_SUITE_END() // VTUExport

#endif // PRECICE_NO_MPI

#ifndef PRECICE_NO_MPI

#include <Eigen/Core>
#include <algorithm>
#include <map>
#include <string>
#include "com/SharedPointer.hpp"
#include "io/Export.hpp"
#include "io/ExportVTKXML.hpp"
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

// void ExportVTKXMLTest:: run()
// {
//   PRECICE_TRACE();
//   typedef utils::Parallel Par;
//   if (Par::getCommunicatorSize() > 3){
//     const std::vector<int> ranksWanted = {0, 1, 2, 3};
//     MPI_Comm comm = Par::getRestrictedCommunicator(ranksWanted);
//     if (Par::getProcessRank() <= 3){
//       Par::setGlobalCommunicator(comm);
//       testMethod(testExportPolygonalMesh);
//       testMethod(testExportTriangulatedMesh);
//       testMethod(testExportQuadMesh);
//       Par::setGlobalCommunicator(Par::getCommunicatorWorld());
//     }
//   }
// }

BOOST_AUTO_TEST_SUITE(IOTests)

using namespace precice;

BOOST_AUTO_TEST_SUITE(VTKXMLExport)

BOOST_AUTO_TEST_CASE(ExportPolygonalMesh)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  int        dim           = 2;
  bool       invertNormals = false;
  mesh::Mesh mesh("MyMesh", dim, invertNormals, testing::nextMeshID());

  if (utils::Parallel::getProcessRank() == 0) {
    mesh::Vertex &  v1      = mesh.createVertex(Eigen::VectorXd::Zero(dim));
    mesh::Vertex &  v2      = mesh.createVertex(Eigen::VectorXd::Constant(dim, 1));
    Eigen::VectorXd coords3 = Eigen::VectorXd::Zero(dim);
    coords3[0]              = 1.0;
    mesh::Vertex &v3        = mesh.createVertex(coords3);

    mesh.createEdge(v1, v2);
    mesh.createEdge(v2, v3);
    mesh.createEdge(v3, v1);
    mesh.getVertexDistribution()[0] = {0, 1, 2};
    mesh.getVertexDistribution()[1] = {};
    mesh.getVertexDistribution()[2] = {3, 4, 5};
    mesh.getVertexDistribution()[3] = {6};
  } else if (utils::Parallel::getProcessRank() == 1) {

  } else if (utils::Parallel::getProcessRank() == 2) {
    mesh::Vertex &  v1      = mesh.createVertex(Eigen::VectorXd::Constant(dim, 1));
    mesh::Vertex &  v2      = mesh.createVertex(Eigen::VectorXd::Constant(dim, 2));
    Eigen::VectorXd coords3 = Eigen::VectorXd::Zero(dim);
    coords3[1]              = 1.0;
    mesh::Vertex &v3        = mesh.createVertex(coords3);

    mesh.createEdge(v1, v2);
    mesh.createEdge(v2, v3);
    mesh.createEdge(v3, v1);
  } else if (utils::Parallel::getProcessRank() == 3) {
    mesh.createVertex(Eigen::VectorXd::Constant(dim, 3.0));
  }

  mesh.computeState();

  bool             exportNormals = true;
  io::ExportVTKXML exportVTKXML(exportNormals);
  std::string      filename = "io-ExportVTKXMLTest-testExportPolygonalMesh";
  std::string      location = "";
  exportVTKXML.doExport(filename, location, mesh);
}

BOOST_AUTO_TEST_CASE(ExportTriangulatedMesh)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  int        dim           = 3;
  bool       invertNormals = false;
  mesh::Mesh mesh("MyMesh", dim, invertNormals, testing::nextMeshID());

  if (utils::Parallel::getProcessRank() == 0) {
    mesh::Vertex &  v1      = mesh.createVertex(Eigen::VectorXd::Zero(dim));
    mesh::Vertex &  v2      = mesh.createVertex(Eigen::VectorXd::Constant(dim, 1));
    Eigen::VectorXd coords3 = Eigen::VectorXd::Zero(dim);
    coords3[0]              = 1.0;
    mesh::Vertex &v3        = mesh.createVertex(coords3);

    mesh::Edge &e1 = mesh.createEdge(v1, v2);
    mesh::Edge &e2 = mesh.createEdge(v2, v3);
    mesh::Edge &e3 = mesh.createEdge(v3, v1);
    mesh.createTriangle(e1, e2, e3);

    mesh.getVertexDistribution()[0] = {0, 1, 2};
    mesh.getVertexDistribution()[1] = {};
    mesh.getVertexDistribution()[2] = {3, 4, 5};
    mesh.getVertexDistribution()[3] = {6};
  } else if (utils::Parallel::getProcessRank() == 1) {

  } else if (utils::Parallel::getProcessRank() == 2) {
    mesh::Vertex &  v1      = mesh.createVertex(Eigen::VectorXd::Constant(dim, 1));
    mesh::Vertex &  v2      = mesh.createVertex(Eigen::VectorXd::Constant(dim, 2));
    Eigen::VectorXd coords3 = Eigen::VectorXd::Zero(dim);
    coords3[1]              = 1.0;
    mesh::Vertex &v3        = mesh.createVertex(coords3);

    mesh::Edge &e1 = mesh.createEdge(v1, v2);
    mesh::Edge &e2 = mesh.createEdge(v2, v3);
    mesh::Edge &e3 = mesh.createEdge(v3, v1);
    mesh.createTriangle(e1, e2, e3);
  } else if (utils::Parallel::getProcessRank() == 3) {
    mesh.createVertex(Eigen::VectorXd::Constant(dim, 3.0));
  }

  mesh.computeState();

  bool             exportNormals = false;
  io::ExportVTKXML exportVTKXML(exportNormals);
  std::string      filename = "io-ExportVTKXMLTest-testExportTriangulatedMesh";
  std::string      location = "";
  exportVTKXML.doExport(filename, location, mesh);
}

BOOST_AUTO_TEST_CASE(ExportQuadMesh)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  using namespace mesh;
  int        dim           = 3;
  bool       invertNormals = false;
  mesh::Mesh mesh("QuadMesh", dim, invertNormals, testing::nextMeshID());

  if (utils::Parallel::getProcessRank() == 0) {
    mesh.getVertexDistribution()[0] = {};
    mesh.getVertexDistribution()[1] = {0, 1, 2, 3};
    mesh.getVertexDistribution()[2] = {4, 5, 6, 7, 8, 9};
    mesh.getVertexDistribution()[3] = {};
  } else if (utils::Parallel::getProcessRank() == 1) {
    // z=0 plane
    Vertex &v0 = mesh.createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
    Vertex &v1 = mesh.createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
    // z=1 plane
    Vertex &v4 = mesh.createVertex(Eigen::Vector3d(0.0, 0.0, 1.0));
    Vertex &v5 = mesh.createVertex(Eigen::Vector3d(1.0, 0.0, 1.0));
    // z=0 plane
    Edge &e0 = mesh.createEdge(v0, v1);
    // z=1 plane
    Edge &e4 = mesh.createEdge(v4, v5);
    // inbetween edges
    Edge &e8 = mesh.createEdge(v0, v4);
    Edge &e9 = mesh.createEdge(v1, v5);
    // x-z plane
    mesh.createQuad(e0, e9, e4, e8);
  } else if (utils::Parallel::getProcessRank() == 2) {
    // z=0 plane
    Vertex &v1 = mesh.createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
    Vertex &v2 = mesh.createVertex(Eigen::Vector3d(1.0, 1.0, 0.0));
    Vertex &v3 = mesh.createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
    // z=1 plane
    Vertex &v5 = mesh.createVertex(Eigen::Vector3d(1.0, 0.0, 1.0));
    Vertex &v6 = mesh.createVertex(Eigen::Vector3d(1.0, 1.0, 1.0));
    Vertex &v7 = mesh.createVertex(Eigen::Vector3d(0.0, 1.0, 1.0));
    // z=0 plane
    Edge &e1 = mesh.createEdge(v1, v2);
    Edge &e2 = mesh.createEdge(v2, v3);
    // z=1 plane
    Edge &e5 = mesh.createEdge(v5, v6);
    Edge &e6 = mesh.createEdge(v6, v7);
    // inbetween edges
    Edge &e9  = mesh.createEdge(v1, v5);
    Edge &e10 = mesh.createEdge(v2, v6);
    Edge &e11 = mesh.createEdge(v3, v7);
    // x-z plane
    mesh.createQuad(e11, e6, e10, e2);
    // y-z plane
    mesh.createQuad(e9, e1, e10, e5);
  } else if (utils::Parallel::getProcessRank() == 3) {
  }

  mesh.computeState();

  bool             exportNormals = false;
  io::ExportVTKXML exportVTKXML(exportNormals);
  std::string      filename = "io-ExportVTKXMLTest-testExportQuadMesh";
  std::string      location = "";
  exportVTKXML.doExport(filename, location, mesh);
}

BOOST_AUTO_TEST_SUITE_END() // IOTests
BOOST_AUTO_TEST_SUITE_END() // VTKXMLExport

#endif // PRECICE_NO_MPI

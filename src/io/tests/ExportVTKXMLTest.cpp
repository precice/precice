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
    coords3(0)              = 1.0;
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
    coords3(1)              = 1.0;
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
    coords3(0)              = 1.0;
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
    coords3(1)              = 1.0;
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

BOOST_AUTO_TEST_SUITE_END() // IOTests
BOOST_AUTO_TEST_SUITE_END() // VTKXMLExport

#endif // PRECICE_NO_MPI

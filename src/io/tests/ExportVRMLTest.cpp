#include <string>
#include "io/ExportVRML.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(IOTests)

using namespace precice;

BOOST_AUTO_TEST_CASE(ExportSimpleMesh, * testing::OnMaster())
{
  int dim = 2;
  mesh::Mesh::resetGeometryIDsGlobally();
  bool       flipNormals = false;
  mesh::Mesh mesh("test-cuboid", dim, flipNormals);

  mesh::Vertex &v1 = mesh.createVertex(Eigen::Vector2d(1.0, 1.0));
  mesh::Vertex &v2 = mesh.createVertex(Eigen::Vector2d(2.0, 1.0));
  mesh::Vertex &v3 = mesh.createVertex(Eigen::Vector2d(1.0, 2.0));
  mesh::Vertex &v4 = mesh.createVertex(Eigen::Vector2d(2.0, 2.0));
  mesh.createEdge(v1, v2);
  mesh.createEdge(v1, v3);
  mesh.createEdge(v2, v4);
  mesh.createEdge(v3, v4);

  mesh.createData("Data", dim);
  mesh.allocateDataValues();

  std::ostringstream stream;
  stream << "io-ExportVRMLTest-testExportCuboid-" << dim << "d.wrl";
  io::ExportVRML ex(false);
  std::string    location = "";
  ex.doExport(stream.str(), location, mesh);
}

BOOST_AUTO_TEST_SUITE_END() // IOTests

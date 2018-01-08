#include <string>
#include "io/ExportVRML.hpp"
#include "io/ExportVTK.hpp"
#include "io/ImportVRML.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Vertex.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(IOTests)

using namespace precice;

BOOST_AUTO_TEST_CASE(InternallyCreatedMesh, * testing::OnMaster())
{
  mesh::Mesh::resetGeometryIDsGlobally();
  bool       flipNormals = false;
  int        dim         = 2;
  mesh::Mesh mesh("TestCuboid", dim, flipNormals);
  // Set geoemtry sub-ids
  std::string nameSubID0("side-0");
  mesh.setSubID(nameSubID0);

  std::string nameSubID1 = "side-1";
  mesh.setSubID(nameSubID1);

  std::string nameSubID2 = "side-2";
  mesh.setSubID(nameSubID2);
  mesh::PtrData data = mesh.createData("TestData", 2);

  mesh::Vertex &v1 = mesh.createVertex(Eigen::Vector2d(1.0, 1.0));
  mesh::Vertex &v2 = mesh.createVertex(Eigen::Vector2d(2.0, 1.0));
  mesh::Vertex &v3 = mesh.createVertex(Eigen::Vector2d(1.0, 2.0));
  mesh::Vertex &v4 = mesh.createVertex(Eigen::Vector2d(2.0, 2.0));
  mesh.createEdge(v1, v2);
  mesh.createEdge(v1, v3);
  mesh.createEdge(v2, v4);
  mesh.createEdge(v3, v4);
  mesh.allocateDataValues();

  // Query mesh information
  size_t          vertexCount = mesh.vertices().size();
  size_t          edgeCount   = mesh.edges().size();
  Eigen::Vector2d coordsVertexOne(mesh.vertices()[0].getCoords());
  Eigen::Vector2d coordsVertexN(mesh.vertices()[vertexCount - 1].getCoords());
  Eigen::VectorXd dataOne   = Eigen::VectorXd::Constant(2, 2.0);
  Eigen::VectorXd dataN     = Eigen::VectorXd::Constant(2, 4.0);
  int             vertex0ID = mesh.vertices()[0].getID();
  int             vertexNID = mesh.vertices()[vertexCount - 1].getID();

  data->values().segment(vertex0ID * 2, 2) = dataOne;
  data->values().segment(vertexNID * 2, 2) = dataN;
  std::string filename("io-ExportandReimportVRMLTest-testInternallyCreatedMesh.wrl");

  // Export mesh
  bool           exportNormals = false;
  io::ExportVRML exportMesh(exportNormals);
  exportMesh.doExportCheckpoint(filename, mesh);

  // Reimport mesh
  mesh::Mesh    reimportedMesh("TestCuboid", 2, false);
  mesh::PtrData reimportedData = reimportedMesh.createData(data->getName(), data->getDimensions());
  reimportedMesh.setSubID(nameSubID0);
  reimportedMesh.setSubID(nameSubID1);
  reimportedMesh.setSubID(nameSubID2);

  std::string    location = "";
  io::ImportVRML importMesh(location);
  importMesh.doImportCheckpoint(filename, reimportedMesh, true);

  // Validate mesh information
  BOOST_TEST(reimportedMesh.vertices().size() == vertexCount);
  BOOST_TEST(reimportedMesh.edges().size() == edgeCount);
  BOOST_TEST(testing::equals(reimportedMesh.vertices()[0].getCoords(), coordsVertexOne));
  BOOST_TEST(testing::equals(reimportedMesh.vertices()[vertexCount - 1].getCoords(), coordsVertexN));
  auto &values   = reimportedData->values();
  int   vertexID = reimportedMesh.vertices()[0].getID();
  BOOST_TEST(values(vertexID), dataOne[0]);
  BOOST_TEST(values(vertexID), dataOne[1]);
  Eigen::Vector2d readDataOne;
  Eigen::Vector2d readDataN;
  for (size_t i = 0; i < 2; i++) {
    readDataOne[i] = data->values()[vertex0ID * 2 + i];
    readDataN[i]   = data->values()[vertexNID * 2 + i];
  }
  BOOST_TEST(math::equals(readDataOne, Eigen::VectorXd(dataOne)));
  BOOST_TEST(math::equals(readDataN, Eigen::VectorXd(dataN)));

  mesh::Mesh reimportedMesh1("TestCuboid", 2, false);
  importMesh.doImport(filename, reimportedMesh1);
  BOOST_TEST(reimportedMesh1.data().empty());
  BOOST_TEST(reimportedMesh1.propertyContainers().empty());
}

BOOST_AUTO_TEST_SUITE_END() // IOTests

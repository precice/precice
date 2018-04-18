#ifndef PRECICE_NO_SPIRIT2

#include <Eigen/Core>
#include "io/ExportVTK.hpp"
#include "io/ImportVRML.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Vertex.hpp"
#include "testing/Testing.hpp"
#include "utils/Globals.hpp"

BOOST_AUTO_TEST_SUITE(IOTests)

using namespace precice;
using namespace precice::io;

struct SetPathFixture {
  std::string testPath;

  SetPathFixture()
      : testPath(utils::getPathToSources() + "/io/tests/")
  {
  }
};

BOOST_FIXTURE_TEST_SUITE(VRMLImport, SetPathFixture)

BOOST_AUTO_TEST_CASE(ImportSquare)
{
  int           dim = 2;
  mesh::Mesh    mesh("MyMesh", dim, false);
  mesh::PtrData dataForces       = mesh.createData("Forces", dim);
  int           dataIDForces     = dataForces->getID();
  mesh::PtrData dataVelocities   = mesh.createData("Velocities", dim);
  int           dataIDVelocities = dataVelocities->getID();
  mesh.setSubID("test-subid");

  ImportVRML in(testPath);
  in.doImportCheckpoint("ImportVRMLTest2D.wrl", mesh, true);

  BOOST_TEST(mesh.vertices().size() == 4);
  BOOST_TEST(mesh.edges().size() == 4);

  // BOOST_TEST vertex coordinates
  BOOST_TEST(testing::equals(mesh.vertices()[0].getCoords(), Eigen::Vector2d(0.0, 0.0)));
  BOOST_TEST(testing::equals(mesh.vertices()[1].getCoords(), Eigen::Vector2d(1.0, 0.0)));
  BOOST_TEST(testing::equals(mesh.vertices()[2].getCoords(), Eigen::Vector2d(1.0, 1.0)));
  BOOST_TEST(testing::equals(mesh.vertices()[3].getCoords(), Eigen::Vector2d(0.0, 1.0)));

  // BOOST_TEST edge vertex coordinates
  BOOST_TEST(testing::equals(mesh.edges()[0].vertex(0).getCoords(), Eigen::Vector2d(0.0, 0.0)));
  BOOST_TEST(testing::equals(mesh.edges()[0].vertex(1).getCoords(), Eigen::Vector2d(1.0, 0.0)));
  BOOST_TEST(testing::equals(mesh.edges()[1].vertex(0).getCoords(), Eigen::Vector2d(1.0, 0.0)));
  BOOST_TEST(testing::equals(mesh.edges()[1].vertex(1).getCoords(), Eigen::Vector2d(1.0, 1.0)));
  BOOST_TEST(testing::equals(mesh.edges()[2].vertex(0).getCoords(), Eigen::Vector2d(1.0, 1.0)));
  BOOST_TEST(testing::equals(mesh.edges()[2].vertex(1).getCoords(), Eigen::Vector2d(0.0, 1.0)));
  BOOST_TEST(testing::equals(mesh.edges()[3].vertex(0).getCoords(), Eigen::Vector2d(0.0, 1.0)));
  BOOST_TEST(testing::equals(mesh.edges()[3].vertex(1).getCoords(), Eigen::Vector2d(0.0, 0.0)));

  // Validate data sets
  Eigen::VectorXd &forces = mesh.data(dataIDForces)->values();
  BOOST_TEST(forces(0) == 1.0);
  BOOST_TEST(forces(1) == 1.0);
  BOOST_TEST(forces(6) == 4.0);
  BOOST_TEST(forces(7) == 4.0);
  Eigen::VectorXd &velocities = mesh.data(dataIDVelocities)->values();
  BOOST_TEST(velocities(0) == 4.0);
  BOOST_TEST(velocities(1) == 4.0);
  BOOST_TEST(velocities(6) == 1.0);
  BOOST_TEST(velocities(7) == 1.0);

  // Compute mesh state and export to vtk file
  mesh.computeState();
  bool        exportNormals = true;
  ExportVTK   out(exportNormals);
  std::string location = "";
  out.doExport("io-VRMLImport-Square", location, mesh);
}

BOOST_AUTO_TEST_CASE(ImportCube)
{
  int           dim = 3;
  mesh::Mesh    mesh("MyMesh", dim, false);
  mesh::PtrData dataForces     = mesh.createData("Forces", dim);
  mesh::PtrData dataVelocities = mesh.createData("Velocities", dim);
  ImportVRML    in(testPath);
  in.doImportCheckpoint("ImportVRMLTest-Cube.wrl", mesh, true);
  mesh.computeState();
  bool        exportNormals = true;
  ExportVTK   out(exportNormals);
  std::string location = "";
  out.doExport("io-ImportVRML-Cube", location, mesh);
}

BOOST_AUTO_TEST_CASE(ImportSphere)
{
  int        dim = 3;
  mesh::Mesh mesh("MyMesh", dim, false);
  ImportVRML in(testPath);
  in.doImport("ImportVRMLTest-Sphere.wrl", mesh);
  mesh.computeState();
  bool        exportNormals = true;
  ExportVTK   out(exportNormals);
  std::string location = "";
  out.doExport("io-VRMLImport-Sphere", location, mesh);
}

BOOST_AUTO_TEST_CASE(ImportApe)
{
  int        dim = 3;
  mesh::Mesh mesh("MyMesh", dim, false);
  ImportVRML in(testPath);
  in.doImport("ImportVRMLTest-Ape.wrl", mesh);
  mesh.computeState();
  bool        exportNormals = true;
  ExportVTK   out(exportNormals);
  std::string location = "";
  out.doExport("io-VRMLImport-Ape", location, mesh);
}

BOOST_AUTO_TEST_CASE(ImportBunny)
{
  int        dim = 3;
  mesh::Mesh mesh("MyMesh", dim, false);
  ImportVRML in(testPath);
  in.doImport("ImportVRMLTest-Bunny.wrl", mesh);
  mesh.computeState();
  bool        exportNormals = true;
  ExportVTK   out(exportNormals);
  std::string location = "";
  out.doExport("io-VRMLImport-Bunny", location, mesh);
}

BOOST_AUTO_TEST_CASE(ImportDragon)
{
  int        dim = 3;
  mesh::Mesh mesh("MyMesh", dim, false);
  ImportVRML in(testPath);
  in.doImport("ImportVRMLTest-Dragon.wrl", mesh);
  mesh.computeState();
  bool        exportNormals = true;
  ExportVTK   out(exportNormals);
  std::string location = "";
  out.doExport("io-VRMLImport-Dragon", location, mesh);
}

BOOST_AUTO_TEST_CASE(ImportReactorPipe)
{
  int        dim = 3;
  mesh::Mesh mesh("MyMesh", dim, false);
  ImportVRML in(testPath);
  in.doImport("ImportVRMLTest-ReactorPipeCut.wrl", mesh);
  mesh.computeState();
  bool        exportNormals = true;
  ExportVTK   out(exportNormals);
  std::string location = "";
  out.doExport("io-VRMLImport-ReactorPipeCut", location, mesh);
}

BOOST_AUTO_TEST_SUITE_END() // ImportVRML
BOOST_AUTO_TEST_SUITE_END() // IOTests

#endif // not PRECICE_NO_SPIRIT2

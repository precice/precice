#include "ImportVRMLTest.hpp"
#include "io/ImportVRML.hpp"
#include "io/ExportVTK.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Data.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "utils/Globals.hpp"
#include <Eigen/Core>
#include "math/math.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::io::tests::ImportVRMLTest)

namespace precice {
namespace io {
namespace tests {

logging::Logger ImportVRMLTest:: _log ( "precice::io::ImportVRMLTest" );

ImportVRMLTest::ImportVRMLTest ()
:
   TestCase ( "io::ImportVRMLTest" ),
   _pathToTests ()
{}

void ImportVRMLTest:: setUp ()
{
  _pathToTests = utils::Globals::getPathToSources() + "/io/tests/";
}

void ImportVRMLTest:: run ()
{
# ifndef PRECICE_NO_SPIRIT2
  PRECICE_MASTER_ONLY {
    testMethod(testImportSquare);
    testMethod(testImportCube);
    testMethod(testImportSphere);
    testMethod(testImportApe);
    testMethod(testImportBunny);
    testMethod(testImportDragon);
    testMethod(testImportReactorPipe);
  }
# endif // not PRECICE_NO_SPIRIT2
}

void ImportVRMLTest:: testImportSquare()
{
  TRACE();
  int dim = 2;
  mesh::Mesh mesh("MyMesh", dim, false);
  mesh::PtrData dataForces = mesh.createData("Forces", dim);
  int dataIDForces = dataForces->getID();
  mesh::PtrData dataVelocities = mesh.createData("Velocities", dim);
  int dataIDVelocities = dataVelocities->getID();
  mesh.setSubID("test-subid");

  ImportVRML in(_pathToTests);
  in.doImportCheckpoint("ImportVRMLTest2D.wrl", mesh, true);

  validateEquals(mesh.vertices().size(), 4);
  validateEquals(mesh.edges().size(), 4);

  // Validate vertex coordinates
  validate(math::equals(mesh.vertices()[0].getCoords(),Eigen::Vector2d(0.0, 0.0)));
  validate(math::equals(mesh.vertices()[1].getCoords(),Eigen::Vector2d(1.0, 0.0)));
  validate(math::equals(mesh.vertices()[2].getCoords(),Eigen::Vector2d(1.0, 1.0)));
  validate(math::equals(mesh.vertices()[3].getCoords(),Eigen::Vector2d(0.0, 1.0)));

  // Validate edge vertex coordinates
  validate(math::equals(mesh.edges()[0].vertex(0).getCoords(), Eigen::Vector2d(0.0, 0.0)));
  validate(math::equals(mesh.edges()[0].vertex(1).getCoords(), Eigen::Vector2d(1.0, 0.0)));
  validate(math::equals(mesh.edges()[1].vertex(0).getCoords(), Eigen::Vector2d(1.0, 0.0)));
  validate(math::equals(mesh.edges()[1].vertex(1).getCoords(), Eigen::Vector2d(1.0, 1.0)));
  validate(math::equals(mesh.edges()[2].vertex(0).getCoords(), Eigen::Vector2d(1.0, 1.0)));
  validate(math::equals(mesh.edges()[2].vertex(1).getCoords(), Eigen::Vector2d(0.0, 1.0)));
  validate(math::equals(mesh.edges()[3].vertex(0).getCoords(), Eigen::Vector2d(0.0, 1.0)));
  validate(math::equals(mesh.edges()[3].vertex(1).getCoords(), Eigen::Vector2d(0.0, 0.0)));

  // Validate data sets
  Eigen::VectorXd& forces = mesh.data(dataIDForces)->values();
  validateNumericalEquals(forces(0), 1.0);
  validateNumericalEquals(forces(1), 1.0);
  validateNumericalEquals(forces(6), 4.0);
  validateNumericalEquals(forces(7), 4.0);
  Eigen::VectorXd& velocities = mesh.data(dataIDVelocities)->values();
  validateNumericalEquals(velocities(0), 4.0);
  validateNumericalEquals(velocities(1), 4.0);
  validateNumericalEquals(velocities(6), 1.0);
  validateNumericalEquals(velocities(7), 1.0);

  // Compute mesh state and export to vtk file
  mesh.computeState();
  bool exportNormals = true;
  ExportVTK out(exportNormals);
  std::string location = "";
  out.doExport("ImportVRMLTest2D",location, mesh);
}

void ImportVRMLTest:: testImportCube()
{
  TRACE();
  using namespace mesh;
  int dim = 3;
  Mesh mesh("MyMesh", dim, false);
  mesh::PtrData dataForces = mesh.createData("Forces", dim);
  mesh::PtrData dataVelocities = mesh.createData("Velocities", dim);
  ImportVRML in(_pathToTests);
  in.doImportCheckpoint("ImportVRMLTest-Cube.wrl", mesh, true);
  mesh.computeState();
  bool exportNormals = true;
  ExportVTK out(exportNormals);
  std::string location = "";
  out.doExport("ImportVRMLTest-Cube", location, mesh);
}

void ImportVRMLTest:: testImportSphere()
{
  TRACE();
  int dim = 3;
  mesh::Mesh mesh("MyMesh", dim, false);
  ImportVRML in(_pathToTests);
  in.doImport("ImportVRMLTest-Sphere.wrl", mesh);
  mesh.computeState();
  bool exportNormals = true;
  ExportVTK out(exportNormals);
  std::string location = "";
  out.doExport("ImportVRMLTest-Sphere",location, mesh);
}

void ImportVRMLTest:: testImportApe ()
{
  TRACE();
  int dim = 3;
  mesh::Mesh mesh("MyMesh", dim, false);
  ImportVRML in(_pathToTests);
  in.doImport("ImportVRMLTest-Ape.wrl", mesh);
  mesh.computeState();
  bool exportNormals = true;
  ExportVTK out(exportNormals);
  std::string location = "";
  out.doExport("ImportVRMLTest-Ape", location, mesh);
}

void ImportVRMLTest:: testImportBunny()
{
  TRACE();
  int dim = 3;
  mesh::Mesh mesh("MyMesh", dim, false);
  ImportVRML in(_pathToTests);
  in.doImport("ImportVRMLTest-Bunny.wrl", mesh);
  mesh.computeState();
  bool exportNormals = true;
  ExportVTK out(exportNormals);
  std::string location = "";
  out.doExport("ImportVRMLTest-Bunny", location, mesh);
}

void ImportVRMLTest:: testImportDragon()
{
  TRACE();
  int dim = 3;
  mesh::Mesh mesh("MyMesh", dim, false);
  ImportVRML in(_pathToTests);
  in.doImport("ImportVRMLTest-Dragon.wrl", mesh);
  mesh.computeState();
  bool exportNormals = true;
  ExportVTK out(exportNormals);
  std::string location = "";
  out.doExport("ImportVRMLTest-Dragon", location, mesh);
}

void ImportVRMLTest:: testImportReactorPipe()
{
  TRACE();
  int dim = 3;
  mesh::Mesh mesh("MyMesh", dim, false);
  ImportVRML in ( _pathToTests );
  in.doImport("ImportVRMLTest-ReactorPipeCut.wrl", mesh);
  mesh.computeState();
  bool exportNormals = true;
  ExportVTK out(exportNormals);
  std::string location = "";
  out.doExport("ImportVRMLTest-ReactorPipeCut", location, mesh);
}

}}} // namespace precice, io, tests


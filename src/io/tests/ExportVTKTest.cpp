#include <Eigen/Core>
#include <algorithm>
#include <filesystem>
#include <string>
#include "io/Export.hpp"
#include "io/ExportVTK.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

namespace precice::mesh {
class Edge;
class Vertex;
} // namespace precice::mesh

BOOST_AUTO_TEST_SUITE(IOTests)

BOOST_AUTO_TEST_SUITE(VTKExport)

using namespace precice;

PRECICE_TEST_SETUP(""_on(1_rank).setupIntraComm())
BOOST_AUTO_TEST_CASE(ExportScalar)
{
  PRECICE_TEST();
  mesh::Mesh mesh("Mesh", 2, testing::nextMeshID());
  mesh.createVertex(Eigen::Vector2d::Zero());
  mesh.createVertex(Eigen::Vector2d::Constant(1));
  mesh::PtrData data = mesh.createData("data", 1, 0_dataID);
  data->setSampleAtTime(0, time::Sample{1, 2}.setZero());

  io::ExportVTK exportCSV{"io-VTKExport", ".", mesh, io::Export::ExportKind::TimeWindows, 1, context.rank, context.size};
  exportCSV.doExport(0, 0.0);
  testing::expectFiles("Mesh-io-VTKExport.init.vtk");
}

PRECICE_TEST_SETUP(""_on(1_rank).setupIntraComm())
BOOST_AUTO_TEST_CASE(ExportVector)
{
  PRECICE_TEST();
  mesh::Mesh mesh("Mesh", 2, testing::nextMeshID());
  mesh.createVertex(Eigen::Vector2d::Zero());
  mesh.createVertex(Eigen::Vector2d::Constant(1));
  mesh::PtrData data = mesh.createData("data", 2, 0_dataID);
  data->setSampleAtTime(0, time::Sample{2, 2}.setZero());

  io::ExportVTK exportVTK{"io-VTKExport", ".", mesh, io::Export::ExportKind::TimeWindows, 1, context.rank, context.size};
  exportVTK.doExport(0, 0.0);
  testing::expectFiles("Mesh-io-VTKExport.init.vtk");
}

PRECICE_TEST_SETUP(""_on(1_rank).setupIntraComm())
BOOST_AUTO_TEST_CASE(ExportMissing)
{
  PRECICE_TEST();
  mesh::Mesh mesh("Mesh", 2, testing::nextMeshID());
  mesh.createVertex(Eigen::Vector2d::Zero());
  mesh.createVertex(Eigen::Vector2d::Constant(1));
  mesh::PtrData data = mesh.createData("data", 2, 0_dataID);
  // no sample
  io::ExportVTK exportVTK{"io-VTKExport", ".", mesh, io::Export::ExportKind::TimeWindows, 1, context.rank, context.size};
  exportVTK.doExport(0, 0.0);
  testing::expectFiles("Mesh-io-VTKExport.init.vtk");
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ExportDataWithGradient)
{
  PRECICE_TEST();
  int dimensions = 2;
  // Create mesh to map from
  mesh::Mesh    mesh("Mesh", dimensions, testing::nextMeshID());
  mesh::PtrData dataScalar = mesh.createData("dataScalar", 1, 0_dataID);
  mesh::PtrData dataVector = mesh.createData("dataVector", 2, 1_dataID);
  dataScalar->requireDataGradient();
  dataVector->requireDataGradient();
  mesh.createVertex(Eigen::Vector2d::Constant(0.0));
  mesh.createVertex(Eigen::Vector2d::Constant(1.0));

  time::Sample scalar(1, 2, dimensions);
  scalar.values.setLinSpaced(0, 1);
  scalar.gradients.setOnes();
  dataScalar->setSampleAtTime(0, scalar);

  time::Sample vectorial(dimensions, 2, dimensions);
  vectorial.values.setLinSpaced(0, 1);
  vectorial.gradients.setOnes();
  dataVector->setSampleAtTime(0, vectorial);

  io::ExportVTK exportVTK{"io-VTKExport", ".", mesh, io::Export::ExportKind::TimeWindows, 1, context.rank, context.size};
  exportVTK.doExport(0, 0.0);
  exportVTK.doExport(1, 1.0);
  testing::expectFiles("Mesh-io-VTKExport.init.vtk", "Mesh-io-VTKExport.dt1.vtk");
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ExportPolygonalMesh)
{
  PRECICE_TEST();
  int           dim = 2;
  mesh::Mesh    mesh("Mesh", dim, testing::nextMeshID());
  mesh::Vertex &v1 = mesh.createVertex(Eigen::Vector2d::Constant(0.0));
  mesh::Vertex &v2 = mesh.createVertex(Eigen::Vector2d::Constant(1.0));
  mesh::Vertex &v3 = mesh.createVertex(Eigen::Vector2d{1.0, 0.0});

  mesh.createEdge(v1, v2);
  mesh.createEdge(v2, v3);
  mesh.createEdge(v3, v1);

  io::ExportVTK exportVTK{"io-VTKExport", ".", mesh, io::Export::ExportKind::TimeWindows, 1, context.rank, context.size};
  exportVTK.doExport(0, 0.0);
  exportVTK.doExport(1, 1.0);
  testing::expectFiles("Mesh-io-VTKExport.init.vtk", "Mesh-io-VTKExport.dt1.vtk");
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ExportTriangulatedMesh)
{
  PRECICE_TEST();
  int           dim = 3;
  mesh::Mesh    mesh("Mesh", dim, testing::nextMeshID());
  mesh::Vertex &v1 = mesh.createVertex(Eigen::Vector3d::Constant(0.0));
  mesh::Vertex &v2 = mesh.createVertex(Eigen::Vector3d::Constant(1.0));
  mesh::Vertex &v3 = mesh.createVertex(Eigen::Vector3d{1.0, 0.0, 0.0});

  mesh::Edge &e1 = mesh.createEdge(v1, v2);
  mesh::Edge &e2 = mesh.createEdge(v2, v3);
  mesh::Edge &e3 = mesh.createEdge(v3, v1);
  mesh.createTriangle(e1, e2, e3);

  io::ExportVTK exportVTK{"io-VTKExport", ".", mesh, io::Export::ExportKind::TimeWindows, 1, context.rank, context.size};
  exportVTK.doExport(0, 0.0);
  exportVTK.doExport(1, 1.0);
  testing::expectFiles("Mesh-io-VTKExport.init.vtk", "Mesh-io-VTKExport.dt1.vtk");
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ExportTetrahedron)
{
  PRECICE_TEST();
  int           dim = 3;
  mesh::Mesh    mesh("Mesh", dim, testing::nextMeshID());
  mesh::Vertex &v1 = mesh.createVertex(Eigen::Vector3d::Constant(0.0));
  mesh::Vertex &v2 = mesh.createVertex(Eigen::Vector3d{1.0, 0.0, 0.0});
  mesh::Vertex &v3 = mesh.createVertex(Eigen::Vector3d{0.0, 1.0, 0.0});
  mesh::Vertex &v4 = mesh.createVertex(Eigen::Vector3d{0.0, 0.0, 1.0});

  mesh.createTetrahedron(v1, v2, v3, v4);

  io::ExportVTK exportVTK{"io-VTKExport", ".", mesh, io::Export::ExportKind::TimeWindows, 1, context.rank, context.size};
  exportVTK.doExport(0, 0.0);
  exportVTK.doExport(1, 1.0);
  testing::expectFiles("Mesh-io-VTKExport.init.vtk", "Mesh-io-VTKExport.dt1.vtk");
}

BOOST_AUTO_TEST_SUITE_END() // ExportVTK
BOOST_AUTO_TEST_SUITE_END() // IOTests

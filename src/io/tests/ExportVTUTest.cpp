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

namespace precice::mesh {
class Edge;
class Vertex;
} // namespace precice::mesh

BOOST_AUTO_TEST_SUITE(IOTests)

using namespace precice;

BOOST_AUTO_TEST_SUITE(VTUExport)

PRECICE_TEST_SETUP(""_on(1_rank).setupIntraComm())
BOOST_AUTO_TEST_CASE(ExportScalar)
{
  PRECICE_TEST();
  mesh::Mesh mesh("Mesh", 2, testing::nextMeshID());
  mesh.createVertex(Eigen::Vector2d::Zero());
  mesh.createVertex(Eigen::Vector2d::Constant(1));
  mesh::PtrData data = mesh.createData("data", 1, 0_dataID);
  data->setSampleAtTime(0, time::Sample{1, 2}.setZero());

  io::ExportVTU exportVTU{"io-VTUExport", ".", mesh, io::Export::ExportKind::TimeWindows, 1, context.rank, context.size};
  exportVTU.doExport(0, 0.0);
  testing::expectFiles("Mesh-io-VTUExport.init.vtu");
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

  io::ExportVTU exportVTU{"io-VTUExport", ".", mesh, io::Export::ExportKind::TimeWindows, 1, context.rank, context.size};
  exportVTU.doExport(0, 0.0);
  testing::expectFiles("Mesh-io-VTUExport.init.vtu");
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
  io::ExportVTU exportVTU{"io-VTUExport", ".", mesh, io::Export::ExportKind::TimeWindows, 1, context.rank, context.size};
  exportVTU.doExport(0, 0.0);
  testing::expectFiles("Mesh-io-VTUExport.init.vtu");
}

PRECICE_TEST_SETUP(""_on(2_ranks).setupIntraComm())
BOOST_AUTO_TEST_CASE(ExportScalarParallel)
{
  PRECICE_TEST();
  mesh::Mesh mesh("Mesh", 2, testing::nextMeshID());
  if (context.isPrimary()) {
    mesh.createVertex(Eigen::Vector2d::Zero());
  } else {
    mesh.createVertex(Eigen::Vector2d::Constant(1));
  }
  mesh::PtrData data = mesh.createData("data", 1, 0_dataID);
  data->setSampleAtTime(0, time::Sample{1, 1}.setZero());

  io::ExportVTU exportVTU{"io-VTUExport", ".", mesh, io::Export::ExportKind::TimeWindows, 1, context.rank, context.size};
  exportVTU.doExport(0, 0.0);

  testing::expectFiles(fmt::format("Mesh-io-VTUExport.init_{}.vtu", context.rank));
  if (context.isPrimary()) {
    testing::expectFiles("Mesh-io-VTUExport.init.pvtu");
  }
}

PRECICE_TEST_SETUP(""_on(2_ranks).setupIntraComm())
BOOST_AUTO_TEST_CASE(ExportVectorParallel)
{
  PRECICE_TEST();
  mesh::Mesh mesh("Mesh", 2, testing::nextMeshID());
  if (context.isPrimary()) {
    mesh.createVertex(Eigen::Vector2d::Zero());
  } else {
    mesh.createVertex(Eigen::Vector2d::Constant(1));
  }
  mesh::PtrData data = mesh.createData("data", 2, 0_dataID);
  data->setSampleAtTime(0, time::Sample{2, 1}.setZero());

  io::ExportVTU exportVTU{"io-VTUExport", ".", mesh, io::Export::ExportKind::TimeWindows, 1, context.rank, context.size};
  exportVTU.doExport(0, 0.0);

  testing::expectFiles(fmt::format("Mesh-io-VTUExport.init_{}.vtu", context.rank));
  if (context.isPrimary()) {
    testing::expectFiles("Mesh-io-VTUExport.init.pvtu");
  }
}

PRECICE_TEST_SETUP(""_on(2_ranks).setupIntraComm())
BOOST_AUTO_TEST_CASE(ExportMissingParallel)
{
  PRECICE_TEST();
  mesh::Mesh mesh("Mesh", 2, testing::nextMeshID());
  if (context.isPrimary()) {
    mesh.createVertex(Eigen::Vector2d::Zero());
  } else {
    mesh.createVertex(Eigen::Vector2d::Constant(1));
  }
  mesh::PtrData data = mesh.createData("data", 2, 0_dataID);
  BOOST_TEST_REQUIRE(data->timeStepsStorage().empty());
  // no sample
  io::ExportVTU exportVTU{"io-VTUExport", ".", mesh, io::Export::ExportKind::TimeWindows, 1, context.rank, context.size};
  exportVTU.doExport(0, 0.0);

  testing::expectFiles(fmt::format("Mesh-io-VTUExport.init_{}.vtu", context.rank));
  if (context.isPrimary()) {
    testing::expectFiles("Mesh-io-VTUExport.init.pvtu");
  }
}

PRECICE_TEST_SETUP(""_on(2_ranks).setupIntraComm())
BOOST_AUTO_TEST_CASE(ExportScalarAndMissingParallel)
{
  PRECICE_TEST();
  mesh::Mesh mesh("Mesh", 2, testing::nextMeshID());
  if (context.isPrimary()) {
    mesh.createVertex(Eigen::Vector2d::Zero());
  } else {
    mesh.createVertex(Eigen::Vector2d::Constant(1));
  }
  auto missing = mesh.createData("missing", 2, 0_dataID);
  BOOST_TEST_REQUIRE(missing->timeStepsStorage().empty());
  mesh::PtrData data = mesh.createData("data", 2, 1_dataID);
  data->setSampleAtTime(0, time::Sample{2, 1}.setZero());
  // no sample
  io::ExportVTU exportVTU{"io-VTUExport", ".", mesh, io::Export::ExportKind::TimeWindows, 1, context.rank, context.size};
  exportVTU.doExport(0, 0.0);

  testing::expectFiles(fmt::format("Mesh-io-VTUExport.init_{}.vtu", context.rank));
  if (context.isPrimary()) {
    testing::expectFiles("Mesh-io-VTUExport.init.pvtu");
  }
}

PRECICE_TEST_SETUP(""_on(1_rank).setupIntraComm())
BOOST_AUTO_TEST_CASE(ExportDataWithGradient2D)
{
  PRECICE_TEST();

  // Create mesh to map from
  int           dimensions = 2;
  mesh::Mesh    mesh("Mesh", dimensions, testing::nextMeshID());
  mesh::PtrData dataScalar = mesh.createData("dataScalar", 1, 0_dataID);
  mesh::PtrData dataVector = mesh.createData("dataVector", dimensions, 1_dataID);
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

  io::ExportVTU exportVTU{"io-VTUExport", ".", mesh, io::Export::ExportKind::TimeWindows, 1, context.rank, context.size};
  exportVTU.doExport(0, 0.0);
  exportVTU.doExport(1, 1.0);
  testing::expectFiles("Mesh-io-VTUExport.init.vtu", "Mesh-io-VTUExport.dt1.vtu");
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ExportDataWithGradient3D)
{
  PRECICE_TEST();
  int dimensions = 3;
  // Create mesh to map from
  mesh::Mesh    mesh("Mesh", dimensions, testing::nextMeshID());
  mesh::PtrData dataScalar = mesh.createData("dataScalar", 1, 0_dataID);
  mesh::PtrData dataVector = mesh.createData("dataVector", dimensions, 1_dataID);
  dataScalar->requireDataGradient();
  dataVector->requireDataGradient();

  mesh.createVertex(Eigen::Vector3d::Constant(0.0));
  mesh.createVertex(Eigen::Vector3d::Constant(1.0));

  time::Sample scalar(1, 2, dimensions);
  scalar.values.setLinSpaced(0, 1);
  scalar.gradients.setOnes();
  dataScalar->setSampleAtTime(0, scalar);

  time::Sample vectorial(dimensions, 2, dimensions);
  vectorial.values.setLinSpaced(0, 1);
  vectorial.gradients.setOnes();
  dataVector->setSampleAtTime(0, vectorial);

  io::ExportVTU exportVTU{"io-VTUExport", ".", mesh, io::Export::ExportKind::TimeWindows, 1, context.rank, context.size};
  exportVTU.doExport(0, 0.0);
  exportVTU.doExport(1, 1.0);
  testing::expectFiles("Mesh-io-VTUExport.init.vtu", "Mesh-io-VTUExport.dt1.vtu");
}

PRECICE_TEST_SETUP(""_on(1_rank).setupIntraComm())
BOOST_AUTO_TEST_CASE(ExportPolygonalMeshSerial)
{
  PRECICE_TEST();
  int           dim = 2;
  mesh::Mesh    mesh("Mesh", dim, testing::nextMeshID());
  mesh::Vertex &v1 = mesh.createVertex(Eigen::Vector2d::Zero());
  mesh::Vertex &v2 = mesh.createVertex(Eigen::Vector2d::Constant(1));
  mesh::Vertex &v3 = mesh.createVertex(Eigen::Vector2d{1.0, 0.0});

  mesh.createEdge(v1, v2);
  mesh.createEdge(v2, v3);
  mesh.createEdge(v3, v1);

  io::ExportVTU exportVTU{"io-VTUExport", ".", mesh, io::Export::ExportKind::TimeWindows, 1, context.rank, context.size};
  exportVTU.doExport(0, 0.0);
  exportVTU.doExport(1, 1.0);
  testing::expectFiles("Mesh-io-VTUExport.init.vtu", "Mesh-io-VTUExport.dt1.vtu");
}

PRECICE_TEST_SETUP(""_on(4_ranks).setupIntraComm())
BOOST_AUTO_TEST_CASE(ExportPolygonalMesh)
{
  PRECICE_TEST();
  int        dim = 2;
  mesh::Mesh mesh("Mesh", dim, testing::nextMeshID());

  if (context.isRank(0)) {
    mesh::Vertex &v1 = mesh.createVertex(Eigen::Vector2d::Zero());
    mesh::Vertex &v2 = mesh.createVertex(Eigen::Vector2d::Constant(1));
    mesh::Vertex &v3 = mesh.createVertex(Eigen::Vector2d{1.0, 0});

    mesh.createEdge(v1, v2);
    mesh.createEdge(v2, v3);
    mesh.createEdge(v3, v1);
    mesh.setVertexOffsets({3, 3, 6, 7});

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

  io::ExportVTU exportVTU{"io-VTUExport", ".", mesh, io::Export::ExportKind::TimeWindows, 1, context.rank, context.size};
  exportVTU.doExport(0, 0.0);
  exportVTU.doExport(1, 1.0);

  testing::expectFiles(fmt::format("Mesh-io-VTUExport.init_{}.vtu", context.rank), fmt::format("Mesh-io-VTUExport.dt1_{}.vtu", context.rank));
  if (context.isPrimary()) {
    testing::expectFiles("Mesh-io-VTUExport.init.pvtu", "Mesh-io-VTUExport.dt1.pvtu");
  }
}

PRECICE_TEST_SETUP(""_on(4_ranks).setupIntraComm())
BOOST_AUTO_TEST_CASE(ExportTriangulatedMesh)
{
  PRECICE_TEST();
  int        dim = 3;
  mesh::Mesh mesh("Mesh", dim, testing::nextMeshID());

  if (context.isRank(0)) {
    mesh::Vertex &v1 = mesh.createVertex(Eigen::Vector3d::Zero());
    mesh::Vertex &v2 = mesh.createVertex(Eigen::Vector3d::Constant(1));
    mesh::Vertex &v3 = mesh.createVertex(Eigen::Vector3d{1.0, 0.0, 0.0});

    mesh::Edge &e1 = mesh.createEdge(v1, v2);
    mesh::Edge &e2 = mesh.createEdge(v2, v3);
    mesh::Edge &e3 = mesh.createEdge(v3, v1);
    mesh.createTriangle(e1, e2, e3);

    mesh.setVertexOffsets({3, 3, 6, 7});

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

  io::ExportVTU exportVTU{"io-VTUExport", ".", mesh, io::Export::ExportKind::TimeWindows, 1, context.rank, context.size};
  exportVTU.doExport(0, 0.0);
  exportVTU.doExport(1, 1.0);

  testing::expectFiles(fmt::format("Mesh-io-VTUExport.init_{}.vtu", context.rank), fmt::format("Mesh-io-VTUExport.dt1_{}.vtu", context.rank));
  if (context.isPrimary()) {
    testing::expectFiles("Mesh-io-VTUExport.init.pvtu", "Mesh-io-VTUExport.dt1.pvtu");
  }
}

PRECICE_TEST_SETUP(""_on(4_ranks).setupIntraComm())
BOOST_AUTO_TEST_CASE(ExportSplitSquare)
{
  PRECICE_TEST();
  int        dim = 3;
  mesh::Mesh mesh("Mesh", dim, testing::nextMeshID());

  mesh::Vertex &vm = mesh.createVertex(Eigen::Vector3d::Zero());
  if (context.isRank(0)) {
    mesh::Vertex &v1  = mesh.createVertex(Eigen::Vector3d{-1.0, 1.0, 0.0});
    mesh::Vertex &v2  = mesh.createVertex(Eigen::Vector3d{1.0, 1.0, 0.0});
    mesh::Vertex &vo  = mesh.createVertex(Eigen::Vector3d{0.0, 2.0, 0.0});
    mesh::Edge   &em1 = mesh.createEdge(vm, v1);
    mesh::Edge   &e12 = mesh.createEdge(v1, v2);
    mesh::Edge   &e2m = mesh.createEdge(v2, vm);
    mesh.createTriangle(em1, e12, e2m);
    mesh::Edge &eo1 = mesh.createEdge(vo, v1);
    mesh::Edge &e2o = mesh.createEdge(v2, vo);
    mesh.createTriangle(eo1, e12, e2o);
    mesh.setVertexOffsets({3, 6, 9, 12});

  } else if (context.isRank(1)) {
    mesh::Vertex &v1  = mesh.createVertex(Eigen::Vector3d{1.0, -1.0, 0.0});
    mesh::Vertex &v2  = mesh.createVertex(Eigen::Vector3d{-1.0, -1.0, 0.0});
    mesh::Vertex &vo  = mesh.createVertex(Eigen::Vector3d{0.0, -2.0, 0.0});
    mesh::Edge   &em1 = mesh.createEdge(vm, v1);
    mesh::Edge   &e12 = mesh.createEdge(v1, v2);
    mesh::Edge   &e2m = mesh.createEdge(v2, vm);
    mesh.createTriangle(em1, e12, e2m);
    mesh::Edge &eo1 = mesh.createEdge(vo, v1);
    mesh::Edge &e2o = mesh.createEdge(v2, vo);
    mesh.createTriangle(eo1, e12, e2o);
  } else if (context.isRank(2)) {
    mesh::Vertex &v1  = mesh.createVertex(Eigen::Vector3d{-1.0, 1.0, 0.0});
    mesh::Vertex &v2  = mesh.createVertex(Eigen::Vector3d{-1.0, -1.0, 0.0});
    mesh::Vertex &vo  = mesh.createVertex(Eigen::Vector3d{-2.0, 0.0, 0.0});
    mesh::Edge   &em1 = mesh.createEdge(vm, v1);
    mesh::Edge   &e12 = mesh.createEdge(v1, v2);
    mesh::Edge   &e2m = mesh.createEdge(v2, vm);
    mesh.createTriangle(em1, e12, e2m);
    mesh::Edge &eo1 = mesh.createEdge(vo, v1);
    mesh::Edge &e2o = mesh.createEdge(v2, vo);
    mesh.createTriangle(eo1, e12, e2o);
  } else if (context.isRank(3)) {
    mesh::Vertex &v1  = mesh.createVertex(Eigen::Vector3d{1.0, 1.0, 0.0});
    mesh::Vertex &v2  = mesh.createVertex(Eigen::Vector3d{1.0, -1.0, 0.0});
    mesh::Vertex &vo  = mesh.createVertex(Eigen::Vector3d{2.0, 0.0, 0.0});
    mesh::Edge   &em1 = mesh.createEdge(vm, v1);
    mesh::Edge   &e12 = mesh.createEdge(v1, v2);
    mesh::Edge   &e2m = mesh.createEdge(v2, vm);
    mesh.createTriangle(em1, e12, e2m);
    mesh::Edge &eo1 = mesh.createEdge(vo, v1);
    mesh::Edge &e2o = mesh.createEdge(v2, vo);
    mesh.createTriangle(eo1, e12, e2o);
  }

  io::ExportVTU exportVTU{"io-VTUExport", ".", mesh, io::Export::ExportKind::TimeWindows, 1, context.rank, context.size};
  exportVTU.doExport(0, 0.0);
  exportVTU.doExport(1, 1.0);

  testing::expectFiles(fmt::format("Mesh-io-VTUExport.init_{}.vtu", context.rank), fmt::format("Mesh-io-VTUExport.dt1_{}.vtu", context.rank));
  if (context.isPrimary()) {
    testing::expectFiles("Mesh-io-VTUExport.init.pvtu", "Mesh-io-VTUExport.dt1.pvtu");
  }
}

PRECICE_TEST_SETUP(""_on(1_rank).setupIntraComm())
BOOST_AUTO_TEST_CASE(ExportOneTetrahedron)
{
  PRECICE_TEST();
  int           dim = 3;
  mesh::Mesh    mesh("Mesh", dim, testing::nextMeshID());
  mesh::Vertex &v0 = mesh.createVertex(Eigen::Vector3d::Zero());
  mesh::Vertex &v1 = mesh.createVertex(Eigen::Vector3d{1.0, 0.0, 0.0});
  mesh::Vertex &v2 = mesh.createVertex(Eigen::Vector3d{0.0, 1.0, 0.0});
  mesh::Vertex &v3 = mesh.createVertex(Eigen::Vector3d{0.0, 0.0, 1.0});

  mesh.createTetrahedron(v0, v1, v2, v3);

  io::ExportVTU exportVTU{"io-VTUExport", ".", mesh, io::Export::ExportKind::TimeWindows, 1, context.rank, context.size};
  exportVTU.doExport(0, 0.0);
  exportVTU.doExport(1, 1.0);
  testing::expectFiles("Mesh-io-VTUExport.init.vtu", "Mesh-io-VTUExport.dt1.vtu");
}

PRECICE_TEST_SETUP(""_on(4_ranks).setupIntraComm())
BOOST_AUTO_TEST_CASE(ExportPartitionedCube)
{
  PRECICE_TEST();
  // Unit cube is made of 6 tetrahedra. We have 3 ranks with 2 tetra each
  // as well as en empty rank. Empty rank is the 3rd
  int        dim = 3;
  mesh::Mesh mesh("Mesh", dim, testing::nextMeshID());

  if (context.isRank(0)) {
    mesh::Vertex &v000 = mesh.createVertex(Eigen::Vector3d{0.0, 0.0, 0.0});
    mesh::Vertex &v001 = mesh.createVertex(Eigen::Vector3d{0.0, 0.0, 1.0});
    mesh::Vertex &v011 = mesh.createVertex(Eigen::Vector3d{0.0, 1.0, 1.0});
    mesh::Vertex &v111 = mesh.createVertex(Eigen::Vector3d{1.0, 1.0, 1.0});
    mesh::Vertex &v010 = mesh.createVertex(Eigen::Vector3d{0.0, 1.0, 0.0});

    mesh.createTetrahedron(v000, v001, v011, v111);
    mesh.createTetrahedron(v000, v010, v011, v111);
    mesh.setVertexOffsets({4, 8, 8, 12});

  } else if (context.isRank(1)) {
    mesh::Vertex &v000 = mesh.createVertex(Eigen::Vector3d{0.0, 0.0, 0.0});
    mesh::Vertex &v001 = mesh.createVertex(Eigen::Vector3d{0.0, 0.0, 1.0});
    mesh::Vertex &v101 = mesh.createVertex(Eigen::Vector3d{1.0, 0.0, 1.0});
    mesh::Vertex &v111 = mesh.createVertex(Eigen::Vector3d{1.0, 1.0, 1.0});
    mesh::Vertex &v100 = mesh.createVertex(Eigen::Vector3d{1.0, 0.0, 0.0});

    mesh.createTetrahedron(v000, v001, v101, v111);
    mesh.createTetrahedron(v000, v100, v101, v111);
  } else if (context.isRank(3)) {
    mesh::Vertex &v000 = mesh.createVertex(Eigen::Vector3d{0.0, 0.0, 0.0});
    mesh::Vertex &v010 = mesh.createVertex(Eigen::Vector3d{0.0, 1.0, 0.0});
    mesh::Vertex &v100 = mesh.createVertex(Eigen::Vector3d{1.0, 0.0, 0.0});
    mesh::Vertex &v111 = mesh.createVertex(Eigen::Vector3d{1.0, 1.0, 1.0});
    mesh::Vertex &v110 = mesh.createVertex(Eigen::Vector3d{1.0, 1.0, 0.0});

    mesh.createTetrahedron(v000, v010, v110, v111);
    mesh.createTetrahedron(v000, v100, v110, v111);
  }

  io::ExportVTU exportVTU{"io-VTUExport", ".", mesh, io::Export::ExportKind::TimeWindows, 1, context.rank, context.size};
  exportVTU.doExport(0, 0.0);
  exportVTU.doExport(1, 1.0);

  testing::expectFiles(fmt::format("Mesh-io-VTUExport.init_{}.vtu", context.rank), fmt::format("Mesh-io-VTUExport.dt1_{}.vtu", context.rank));
  if (context.isPrimary()) {
    testing::expectFiles("Mesh-io-VTUExport.init.pvtu", "Mesh-io-VTUExport.dt1.pvtu");
  }
}

BOOST_AUTO_TEST_SUITE_END() // IOTests
BOOST_AUTO_TEST_SUITE_END() // VTUExport

#endif // PRECICE_NO_MPI

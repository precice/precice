#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include <string_view>
#include "logging/LogMacros.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/NearestNeighborMapping.hpp"
#include "math/constants.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Utils.hpp"
#include "mesh/Vertex.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "time/Sample.hpp"

using namespace precice;
using namespace precice::mesh;

BOOST_AUTO_TEST_SUITE(MappingTests)
BOOST_AUTO_TEST_SUITE(NearestNeighborMapping)

BOOST_AUTO_TEST_CASE(ConsistentNonIncremental)
{
  PRECICE_TEST(1_rank);
  int dimensions = 2;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  Vertex &inVertex0 = inMesh->createVertex(Eigen::Vector2d::Constant(0.0));
  Vertex &inVertex1 = inMesh->createVertex(Eigen::Vector2d::Constant(1.0));

  Eigen::VectorXd inValuesScalar = Eigen::VectorXd::Zero(2);
  Eigen::VectorXd inValuesVector = Eigen::VectorXd::Zero(4);
  inValuesScalar << 1.0, 2.0;
  inValuesVector << 1.0, 2.0, 3.0, 4.0;

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  Vertex &outVertex0 = outMesh->createVertex(Eigen::Vector2d::Constant(0.0));
  Vertex &outVertex1 = outMesh->createVertex(Eigen::Vector2d::Constant(1.0));

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::NearestNeighborMapping mapping(mapping::Mapping::CONSISTENT, dimensions);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Map data with coinciding vertices, has to result in equal values.
  mapping.computeMapping();
  Eigen::VectorXd outValuesScalar = Eigen::VectorXd::Zero(2);
  time::Sample    inSample(1, inValuesScalar);
  mapping.map(inSample, outValuesScalar);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValuesScalar(0) == inValuesScalar(0));
  BOOST_TEST(outValuesScalar(1) == inValuesScalar(1));
  Eigen::VectorXd outValuesVector = Eigen::VectorXd::Zero(4);
  inSample                        = time::Sample(2, inValuesVector);
  mapping.map(inSample, outValuesVector);
  BOOST_CHECK(equals(inValuesVector, outValuesVector));

  // Map data with almost coinciding vertices, has to result in equal values.
  inVertex0.setCoords(outVertex0.getCoords() + Eigen::Vector2d::Constant(0.1));
  inVertex1.setCoords(outVertex1.getCoords() + Eigen::Vector2d::Constant(0.1));
  mapping.computeMapping();
  inSample = time::Sample(1, inValuesScalar);
  mapping.map(inSample, outValuesScalar);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValuesScalar(0) == inValuesScalar(0));
  BOOST_TEST(outValuesScalar(1) == inValuesScalar(1));
  inSample = time::Sample(2, inValuesVector);
  mapping.map(inSample, outValuesVector);
  BOOST_CHECK(equals(inValuesVector, outValuesVector));

  // Map data with exchanged vertices, has to result in exchanged values.
  inVertex0.setCoords(outVertex1.getCoords());
  inVertex1.setCoords(outVertex0.getCoords());
  mapping.computeMapping();
  inSample = time::Sample(1, inValuesScalar);
  mapping.map(inSample, outValuesScalar);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValuesScalar(1) == inValuesScalar(0));
  BOOST_TEST(outValuesScalar(0) == inValuesScalar(1));
  inSample = time::Sample(2, inValuesVector);
  mapping.map(inSample, outValuesVector);
  Eigen::Vector4d expected(3.0, 4.0, 1.0, 2.0);
  BOOST_CHECK(equals(expected, outValuesVector));

  // Map data with coinciding output vertices, has to result in same values.
  outVertex1.setCoords(outVertex0.getCoords());
  mapping.computeMapping();
  inSample = time::Sample(1, inValuesScalar);
  mapping.map(inSample, outValuesScalar);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValuesScalar(1) == inValuesScalar(1));
  BOOST_TEST(outValuesScalar(0) == inValuesScalar(1));
  inSample = time::Sample(2, inValuesVector);
  mapping.map(inSample, outValuesVector);
  expected << 3.0, 4.0, 3.0, 4.0;
  BOOST_CHECK(equals(expected, outValuesVector));
}

BOOST_AUTO_TEST_CASE(ConservativeNonIncremental)
{
  PRECICE_TEST(1_rank);
  int dimensions = 2;

  // Create mesh to map from
  PtrMesh         inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  Vertex &        inVertex0 = inMesh->createVertex(Eigen::Vector2d::Constant(0.0));
  Vertex &        inVertex1 = inMesh->createVertex(Eigen::Vector2d::Constant(1.0));
  Eigen::VectorXd inValues  = Eigen::VectorXd::Zero(2);
  inValues(0)               = 1.0;
  inValues(1)               = 2.0;

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  Vertex &outVertex0 = outMesh->createVertex(Eigen::Vector2d::Constant(0.0));
  Vertex &outVertex1 = outMesh->createVertex(Eigen::Vector2d::Constant(1.0));

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::NearestNeighborMapping mapping(mapping::Mapping::CONSERVATIVE, dimensions);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Map data with coinciding vertices, has to result in equal values.
  mapping.computeMapping();
  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(2);
  time::Sample    inSample(1, inValues);
  mapping.map(inSample, outValues);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValues(0) == inValues(0));
  BOOST_TEST(outValues(1) == inValues(1));
  outValues = Eigen::VectorXd::Constant(outValues.size(), 0.0);

  // Map data with almost coinciding vertices, has to result in equal values.
  inVertex0.setCoords(outVertex0.getCoords() + Eigen::Vector2d::Constant(0.1));
  inVertex1.setCoords(outVertex1.getCoords() + Eigen::Vector2d::Constant(0.1));
  mapping.computeMapping();
  inSample = time::Sample(1, inValues);
  mapping.map(inSample, outValues);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValues(0) == inValues(0));
  BOOST_TEST(outValues(1) == inValues(1));
  outValues = Eigen::VectorXd::Constant(outValues.size(), 0.0);

  // Map data with exchanged vertices, has to result in exchanged values.
  inVertex0.setCoords(outVertex1.getCoords());
  inVertex1.setCoords(outVertex0.getCoords());
  mapping.computeMapping();
  inSample = time::Sample(1, inValues);
  mapping.map(inSample, outValues);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValues(1) == inValues(0));
  BOOST_TEST(outValues(0) == inValues(1));
  outValues = Eigen::VectorXd::Constant(outValues.size(), 0.0);

  // Map data with coinciding output vertices, has to result in double values.
  outVertex1.setCoords(Eigen::Vector2d::Constant(-1.0));
  mapping.computeMapping();
  inSample = time::Sample(1, inValues);
  mapping.map(inSample, outValues);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValues(0) == inValues(0) + inValues(1));
  BOOST_TEST(outValues(1) == 0.0);
}

BOOST_AUTO_TEST_CASE(ScaledConsistentNonIncremental)
{
  PRECICE_TEST(1_rank);
  int dimensions = 2;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  Vertex &inVertex0 = inMesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  Vertex &inVertex1 = inMesh->createVertex(Eigen::Vector2d{1.0, 0.0});
  Vertex &inVertex2 = inMesh->createVertex(Eigen::Vector2d{3.0, 0.0});
  Vertex &inVertex3 = inMesh->createVertex(Eigen::Vector2d{6.0, 0.0});

  inMesh->createEdge(inVertex0, inVertex1);
  inMesh->createEdge(inVertex1, inVertex2);
  inMesh->createEdge(inVertex2, inVertex3);

  Eigen::VectorXd inValues = Eigen::VectorXd::Zero(4);
  inValues(0)              = 1.0;
  inValues(1)              = 2.0;
  inValues(2)              = 3.0;
  inValues(3)              = 4.0;

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  Vertex &outVertex0 = outMesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  Vertex &outVertex1 = outMesh->createVertex(Eigen::Vector2d(0.8, 0.0));
  Vertex &outVertex2 = outMesh->createVertex(Eigen::Vector2d(3.0, 0.0));
  Vertex &outVertex3 = outMesh->createVertex(Eigen::Vector2d(6.2, 0.0));

  outMesh->createEdge(outVertex0, outVertex1);
  outMesh->createEdge(outVertex1, outVertex2);
  outMesh->createEdge(outVertex2, outVertex3);

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::NearestNeighborMapping mapping(mapping::Mapping::SCALED_CONSISTENT_SURFACE, dimensions);

  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  mapping.computeMapping();

  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(4);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  time::Sample inSample(1, inValues);
  mapping.map(inSample, outValues);

  auto inputIntegral  = mesh::integrateSurface(inMesh, inValues);
  auto outputIntegral = mesh::integrateSurface(outMesh, outValues);

  for (int dim = 0; dim < inputIntegral.size(); ++dim) {
    BOOST_TEST(inputIntegral(dim) == outputIntegral(dim));
  }

  double scaleFactor = outValues(0) / inValues(0);
  BOOST_TEST(scaleFactor != 1.0);

  BOOST_TEST(inValues(0) * scaleFactor == outValues(0));
  BOOST_TEST(inValues(1) * scaleFactor == outValues(1));
  BOOST_TEST(inValues(2) * scaleFactor == outValues(2));
  BOOST_TEST(inValues(3) * scaleFactor == outValues(3));
}

namespace {
template <int N>
PtrMesh create2DLinSpaceMesh(std::string_view name, int dims, double x0, double x1)
{
  static_assert(N > 1, "More than 1 vertex required");
  PtrMesh               mesh(new Mesh("InMesh", dims, testing::nextMeshID()));
  std::vector<Vertex *> vertices;
  double                dx = (x1 - x0) / (double) (N - 1);
  for (int i = 0; i < N; ++i) {
    double x = x0 + dx * i;
    vertices.push_back(&mesh->createVertex(Eigen::Vector2d{x, 0.0}));
  }
  auto first  = vertices.begin();
  auto second = ++vertices.begin();
  while (second != vertices.end()) {
    mesh->createEdge(**first, **second);
    ++first;
    ++second;
  }
  return mesh;
}
} // namespace

BOOST_AUTO_TEST_CASE(ScaledConsistentZeroData)
{
  PRECICE_TEST(1_rank);
  int dimensions = 2;

  // Create mesh to map from
  PtrMesh inMesh = create2DLinSpaceMesh<4>("InMesh", dimensions, 0, 1);

  // Create mesh to map to
  PtrMesh outMesh = create2DLinSpaceMesh<3>("OutMesh", dimensions, 0, 1);

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::NearestNeighborMapping mapping(mapping::Mapping::SCALED_CONSISTENT_SURFACE, dimensions);

  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);
  mapping.computeMapping();
  BOOST_TEST(mapping.hasComputedMapping() == true);

  Eigen::VectorXd inValues = Eigen::VectorXd::Zero(4);
  time::Sample    inSample(1, inValues);
  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(3);
  mapping.map(inSample, outValues);

  BOOST_TEST(!outValues.hasNaN());
  BOOST_TEST(outValues.isZero());

  auto inputIntegral = mesh::integrateSurface(inMesh, inValues);
  BOOST_TEST(!inputIntegral.hasNaN());
  BOOST_TEST(inputIntegral.isZero());

  auto outputIntegral = mesh::integrateSurface(outMesh, outValues);
  BOOST_TEST(!outputIntegral.hasNaN());
  BOOST_TEST(outputIntegral.isZero());
}

BOOST_AUTO_TEST_CASE(ScaledConsistentZeroIntegral)
{
  PRECICE_TEST(1_rank);
  int dimensions = 2;

  // Create mesh to map from
  PtrMesh inMesh = create2DLinSpaceMesh<4>("InMesh", dimensions, 0, 1);

  // Create mesh to map to
  PtrMesh outMesh = create2DLinSpaceMesh<3>("OutMesh", dimensions, 0, 1);

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::NearestNeighborMapping mapping(mapping::Mapping::SCALED_CONSISTENT_SURFACE, dimensions);

  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);
  mapping.computeMapping();
  BOOST_TEST(mapping.hasComputedMapping() == true);

  Eigen::Vector4d inValues{1, 1, -1, -1};
  time::Sample    inSample(1, inValues);
  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(3);
  mapping.map(inSample, outValues);

  BOOST_TEST(!outValues.hasNaN());
  BOOST_TEST((outValues.array() == 0.0).count() == 0);

  auto inputIntegral = mesh::integrateSurface(inMesh, inValues);
  BOOST_TEST(!inputIntegral.hasNaN());
  BOOST_TEST(inputIntegral(0) == 0.0);

  auto outputIntegral = mesh::integrateSurface(outMesh, outValues);
  BOOST_TEST(!outputIntegral.hasNaN());
  BOOST_TEST(outputIntegral(0) == 0.0);
}

BOOST_AUTO_TEST_CASE(ScaledConsistentZeroDataComponent)
{
  PRECICE_TEST(1_rank);
  int dimensions = 2;

  // Create mesh to map from
  PtrMesh inMesh = create2DLinSpaceMesh<4>("InMesh", dimensions, 0, 1);

  // Create mesh to map to
  PtrMesh outMesh = create2DLinSpaceMesh<3>("OutMesh", dimensions, 0, 1);

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::NearestNeighborMapping mapping(mapping::Mapping::SCALED_CONSISTENT_SURFACE, dimensions);

  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);
  mapping.computeMapping();
  BOOST_TEST(mapping.hasComputedMapping() == true);

  Eigen::VectorXd inValues(8);
  inValues << 1, 0, 0.5, 0, 1.5, 0, 1, 0;
  time::Sample    inSample(2, inValues);
  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(6);
  mapping.map(inSample, outValues);

  BOOST_TEST(!outValues.hasNaN());
  BOOST_TEST((outValues.array() == 0.0).count() == 3);

  auto inputIntegral = mesh::integrateSurface(inMesh, inValues);
  BOOST_TEST(!inputIntegral.hasNaN());
  BOOST_TEST(inputIntegral(0) == 1.0);
  BOOST_TEST(inputIntegral(1) == 0.0);

  auto outputIntegral = mesh::integrateSurface(outMesh, outValues);
  BOOST_TEST(!outputIntegral.hasNaN());
  BOOST_TEST(outputIntegral(0) == 1.0);
  BOOST_TEST(outputIntegral(1) == 0.0);
}

BOOST_AUTO_TEST_CASE(ScaledConsistentZeroIntegralComponent)
{
  PRECICE_TEST(1_rank);
  int dimensions = 2;

  // Create mesh to map from
  PtrMesh inMesh = create2DLinSpaceMesh<4>("InMesh", dimensions, 0, 1);

  // Create mesh to map to
  PtrMesh outMesh = create2DLinSpaceMesh<3>("OutMesh", dimensions, 0, 1);

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::NearestNeighborMapping mapping(mapping::Mapping::SCALED_CONSISTENT_SURFACE, dimensions);

  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);
  mapping.computeMapping();
  BOOST_TEST(mapping.hasComputedMapping() == true);

  Eigen::VectorXd inValues(8);
  inValues << 2, 1, 3, 1, 1, -1, 2, -1;
  time::Sample    inSample(2, inValues);
  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(6);
  mapping.map(inSample, outValues);

  BOOST_TEST(!outValues.hasNaN());
  BOOST_TEST((outValues.array() == 0.0).count() == 0);

  auto inputIntegral = mesh::integrateSurface(inMesh, inValues);
  BOOST_TEST(!inputIntegral.hasNaN());
  BOOST_TEST(inputIntegral(0) == 2.0);
  BOOST_TEST(inputIntegral(1) == 0.0);

  auto outputIntegral = mesh::integrateSurface(outMesh, outValues);
  BOOST_TEST(!outputIntegral.hasNaN());
  BOOST_TEST(outputIntegral(0) == 2.0);
  BOOST_TEST(outputIntegral(1) == 0.0);
}

BOOST_AUTO_TEST_CASE(ScaledConsistentVolume2D)
{
  PRECICE_TEST(1_rank);
  int dimensions = 2;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));

  // One square with 3 triangles
  Vertex &inVertex0 = inMesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  Vertex &inVertex1 = inMesh->createVertex(Eigen::Vector2d{1.0, 0.0});
  Vertex &inVertex2 = inMesh->createVertex(Eigen::Vector2d{1.0, 1.0});
  Vertex &inVertex3 = inMesh->createVertex(Eigen::Vector2d{0.0, 1.0});
  Vertex &inVertex4 = inMesh->createVertex(Eigen::Vector2d{0.5, 1.0});

  inMesh->createTriangle(inVertex0, inVertex1, inVertex4);
  inMesh->createTriangle(inVertex0, inVertex3, inVertex4);
  inMesh->createTriangle(inVertex1, inVertex2, inVertex4);

  Eigen::VectorXd inValues(5);
  inValues(0) = 1.0;
  inValues(1) = 2.0;
  inValues(2) = 3.0;
  inValues(3) = 4.0;
  inValues(4) = 5.0;

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));

  // Unit square as 2 triangles
  Vertex &outVertex0 = outMesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  Vertex &outVertex1 = outMesh->createVertex(Eigen::Vector2d(1.0, 0.0));
  Vertex &outVertex2 = outMesh->createVertex(Eigen::Vector2d(1.0, 1.0));
  Vertex &outVertex3 = outMesh->createVertex(Eigen::Vector2d(0.0, 1.0));

  outMesh->createTriangle(outVertex0, outVertex1, outVertex2);
  outMesh->createTriangle(outVertex0, outVertex2, outVertex3);

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::NearestNeighborMapping mapping(mapping::Mapping::SCALED_CONSISTENT_VOLUME, dimensions);

  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  Eigen::VectorXd outValues(5);
  mapping.computeMapping();
  time::Sample inSample(1, inValues);
  mapping.map(inSample, outValues);

  BOOST_TEST(mapping.hasComputedMapping() == true);

  auto inputIntegral  = mesh::integrateVolume(inMesh, inValues);
  auto outputIntegral = mesh::integrateVolume(outMesh, outValues);

  Eigen::VectorXd expectedIntegral(1);
  expectedIntegral << 3.0;

  BOOST_TEST(inputIntegral(0) == expectedIntegral(0));

  for (int dim = 0; dim < inputIntegral.size(); ++dim) {
    BOOST_TEST(inputIntegral(dim) == outputIntegral(dim));
  }

  double scaleFactor = outValues(0) / inValues(0);
  BOOST_TEST(scaleFactor != 1.0);

  BOOST_TEST(math::equals(inValues(0) * scaleFactor, outValues(0)));
  BOOST_TEST(inValues(1) * scaleFactor == outValues(1));
  BOOST_TEST(inValues(2) * scaleFactor == outValues(2));
  BOOST_TEST(inValues(3) * scaleFactor == outValues(3));
}

BOOST_AUTO_TEST_CASE(ScaledConsistentVolume3D)
{
  PRECICE_TEST(1_rank);
  int dimensions = 3;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));

  // One tetra on "out". The "in" has the same but split into two
  Vertex &inVertex0 = inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  Vertex &inVertex1 = inMesh->createVertex(Eigen::Vector3d{0.5, 1.0, 0.0});
  Vertex &inVertex2 = inMesh->createVertex(Eigen::Vector3d{1.0, 0.0, 0.0});
  Vertex &inVertex3 = inMesh->createVertex(Eigen::Vector3d{0.5, 0.0, 1.0});
  Vertex &inVertex4 = inMesh->createVertex(Eigen::Vector3d{0.5, 0.0, 0.0});

  inMesh->createTetrahedron(inVertex0, inVertex1, inVertex3, inVertex4);
  inMesh->createTetrahedron(inVertex1, inVertex2, inVertex3, inVertex4);

  Eigen::VectorXd inValues(5);
  inValues(0) = 1.0;
  inValues(1) = 2.0;
  inValues(2) = 3.0;
  inValues(3) = 4.0;
  inValues(4) = 5.0;

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));

  // One big tetra
  Vertex &outVertex0 = outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  Vertex &outVertex1 = outMesh->createVertex(Eigen::Vector3d{0.5, 1.0, 0.0});
  Vertex &outVertex2 = outMesh->createVertex(Eigen::Vector3d{1.0, 0.0, 0.0});
  Vertex &outVertex3 = outMesh->createVertex(Eigen::Vector3d{0.5, 0.0, 1.0});

  outMesh->createTetrahedron(outVertex0, outVertex1, outVertex2, outVertex3);

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::NearestNeighborMapping mapping(mapping::Mapping::SCALED_CONSISTENT_VOLUME, dimensions);

  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  Eigen::VectorXd outValues(5);
  mapping.computeMapping();
  time::Sample inSample(1, inValues);
  mapping.map(inSample, outValues);

  BOOST_TEST(mapping.hasComputedMapping() == true);

  auto inputIntegral  = mesh::integrateVolume(inMesh, inValues);
  auto outputIntegral = mesh::integrateVolume(outMesh, outValues);

  Eigen::VectorXd expectedIntegral(1);
  expectedIntegral << 6.5 * 1. / 12;

  BOOST_TEST(inputIntegral(0) == expectedIntegral(0));

  for (int dim = 0; dim < inputIntegral.size(); ++dim) {
    BOOST_TEST(inputIntegral(dim) == outputIntegral(dim));
  }

  double scaleFactor = outValues(0) / inValues(0);
  BOOST_TEST(scaleFactor != 1.0);

  BOOST_TEST(math::equals(inValues(0) * scaleFactor, outValues(0)));
  BOOST_TEST(inValues(1) * scaleFactor == outValues(1));
  BOOST_TEST(inValues(2) * scaleFactor == outValues(2));
  BOOST_TEST(inValues(3) * scaleFactor == outValues(3));
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

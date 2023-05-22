#include <Eigen/Core>
#include <algorithm>
#include <memory>
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
  PtrData inDataScalar   = inMesh->createData("InDataScalar", 1, 0_dataID);
  PtrData inDataVector   = inMesh->createData("InDataVector", 2, 1_dataID);
  int     inDataScalarID = inDataScalar->getID();
  int     inDataVectorID = inDataVector->getID();
  Vertex &inVertex0      = inMesh->createVertex(Eigen::Vector2d::Constant(0.0));
  Vertex &inVertex1      = inMesh->createVertex(Eigen::Vector2d::Constant(1.0));
  inMesh->allocateDataValues();
  Eigen::VectorXd &inValuesScalar = inDataScalar->values();
  Eigen::VectorXd &inValuesVector = inDataVector->values();
  inValuesScalar << 1.0, 2.0;
  inValuesVector << 1.0, 2.0, 3.0, 4.0;

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  PtrData outDataScalar   = outMesh->createData("OutDataScalar", 1, 2_dataID);
  PtrData outDataVector   = outMesh->createData("OutDataVector", 2, 3_dataID);
  int     outDataScalarID = outDataScalar->getID();
  int     outDataVectorID = outDataVector->getID();
  Vertex &outVertex0      = outMesh->createVertex(Eigen::Vector2d::Constant(0.0));
  Vertex &outVertex1      = outMesh->createVertex(Eigen::Vector2d::Constant(1.0));
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::NearestNeighborMapping mapping(mapping::Mapping::CONSISTENT, dimensions);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Map data with coinciding vertices, has to result in equal values.
  mapping.computeMapping();
  mapping.map(inDataScalarID, outDataScalarID);
  const Eigen::VectorXd &outValuesScalar = outDataScalar->values();
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValuesScalar(0) == inValuesScalar(0));
  BOOST_TEST(outValuesScalar(1) == inValuesScalar(1));
  mapping.map(inDataVectorID, outDataVectorID);
  const Eigen::VectorXd &outValuesVector = outDataVector->values();
  BOOST_CHECK(equals(inValuesVector, outValuesVector));

  // Map data with almost coinciding vertices, has to result in equal values.
  inVertex0.setCoords(outVertex0.getCoords() + Eigen::Vector2d::Constant(0.1));
  inVertex1.setCoords(outVertex1.getCoords() + Eigen::Vector2d::Constant(0.1));
  mapping.computeMapping();
  mapping.map(inDataScalarID, outDataScalarID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValuesScalar(0) == inValuesScalar(0));
  BOOST_TEST(outValuesScalar(1) == inValuesScalar(1));
  mapping.map(inDataVectorID, outDataVectorID);
  BOOST_CHECK(equals(inValuesVector, outValuesVector));

  // Map data with exchanged vertices, has to result in exchanged values.
  inVertex0.setCoords(outVertex1.getCoords());
  inVertex1.setCoords(outVertex0.getCoords());
  mapping.computeMapping();
  mapping.map(inDataScalarID, outDataScalarID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValuesScalar(1) == inValuesScalar(0));
  BOOST_TEST(outValuesScalar(0) == inValuesScalar(1));
  mapping.map(inDataVectorID, outDataVectorID);
  Eigen::Vector4d expected(3.0, 4.0, 1.0, 2.0);
  BOOST_CHECK(equals(expected, outValuesVector));

  // Map data with coinciding output vertices, has to result in same values.
  outVertex1.setCoords(outVertex0.getCoords());
  mapping.computeMapping();
  mapping.map(inDataScalarID, outDataScalarID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValuesScalar(1) == inValuesScalar(1));
  BOOST_TEST(outValuesScalar(0) == inValuesScalar(1));
  mapping.map(inDataVectorID, outDataVectorID);
  expected << 3.0, 4.0, 3.0, 4.0;
  BOOST_CHECK(equals(expected, outValuesVector));
}

BOOST_AUTO_TEST_CASE(ConservativeNonIncremental)
{
  PRECICE_TEST(1_rank);
  int dimensions = 2;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  PtrData inData    = inMesh->createData("InData", 1, 0_dataID);
  int     inDataID  = inData->getID();
  Vertex &inVertex0 = inMesh->createVertex(Eigen::Vector2d::Constant(0.0));
  Vertex &inVertex1 = inMesh->createVertex(Eigen::Vector2d::Constant(1.0));
  inMesh->allocateDataValues();
  Eigen::VectorXd &inValues = inData->values();
  inValues(0)               = 1.0;
  inValues(1)               = 2.0;

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  PtrData outData    = outMesh->createData("OutData", 1, 1_dataID);
  int     outDataID  = outData->getID();
  Vertex &outVertex0 = outMesh->createVertex(Eigen::Vector2d::Constant(0.0));
  Vertex &outVertex1 = outMesh->createVertex(Eigen::Vector2d::Constant(1.0));
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::NearestNeighborMapping mapping(mapping::Mapping::CONSERVATIVE, dimensions);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Map data with coinciding vertices, has to result in equal values.
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  Eigen::VectorXd &outValues = outData->values();
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValues(0) == inValues(0));
  BOOST_TEST(outValues(1) == inValues(1));
  outValues = Eigen::VectorXd::Constant(outValues.size(), 0.0);

  // Map data with almost coinciding vertices, has to result in equal values.
  inVertex0.setCoords(outVertex0.getCoords() + Eigen::Vector2d::Constant(0.1));
  inVertex1.setCoords(outVertex1.getCoords() + Eigen::Vector2d::Constant(0.1));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValues(0) == inValues(0));
  BOOST_TEST(outValues(1) == inValues(1));
  outValues = Eigen::VectorXd::Constant(outValues.size(), 0.0);

  // Map data with exchanged vertices, has to result in exchanged values.
  inVertex0.setCoords(outVertex1.getCoords());
  inVertex1.setCoords(outVertex0.getCoords());
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValues(1) == inValues(0));
  BOOST_TEST(outValues(0) == inValues(1));
  outValues = Eigen::VectorXd::Constant(outValues.size(), 0.0);

  // Map data with coinciding output vertices, has to result in double values.
  outVertex1.setCoords(Eigen::Vector2d::Constant(-1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
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
  PtrData inData    = inMesh->createData("InData", 1, 0_dataID);
  int     inDataID  = inData->getID();
  Vertex &inVertex0 = inMesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  Vertex &inVertex1 = inMesh->createVertex(Eigen::Vector2d{1.0, 0.0});
  Vertex &inVertex2 = inMesh->createVertex(Eigen::Vector2d{3.0, 0.0});
  Vertex &inVertex3 = inMesh->createVertex(Eigen::Vector2d{6.0, 0.0});

  inMesh->createEdge(inVertex0, inVertex1);
  inMesh->createEdge(inVertex1, inVertex2);
  inMesh->createEdge(inVertex2, inVertex3);

  inMesh->allocateDataValues();
  Eigen::VectorXd &inValues = inData->values();
  inValues(0)               = 1.0;
  inValues(1)               = 2.0;
  inValues(2)               = 3.0;
  inValues(3)               = 4.0;

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  PtrData outData    = outMesh->createData("OutData", 1, 1_dataID);
  int     outDataID  = outData->getID();
  Vertex &outVertex0 = outMesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  Vertex &outVertex1 = outMesh->createVertex(Eigen::Vector2d(0.8, 0.0));
  Vertex &outVertex2 = outMesh->createVertex(Eigen::Vector2d(3.0, 0.0));
  Vertex &outVertex3 = outMesh->createVertex(Eigen::Vector2d(6.2, 0.0));

  outMesh->createEdge(outVertex0, outVertex1);
  outMesh->createEdge(outVertex1, outVertex2);
  outMesh->createEdge(outVertex2, outVertex3);

  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::NearestNeighborMapping mapping(mapping::Mapping::SCALED_CONSISTENT_SURFACE, dimensions);

  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  mapping.computeMapping();
  mapping.map(inDataID, outDataID);

  Eigen::VectorXd &outValues = outData->values();
  BOOST_TEST(mapping.hasComputedMapping() == true);

  auto inputIntegral  = mesh::integrateSurface(inMesh, inData->values());
  auto outputIntegral = mesh::integrateSurface(outMesh, outData->values());

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

BOOST_AUTO_TEST_CASE(ScaledConsistentVolume2D)
{
  PRECICE_TEST(1_rank);
  int dimensions = 2;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  PtrData inData   = inMesh->createData("InData", 1, 0_dataID);
  int     inDataID = inData->getID();

  // One square with 3 triangles
  Vertex &inVertex0 = inMesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  Vertex &inVertex1 = inMesh->createVertex(Eigen::Vector2d{1.0, 0.0});
  Vertex &inVertex2 = inMesh->createVertex(Eigen::Vector2d{1.0, 1.0});
  Vertex &inVertex3 = inMesh->createVertex(Eigen::Vector2d{0.0, 1.0});
  Vertex &inVertex4 = inMesh->createVertex(Eigen::Vector2d{0.5, 1.0});

  inMesh->createTriangle(inVertex0, inVertex1, inVertex4);
  inMesh->createTriangle(inVertex0, inVertex3, inVertex4);
  inMesh->createTriangle(inVertex1, inVertex2, inVertex4);

  inMesh->allocateDataValues();
  Eigen::VectorXd &inValues = inData->values();
  inValues(0)               = 1.0;
  inValues(1)               = 2.0;
  inValues(2)               = 3.0;
  inValues(3)               = 4.0;
  inValues(4)               = 5.0;

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  PtrData outData   = outMesh->createData("OutData", 1, 1_dataID);
  int     outDataID = outData->getID();

  // Unit square as 2 triangles
  Vertex &outVertex0 = outMesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  Vertex &outVertex1 = outMesh->createVertex(Eigen::Vector2d(1.0, 0.0));
  Vertex &outVertex2 = outMesh->createVertex(Eigen::Vector2d(1.0, 1.0));
  Vertex &outVertex3 = outMesh->createVertex(Eigen::Vector2d(0.0, 1.0));

  outMesh->createTriangle(outVertex0, outVertex1, outVertex2);
  outMesh->createTriangle(outVertex0, outVertex2, outVertex3);

  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::NearestNeighborMapping mapping(mapping::Mapping::SCALED_CONSISTENT_VOLUME, dimensions);

  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  mapping.computeMapping();
  mapping.map(inDataID, outDataID);

  Eigen::VectorXd &outValues = outData->values();
  BOOST_TEST(mapping.hasComputedMapping() == true);

  auto inputIntegral  = mesh::integrateVolume(inMesh, inData->values());
  auto outputIntegral = mesh::integrateVolume(outMesh, outData->values());

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
  PtrData inData   = inMesh->createData("InData", 1, 0_dataID);
  int     inDataID = inData->getID();

  // One tetra on "out". The "in" has the same but split into two
  Vertex &inVertex0 = inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  Vertex &inVertex1 = inMesh->createVertex(Eigen::Vector3d{0.5, 1.0, 0.0});
  Vertex &inVertex2 = inMesh->createVertex(Eigen::Vector3d{1.0, 0.0, 0.0});
  Vertex &inVertex3 = inMesh->createVertex(Eigen::Vector3d{0.5, 0.0, 1.0});
  Vertex &inVertex4 = inMesh->createVertex(Eigen::Vector3d{0.5, 0.0, 0.0});

  inMesh->createTetrahedron(inVertex0, inVertex1, inVertex3, inVertex4);
  inMesh->createTetrahedron(inVertex1, inVertex2, inVertex3, inVertex4);

  inMesh->allocateDataValues();
  Eigen::VectorXd &inValues = inData->values();
  inValues(0)               = 1.0;
  inValues(1)               = 2.0;
  inValues(2)               = 3.0;
  inValues(3)               = 4.0;
  inValues(4)               = 5.0;

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  PtrData outData   = outMesh->createData("OutData", 1, 1_dataID);
  int     outDataID = outData->getID();

  // One big tetra
  Vertex &outVertex0 = outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  Vertex &outVertex1 = outMesh->createVertex(Eigen::Vector3d{0.5, 1.0, 0.0});
  Vertex &outVertex2 = outMesh->createVertex(Eigen::Vector3d{1.0, 0.0, 0.0});
  Vertex &outVertex3 = outMesh->createVertex(Eigen::Vector3d{0.5, 0.0, 1.0});

  outMesh->createTetrahedron(outVertex0, outVertex1, outVertex2, outVertex3);

  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::NearestNeighborMapping mapping(mapping::Mapping::SCALED_CONSISTENT_VOLUME, dimensions);

  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  mapping.computeMapping();
  mapping.map(inDataID, outDataID);

  Eigen::VectorXd &outValues = outData->values();
  BOOST_TEST(mapping.hasComputedMapping() == true);

  auto inputIntegral  = mesh::integrateVolume(inMesh, inData->values());
  auto outputIntegral = mesh::integrateVolume(outMesh, outData->values());

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

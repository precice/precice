#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>
#include "logging/Logger.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/PartitionOfUnityMapping.hpp"
#include "mapping/impl/BasisFunctions.hpp"
#include "mapping/impl/MappingDataCache.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Utils.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::mesh;
using namespace precice::mapping;
using namespace precice::testing;
using precice::testing::TestContext;

BOOST_AUTO_TEST_SUITE(MappingTests)
BOOST_AUTO_TEST_SUITE(PartitionOfUnityMapping)

void addGlobalIndex(mesh::PtrMesh &mesh, int offset = 0)
{
  for (mesh::Vertex &v : mesh->vertices()) {
    v.setGlobalIndex(v.getID() + offset);
  }
}

double sumComponentWise(const Eigen::VectorXd &vec, int component, int dataDimension)
{
  PRECICE_ASSERT(component < dataDimension);
  double sum = 0;
  for (unsigned int i = 0; i < vec.size() / dataDimension; ++i) {
    sum += vec[(i * dataDimension) + component];
  }
  return sum;
}

BOOST_AUTO_TEST_SUITE(Serial)

void perform2DTestConsistentMapping(Mapping &mapping)
{
  int dimensions = 2;
  using Eigen::Vector2d;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", 1, 0_dataID);
  int           inDataID = inData->getID();
  inMesh->createVertex(Vector2d(0.0, 0.0));
  inMesh->createVertex(Vector2d(1.0, 0.0));
  inMesh->createVertex(Vector2d(1.0, 1.0));
  inMesh->createVertex(Vector2d(0.0, 1.0));

  inMesh->createVertex(Vector2d(2.0, 0.0));
  inMesh->createVertex(Vector2d(3.0, 0.0));
  inMesh->createVertex(Vector2d(3.0, 1.0));
  inMesh->createVertex(Vector2d(2.0, 1.0));

  inMesh->createVertex(Vector2d(4.0, 0.0));
  inMesh->createVertex(Vector2d(5.0, 0.0));
  inMesh->createVertex(Vector2d(5.0, 1.0));
  inMesh->createVertex(Vector2d(4.0, 1.0));

  inMesh->allocateDataValues();
  addGlobalIndex(inMesh);

  auto &values = inData->values();
  values << 1.0, 2.0, 2.0, 1.0, 3.0, 4.0, 4.0, 3.0, 5.0, 6.0, 6.0, 5.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", 1, 1_dataID);
  int           outDataID = outData->getID();
  mesh::Vertex &vertex    = outMesh->createVertex(Vector2d(0, 0));
  mesh::Vertex &vertex1   = outMesh->createVertex(Vector2d(3.5, 0.5));
  mesh::Vertex &vertex2   = outMesh->createVertex(Vector2d(2.5, 0.5));
  // vertex will be changed
  outMesh->createVertex(Vector2d(0, 0));
  outMesh->createVertex(Vector2d(5, 1));
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  vertex.setCoords(Vector2d(0.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  double value  = outData->values()(0);
  double value1 = outData->values()(1);
  double value2 = outData->values()(2);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.0);
  BOOST_TEST(value1 == 4.5);
  BOOST_TEST(value2 == 3.5);

  vertex.setCoords(Vector2d(0.0, 0.5));
  mapping.clear();
  outData->values().setZero();
  BOOST_TEST(mapping.hasComputedMapping() == false);
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Vector2d(0.0, 1.0));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Vector2d(1.0, 0.0));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Vector2d(1.0, 0.5));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Vector2d(1.0, 1.0));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Vector2d(0.5, 0.0));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.5);

  vertex.setCoords(Vector2d(0.5, 0.5));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.5);

  vertex.setCoords(Vector2d(0.5, 1.0));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.5);
}

void perform2DTestConsistentMappingVector(Mapping &mapping)
{
  int dimensions     = 2;
  int dataDimensions = dimensions;
  using Eigen::Vector2d;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", dataDimensions, 0_dataID);
  int           inDataID = inData->getID();
  inMesh->createVertex(Vector2d(0.0, 0.0));
  inMesh->createVertex(Vector2d(1.0, 0.0));
  inMesh->createVertex(Vector2d(1.0, 1.0));
  inMesh->createVertex(Vector2d(0.0, 1.0));

  inMesh->createVertex(Vector2d(2.0, 0.0));
  inMesh->createVertex(Vector2d(3.0, 0.0));
  inMesh->createVertex(Vector2d(3.0, 1.0));
  inMesh->createVertex(Vector2d(2.0, 1.0));

  inMesh->createVertex(Vector2d(4.0, 0.0));
  inMesh->createVertex(Vector2d(5.0, 0.0));
  inMesh->createVertex(Vector2d(5.0, 1.0));
  inMesh->createVertex(Vector2d(4.0, 1.0));

  inMesh->allocateDataValues();
  addGlobalIndex(inMesh);

  auto &values = inData->values();
  // We select the second component = 2 * first component
  values << 1.0, 2.0, 2.0, 4.0, 2.0, 4.0, 1.0, 2.0, 3.0, 6.0, 4.0, 8.0, 4.0, 8.0, 3.0, 6.0, 5.0, 10.0, 6.0, 12.0, 6.0, 12.0, 5.0, 10.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", dataDimensions, 1_dataID);
  int           outDataID = outData->getID();
  mesh::Vertex &vertex    = outMesh->createVertex(Vector2d(0, 0));
  mesh::Vertex &vertex1   = outMesh->createVertex(Vector2d(3.5, 0.5));
  mesh::Vertex &vertex2   = outMesh->createVertex(Vector2d(2.5, 0.5));
  // vertex will be changed
  outMesh->createVertex(Vector2d(0, 0));
  outMesh->createVertex(Vector2d(5, 1));
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  vertex.setCoords(Vector2d(0.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  DataID index  = 0;
  DataID index1 = 1 * dataDimensions;
  DataID index2 = 2 * dataDimensions;
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values()(index) == 1.0);
  BOOST_TEST(outData->values()(index + 1) == 2.0);
  BOOST_TEST(outData->values()(index1) == 4.5);
  BOOST_TEST(outData->values()(index1 + 1) == 9.0);
  BOOST_TEST(outData->values()(index2) == 3.5);
  BOOST_TEST(outData->values()(index2 + 1) == 7.0);

  vertex.setCoords(Vector2d(0.0, 0.5));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values()(0) == 1.0);
  BOOST_TEST(outData->values()(1) == 2.0);

  vertex.setCoords(Vector2d(0.0, 1.0));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values()(0) == 1.0);
  BOOST_TEST(outData->values()(1) == 2.0);

  vertex.setCoords(Vector2d(1.0, 0.0));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values()(0) == 2.0);
  BOOST_TEST(outData->values()(1) == 4.0);

  vertex.setCoords(Vector2d(1.0, 0.5));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values()(0) == 2.0);
  BOOST_TEST(outData->values()(1) == 4.0);

  vertex.setCoords(Vector2d(1.0, 1.0));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values()(0) == 2.0);
  BOOST_TEST(outData->values()(1) == 4.0);

  vertex.setCoords(Vector2d(0.5, 0.0));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values()(0) == 1.5);
  BOOST_TEST(outData->values()(1) == 3.0);

  vertex.setCoords(Vector2d(0.5, 0.5));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values()(0) == 1.5);
  BOOST_TEST(outData->values()(1) == 3.0);

  vertex.setCoords(Vector2d(0.5, 1.0));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values()(0) == 1.5);
  BOOST_TEST(outData->values()(1) == 3.0);
}

// Tests the automatic reduction of dead axis
void performTestConsistentMapDeadAxis(Mapping &mapping, int dim)
{
  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dim, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", 1, 0_dataID);
  int           inDataID = inData->getID();

  // create 20 vertices with two very close lines
  for (unsigned int i = 0; i < 20; ++i)
    if (dim == 2) {
      inMesh->createVertex(Eigen::Vector2d(0.0, 0.0 + i * dim));
      inMesh->createVertex(Eigen::Vector2d(1e-5, 0.0 + i * dim));
    } else {
      inMesh->createVertex(Eigen::Vector3d(7.0, 7.0, 40 + i * dim));
      inMesh->createVertex(Eigen::Vector3d(7.001, 7.001, 40 + i * dim));
    }

  inMesh->allocateDataValues();
  addGlobalIndex(inMesh);

  auto &values = inData->values();
  for (std::size_t v = 0; v < 20; ++v) {
    values[v * 2]     = v * dim * 1e-5;
    values[v * 2 + 1] = v * dim * 1e-5 + 2;
  }

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dim, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", 1, 1_dataID);
  int           outDataID = outData->getID();

  for (unsigned int i = 0; i < 15; ++i) {
    if (dim == 2)
      outMesh->createVertex(Eigen::Vector2d(.1, 1 + 1 + i * dim));
    else {
      outMesh->createVertex(Eigen::Vector3d(7.4, 7.4, 41 + i * dim));
    }
  }
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  double value  = outData->values()(0);
  double value1 = outData->values()(1);
  double value2 = outData->values()(2);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  if (dim == 2) {
    BOOST_TEST(math::equals(value, 19935.374757794027, 0.2));
    BOOST_TEST(math::equals(value1, 19935.37795225389, 0.2));
    BOOST_TEST(math::equals(value2, 19935.35164938603, 0.2));
  } else {
    BOOST_TEST(value == 829.28063055069435);
    BOOST_TEST(value1 == 819.35577218983303);
    BOOST_TEST(value2 == 829.38388713302811);
  }
}

void perform2DTestConservativeMapping(Mapping &mapping)
{
  int dimensions = 2;
  using Eigen::Vector2d;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", 1, 0_dataID);
  int           inDataID = inData->getID();

  // vertices will be changed
  mesh::Vertex &vertex  = inMesh->createVertex(Vector2d(1.0, 0.0));
  mesh::Vertex &vertex1 = inMesh->createVertex(Vector2d(3.5, 0.0));
  mesh::Vertex &vertex2 = inMesh->createVertex(Vector2d(2.5, 0.5));

  inMesh->createVertex(Vector2d(0, 0));
  inMesh->createVertex(Vector2d(5, 1));
  inMesh->allocateDataValues();
  addGlobalIndex(inMesh);
  inMesh->setGlobalNumberOfVertices(inMesh->nVertices());

  auto &values = inData->values();
  values << 10.0, 0.0, 0.0, 10.0, 0.0;
  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", 1, 1_dataID);
  int           outDataID = outData->getID();

  outMesh->createVertex(Vector2d(0.0, 0.0));
  outMesh->createVertex(Vector2d(1.0, 0.0));
  outMesh->createVertex(Vector2d(1.0, 1.0));
  outMesh->createVertex(Vector2d(0.0, 1.0));

  outMesh->createVertex(Vector2d(2.0, 0.0));
  outMesh->createVertex(Vector2d(3.0, 0.0));
  outMesh->createVertex(Vector2d(3.0, 1.0));
  outMesh->createVertex(Vector2d(2.0, 1.0));

  outMesh->createVertex(Vector2d(4.0, 0.0));
  outMesh->createVertex(Vector2d(5.0, 0.0));
  outMesh->createVertex(Vector2d(5.0, 1.0));
  outMesh->createVertex(Vector2d(4.0, 1.0));

  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);
  outMesh->setGlobalNumberOfVertices(outMesh->nVertices());

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Test that we get point-wise exact data values on output vertices 1 and 2
  // remaining vertices are zero
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);

  BOOST_TEST(outData->values().sum() == inData->values().sum());
  BOOST_TEST(outData->values()(0) == 10);
  BOOST_TEST(outData->values()(1) == 10);
  BOOST_TEST(outData->values()(2) == 0);
  BOOST_TEST(outData->values()(3) == 0);

  // Test that we have a symmetric solution if we have a force acting
  // in between two output vertices (here 5 and 8), there is actually
  // no guarantee, but it holds for this setup
  values << 0.0, 10.0, 0.0, 0.0, 0.0;

  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(outData->values().sum() == inData->values().sum());
  BOOST_TEST(outData->values()(5) == outData->values()(8));

  // Test the conservation property if we have everywhere non-zero input data
  values << 5.0, 10.0, 7.0, 3.0, 4.0;

  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);

  BOOST_TEST(outData->values().sum() == inData->values().sum());
  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Test the conservation property if we have everywhere non-zero input data
  values << 3.0, 4.0, 5.0, 7.0, 9.0;

  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(outData->values().sum() == inData->values().sum());
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values()(5) == 3.2931869752057001);
}

void perform2DTestConservativeMappingVector(Mapping &mapping)
{
  int dimensions     = 2;
  int dataDimensions = dimensions;
  using Eigen::Vector2d;

  // Create mesh to map to
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", dataDimensions, 0_dataID);
  DataID        inDataID = inData->getID();
  mesh::Vertex &vertex   = inMesh->createVertex(Vector2d(0, 0));
  mesh::Vertex &vertex1  = inMesh->createVertex(Vector2d(3.5, 0.5));
  mesh::Vertex &vertex2  = inMesh->createVertex(Vector2d(2.5, 0.5));
  // vertex will be changed
  inMesh->createVertex(Vector2d(0, 0));
  inMesh->createVertex(Vector2d(5, 1));

  inMesh->allocateDataValues();
  addGlobalIndex(inMesh);

  auto &values = inData->values();
  // We select the second component = 2 * first component
  values << 1.0, 2.0, 2.0, 4.0, 2.0, 4.0, 1.0, 2.0, 3.0, 6.0;

  // Create mesh to map from
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", dataDimensions, 1_dataID);
  DataID        outDataID = outData->getID();

  outMesh->createVertex(Vector2d(0.0, 0.0));
  outMesh->createVertex(Vector2d(1.0, 0.0));
  outMesh->createVertex(Vector2d(1.0, 1.0));
  outMesh->createVertex(Vector2d(0.0, 1.0));

  outMesh->createVertex(Vector2d(2.0, 0.0));
  outMesh->createVertex(Vector2d(3.0, 0.0));
  outMesh->createVertex(Vector2d(3.0, 1.0));
  outMesh->createVertex(Vector2d(2.0, 1.0));

  outMesh->createVertex(Vector2d(4.0, 0.0));
  outMesh->createVertex(Vector2d(5.0, 0.0));
  outMesh->createVertex(Vector2d(5.0, 1.0));
  outMesh->createVertex(Vector2d(4.0, 1.0));

  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Check that we preserve data component-wise
  vertex.setCoords(Vector2d(1.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);

  // component 0
  double insumc0  = sumComponentWise(inData->values(), 0, dataDimensions);
  double outsumc0 = sumComponentWise(outData->values(), 0, dataDimensions);
  BOOST_TEST(insumc0 == outsumc0);

  // component 1
  double insumc1  = sumComponentWise(inData->values(), 1, dataDimensions);
  double outsumc1 = sumComponentWise(outData->values(), 1, dataDimensions);

  BOOST_TEST(insumc0 * 2 == insumc1);
  BOOST_TEST(outsumc0 * 2 == outsumc1);
  BOOST_TEST(insumc1 == outsumc1);

  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Check for specific sum values for each component
  values << 1.0, 0.0, 12.0, 0.0, 2.0, 0.0, 1.0, 0.0, 3.0, 0.0;

  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);

  // component 0
  outsumc0 = sumComponentWise(outData->values(), 0, dataDimensions);
  BOOST_TEST(19.0 == outsumc0);

  // we expect no contribution for the zeroth component
  // component 1
  outsumc1 = sumComponentWise(outData->values(), 1, dataDimensions);
  BOOST_TEST(0.0 == outsumc1);

  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Check for specific sum values for each component
  values << 0.0, 2.0, 0.0, 4.0, 0.0, 8.0, 0.0, 7.0, 0.0, 0.0;

  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);

  // component 0
  outsumc0 = sumComponentWise(outData->values(), 0, dataDimensions);
  BOOST_TEST(0.0 == outsumc0);

  // we expect no contribution for the zeroth component
  // component 1
  outsumc1 = sumComponentWise(outData->values(), 1, dataDimensions);
  BOOST_TEST(21.0 == outsumc1);

  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Check for the exact reproduction of matching vertices
  values << 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 27.0;

  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);

  // component 0
  outsumc0 = sumComponentWise(outData->values(), 0, dataDimensions);
  BOOST_TEST(10.0 == outsumc0);
  BOOST_TEST(outData->values()(2) == 10);

  // component 1
  outsumc1 = sumComponentWise(outData->values(), 1, dataDimensions);
  BOOST_TEST(27.0 == outsumc1);
  BOOST_TEST(outData->values()(21) == 27);

  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Check for the interpolation at some vertex
  values << 3.0, 6.0, 2.0, 4.0, 12.0, 24.0, 1.0, 2.0, 3.0, 6.0;

  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);

  BOOST_TEST(outData->values().sum() == inData->values().sum());

  double expectedValue = 3.5933508322619825;
  // component 0
  BOOST_TEST(outData->values()(12) == expectedValue);

  // component 1
  BOOST_TEST(outData->values()(13) == 2 * expectedValue);
}

void perform3DTestConsistentMapping(Mapping &mapping)
{
  int dimensions = 3;
  using Eigen::Vector3d;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", 1, 0_dataID);
  int           inDataID = inData->getID();
  inMesh->createVertex(Vector3d(0.0, 0.0, 0.0));
  inMesh->createVertex(Vector3d(1.0, 0.0, 0.0));
  inMesh->createVertex(Vector3d(1.0, 1.0, 0.0));
  inMesh->createVertex(Vector3d(0.0, 1.0, 0.0));

  inMesh->createVertex(Vector3d(0.0, 0.0, 1.0));
  inMesh->createVertex(Vector3d(1.0, 0.0, 1.0));
  inMesh->createVertex(Vector3d(1.0, 1.0, 1.0));
  inMesh->createVertex(Vector3d(0.0, 1.0, 1.0));

  inMesh->createVertex(Vector3d(2.0, 0.0, 0.0));
  inMesh->createVertex(Vector3d(3.0, 0.0, 0.0));
  inMesh->createVertex(Vector3d(3.0, 1.0, 0.0));
  inMesh->createVertex(Vector3d(2.0, 1.0, 0.0));

  inMesh->createVertex(Vector3d(2.0, 0.0, 1.0));
  inMesh->createVertex(Vector3d(3.0, 0.0, 1.0));
  inMesh->createVertex(Vector3d(3.0, 1.0, 1.0));
  inMesh->createVertex(Vector3d(2.0, 1.0, 1.0));

  inMesh->createVertex(Vector3d(4.0, 0.0, 0.0));
  inMesh->createVertex(Vector3d(5.0, 0.0, 0.0));
  inMesh->createVertex(Vector3d(5.0, 1.0, 0.0));
  inMesh->createVertex(Vector3d(4.0, 1.0, 0.0));

  inMesh->createVertex(Vector3d(4.0, 0.0, 1.0));
  inMesh->createVertex(Vector3d(5.0, 0.0, 1.0));
  inMesh->createVertex(Vector3d(5.0, 1.0, 1.0));
  inMesh->createVertex(Vector3d(4.0, 1.0, 1.0));

  inMesh->allocateDataValues();
  addGlobalIndex(inMesh);

  auto &values = inData->values();
  // Set the values in the parallel (z) plane 3*values
  values << 1.0, 2.0, 2.0, 1.0, 3.0, 6.0, 6.0, 3.0, 3.0, 4.0, 4.0, 3.0, 9.0, 12.0, 12.0, 9.0, 5.0, 6.0, 6.0, 5.0, 15.0, 18.0, 18.0, 15.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", 1, 1_dataID);
  int           outDataID = outData->getID();
  mesh::Vertex &vertex    = outMesh->createVertex(Vector3d(0, 0, 0));
  mesh::Vertex &vertex1   = outMesh->createVertex(Vector3d(3.5, 0.5, 0.5));
  mesh::Vertex &vertex2   = outMesh->createVertex(Vector3d(2.5, 0.5, 1.0));
  outMesh->createVertex(Vector3d(0, 0, 0.5));
  outMesh->createVertex(Vector3d(5, 1, 0.5));
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  vertex.setCoords(Vector3d(0.0, 0.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  double value  = outData->values()(0);
  double value1 = outData->values()(1);
  double value2 = outData->values()(2);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.0);
  BOOST_TEST(value1 == 9.0);
  BOOST_TEST(value2 == 10.5);

  vertex.setCoords(Vector3d(0.0, 0.5, 0.5));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Vector3d(0.0, 1.0, 1.0));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 3.0);

  vertex.setCoords(Vector3d(1.0, 0.0, 0.0));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Vector3d(1.0, 0.5, 0.5));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 4.0);

  vertex.setCoords(Vector3d(1.0, 1.0, 1.0));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 6.0);

  vertex.setCoords(Vector3d(0.5, 0.0, 0.5));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 3.0);

  vertex.setCoords(Vector3d(0.5, 0.5, 1.0));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 4.5);

  vertex.setCoords(Vector3d(0.5, 1.0, 0.0));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.5);
}

// uses scalar data in 3D
void perform3DTestJustInTimeMappingWithPolynomial(Mapping &mapping)
{
  const int dimensions     = 3;
  const int dataComponents = 1;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", dataComponents, 0_dataID);
  int           inDataID = inData->getID();
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(1.0, 1.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));

  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 1.0));
  inMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 1.0));
  inMesh->createVertex(Eigen::Vector3d(1.0, 1.0, 1.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 1.0));

  inMesh->createVertex(Eigen::Vector3d(2.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(3.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(3.0, 1.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(2.0, 1.0, 0.0));

  inMesh->createVertex(Eigen::Vector3d(2.0, 0.0, 1.0));
  inMesh->createVertex(Eigen::Vector3d(3.0, 0.0, 1.0));
  inMesh->createVertex(Eigen::Vector3d(3.0, 1.0, 1.0));
  inMesh->createVertex(Eigen::Vector3d(2.0, 1.0, 1.0));

  inMesh->createVertex(Eigen::Vector3d(4.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(5.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(5.0, 1.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(4.0, 1.0, 0.0));

  inMesh->createVertex(Eigen::Vector3d(4.0, 0.0, 1.0));
  inMesh->createVertex(Eigen::Vector3d(5.0, 0.0, 1.0));
  inMesh->createVertex(Eigen::Vector3d(5.0, 1.0, 1.0));
  inMesh->createVertex(Eigen::Vector3d(4.0, 1.0, 1.0));

  inMesh->allocateDataValues();
  addGlobalIndex(inMesh);

  auto &values = inData->values();

  // Set the values in the parallel (z) plane 3*values
  values << 1.0, 2.0, 2.0, 1.0, 3.0, 6.0, 6.0, 3.0, 3.0, 4.0, 4.0, 3.0, 9.0, 12.0, 12.0, 9.0, 5.0, 6.0, 6.0, 5.0, 15.0, 18.0, 18.0, 15.0;

  // the dummy target mesh
  mesh::PtrMesh toMesh = mesh::MeshConfiguration::getJustInTimeMappingMesh(dimensions);

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, toMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);
  // compute the mapping (affects only the inMesh)
  mapping.computeMapping();
  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Now, we can setup the MappingDataCache with the structures
  // from computeMapping
  mapping::impl::MappingDataCache cache(dataComponents);
  mapping.initializeMappingDataCache(cache);

  // computeMapping and the initializeMappingDataCache are only required for changes in the input coords
  // whereas the cache update is for new data values
  mapping.updateMappingDataCache(cache, values);
  // we can also give the cache a time stamp (not strictly necessary)
  double stamp = 1;
  cache.setTimeStamp(stamp);

  // Test infrastructure
  Eigen::MatrixXd coords(dimensions, 1);
  Eigen::MatrixXd result(dataComponents, 1);

  // Now we can evaluate at any point we want
  coords.setZero();
  // Check that we have the right time
  BOOST_TEST(cache.hasDataAtTimeStamp(stamp));
  mapping.mapConsistentAt(coords, cache, result);
  BOOST_TEST(result(0, 0) == 1.0);

  coords.col(0) = Eigen::RowVector3d(3.5, 0.5, 0.5);
  // Check that we have the right time
  BOOST_TEST(cache.hasDataAtTimeStamp(stamp));
  mapping.mapConsistentAt(coords, cache, result);
  BOOST_TEST(result(0, 0) == 9.0);

  coords.col(0) = Eigen::RowVector3d(2.5, 0.5, 1.0);
  // Check that we have the right time
  BOOST_TEST(cache.hasDataAtTimeStamp(stamp));
  mapping.mapConsistentAt(coords, cache, result);
  BOOST_TEST(result(0, 0) == 10.5);

  // We can also pass multiple vertices at once
  coords.resize(dimensions, 3);
  result.resize(dataComponents, 3);

  coords.col(0) = Eigen::RowVector3d(0.0, 0.5, 0.5);
  coords.col(1) = Eigen::RowVector3d(0.0, 1.0, 1.0);
  coords.col(2) = Eigen::RowVector3d(1.0, 0.0, 0.0);
  // Check that we have the right time
  BOOST_TEST(cache.hasDataAtTimeStamp(stamp));
  mapping.mapConsistentAt(coords, cache, result);
  BOOST_TEST(result(0, 0) == 2.0);
  BOOST_TEST(result(0, 1) == 3.0);
  BOOST_TEST(result(0, 2) == 2.0);

  // Now we alter the data (double values)
  values *= 2;
  // Mapping remains the same
  BOOST_TEST(mapping.hasComputedMapping() == true);
  // Step 1: invalidate cache (happens in the DataContext)
  cache.resetTimeStamp();
  BOOST_TEST(!cache.hasDataAtTimeStamp(stamp));
  // Step 2: update the cache
  mapping.updateMappingDataCache(cache, values);
  // Step 3: mark the cache
  stamp = 2.0;
  cache.setTimeStamp(stamp);
  BOOST_TEST(cache.hasDataAtTimeStamp(stamp));

  // Ready for new mappings
  mapping.mapConsistentAt(coords, cache, result);
  BOOST_TEST(result(0, 0) == 4.0);
  BOOST_TEST(result(0, 1) == 6.0);
  BOOST_TEST(result(0, 2) == 4.0);

  coords.col(0) = Eigen::RowVector3d(0.5, 0.0, 0.5);
  coords.col(1) = Eigen::RowVector3d(0.5, 0.5, 1.0);
  coords.col(2) = Eigen::RowVector3d(0.5, 1.0, 0.0);
  // Check that we have the right time
  BOOST_TEST(cache.hasDataAtTimeStamp(stamp));
  mapping.mapConsistentAt(coords, cache, result);
  BOOST_TEST(result(0, 0) == 6.0);
  BOOST_TEST(result(0, 1) == 9.0);
  BOOST_TEST(result(0, 2) == 3.0);

  // If we clear the mapping, the we have to recompute it
  mapping.clear();
  // We just clear the mesh index here instead of calling mesh.clear
  // to not repeat all the vertices above. Recomputing the index is require
  // when altering the mesh
  inMesh->index().clear();
  // We extend the mesh a bit
  inMesh->createVertex(Eigen::Vector3d(6.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(7.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(7.0, 1.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(6.0, 1.0, 0.0));
  inMesh->allocateDataValues();
  values << 1.0, 2.0, 2.0, 1.0, 3.0, 6.0, 6.0, 3.0, 3.0, 4.0, 4.0, 3.0, 9.0, 12.0, 12.0, 9.0, 5.0, 6.0, 6.0, 5.0, 15.0, 18.0, 18.0, 15.0, 7.0, 8.0, 8.0, 7.0;
  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, toMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);
  // compute the mapping (affects only the inMesh)
  mapping.computeMapping();
  BOOST_TEST(mapping.hasComputedMapping() == true);

  cache.resetTimeStamp();
  mapping.initializeMappingDataCache(cache);
  mapping.updateMappingDataCache(cache, values);
  cache.setTimeStamp(stamp);
  coords.col(0) = Eigen::RowVector3d(6.5, 1.0, 0.0);
  mapping.mapConsistentAt(coords, cache, result);
  BOOST_TEST(result(0, 0) == 7.5);
  BOOST_TEST(result(0, 1) == 4.5);
  BOOST_TEST(result(0, 2) == 1.5);
}

// uses vectorial data in 2D
void perform2DTestJustInTimeMappingNoPolynomial(Mapping &mapping)
{
  int dimensions     = 2;
  int dataComponents = 2;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", dataComponents, 0_dataID);
  int           inDataID = inData->getID();
  inMesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  inMesh->createVertex(Eigen::Vector2d(1.0, 0.0));
  inMesh->createVertex(Eigen::Vector2d(1.0, 1.0));
  inMesh->createVertex(Eigen::Vector2d(0.0, 1.0));

  inMesh->createVertex(Eigen::Vector2d(2.0, 0.0));
  inMesh->createVertex(Eigen::Vector2d(3.0, 0.0));
  inMesh->createVertex(Eigen::Vector2d(3.0, 1.0));
  inMesh->createVertex(Eigen::Vector2d(2.0, 1.0));

  inMesh->createVertex(Eigen::Vector2d(4.0, 0.0));
  inMesh->createVertex(Eigen::Vector2d(5.0, 0.0));
  inMesh->createVertex(Eigen::Vector2d(5.0, 1.0));
  inMesh->createVertex(Eigen::Vector2d(4.0, 1.0));

  inMesh->allocateDataValues();
  addGlobalIndex(inMesh);

  auto &values = inData->values();

  // Set the values in the parallel (z) plane 3*values
  values << 1.0, 2.0, 2.0, 1.0, 3.0, 6.0, 6.0, 3.0, 3.0, 4.0, 4.0, 3.0, 9.0, 12.0, 12.0, 9.0, 5.0, 6.0, 6.0, 5.0, 15.0, 18.0, 18.0, 15.0;

  // the dummy target mesh
  mesh::PtrMesh toMesh = mesh::MeshConfiguration::getJustInTimeMappingMesh(dimensions);

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, toMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);
  // compute the mapping (affects only the inMesh)
  mapping.computeMapping();
  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Now, we can setup the MappingDataCache with the structures
  // from computeMapping
  mapping::impl::MappingDataCache cache(dataComponents);
  mapping.initializeMappingDataCache(cache);

  // computeMapping and the initializeMappingDataCache are only required for changes in the input coords
  // whereas the cache update is for new data values
  mapping.updateMappingDataCache(cache, values);
  // we can also give the cache a time stamp (not strictly necessary)
  double stamp = 1;
  cache.setTimeStamp(stamp);

  // Test infrastructure
  Eigen::MatrixXd coords(dimensions, 1);
  Eigen::MatrixXd result(dataComponents, 1);

  // Now we can evaluate at any point we want
  coords.setZero();
  // Check that we have the right time
  BOOST_TEST(cache.hasDataAtTimeStamp(stamp));
  mapping.mapConsistentAt(coords, cache, result);
  BOOST_TEST(result(0, 0) == 1.0);
  BOOST_TEST(result(1, 0) == 2.0);

  // Second last point given
  coords.col(0) = Eigen::RowVector2d(5.0, 1.0);
  // Check that we have the right time
  BOOST_TEST(cache.hasDataAtTimeStamp(stamp));
  mapping.mapConsistentAt(coords, cache, result);
  BOOST_TEST(result(0, 0) == 15.0);
  BOOST_TEST(result(1, 0) == 18.0);

  // Between the first and the second point
  coords.col(0) = Eigen::RowVector2d(0.5, 0.0);
  // Check that we have the right time
  BOOST_TEST(cache.hasDataAtTimeStamp(stamp));
  mapping.mapConsistentAt(coords, cache, result);
  BOOST_TEST(result(0, 0) == 1.71978075, boost::test_tools::tolerance(1e-7));
  BOOST_TEST(result(1, 0) == 1.71978075, boost::test_tools::tolerance(1e-7));

  coords.col(0) = Eigen::RowVector2d(3.5, 0.5);
  // Check that we have the right time
  BOOST_TEST(cache.hasDataAtTimeStamp(stamp));
  mapping.mapConsistentAt(coords, cache, result);
  BOOST_TEST(result(0, 0) == 10.39684625, boost::test_tools::tolerance(1e-7));
  BOOST_TEST(result(1, 0) == 10.48101386, boost::test_tools::tolerance(1e-7));

  // We can also pass multiple vertices at once
  coords.resize(dimensions, 2);
  result.resize(dataComponents, 2);

  coords.col(0) = Eigen::RowVector2d(4.0, 0.0);
  coords.col(1) = Eigen::RowVector2d(5.0, 0.0);
  // // Check that we have the right time
  BOOST_TEST(cache.hasDataAtTimeStamp(stamp));
  mapping.mapConsistentAt(coords, cache, result);
  BOOST_TEST(result(0, 0) == 5.0);
  BOOST_TEST(result(1, 0) == 6.0);
  BOOST_TEST(result(0, 1) == 6.0);
  BOOST_TEST(result(1, 1) == 5.0);

  // Now we alter the data (double values)
  values *= 0.5;
  // Mapping remains the same
  BOOST_TEST(mapping.hasComputedMapping() == true);
  // Step 1: invalidate cache (happens in the DataContext)
  cache.resetTimeStamp();
  BOOST_TEST(!cache.hasDataAtTimeStamp(stamp));
  // Step 2: update the cache
  mapping.updateMappingDataCache(cache, values);
  // Step 3: mark the cache
  stamp = 2.0;
  cache.setTimeStamp(stamp);
  BOOST_TEST(cache.hasDataAtTimeStamp(stamp));

  // Ready for new mappings
  mapping.mapConsistentAt(coords, cache, result);
  BOOST_TEST(result(0, 0) == 2.5);
  BOOST_TEST(result(1, 0) == 3.0);
  BOOST_TEST(result(0, 1) == 3.0);
  BOOST_TEST(result(1, 1) == 2.5);

  // If we clear the mapping, the we have to recompute it
  mapping.clear();
  // We just clear the mesh index here instead of calling mesh.clear
  // to not repeat all the vertices above. Recomputing the index is require
  // when altering the mesh
  inMesh->index().clear();
  // We extend the mesh a bit
  inMesh->createVertex(Eigen::Vector2d(6.0, 0.0));
  inMesh->createVertex(Eigen::Vector2d(7.0, 0.0));
  inMesh->createVertex(Eigen::Vector2d(7.0, 1.0));
  inMesh->createVertex(Eigen::Vector2d(6.0, 1.0));
  inMesh->allocateDataValues();
  values << 1.0, 2.0, 2.0, 1.0, 3.0, 6.0, 6.0, 3.0, 3.0, 4.0, 4.0, 3.0, 9.0, 12.0, 12.0, 9.0, 5.0, 6.0, 6.0, 5.0, 15.0, 18.0, 18.0, 15.0, 7.0, 8.0, 8.0, 7.0, 20.0, 22.0, 22.0, 20.0;
  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, toMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);
  // compute the mapping (affects only the inMesh)
  mapping.computeMapping();
  BOOST_TEST(mapping.hasComputedMapping() == true);

  cache.resetTimeStamp();
  mapping.initializeMappingDataCache(cache);
  mapping.updateMappingDataCache(cache, values);
  cache.setTimeStamp(stamp);

  mapping.mapConsistentAt(coords, cache, result);
  BOOST_TEST(result(0, 0) == 5.0);
  BOOST_TEST(result(1, 0) == 6.0);
  BOOST_TEST(result(0, 1) == 6.0);
  BOOST_TEST(result(1, 1) == 5.0);

  coords.col(0) = Eigen::RowVector2d(7.0, 1.0);
  coords.col(1) = Eigen::RowVector2d(6.0, 0.0);
  mapping.mapConsistentAt(coords, cache, result);
  BOOST_TEST(result(0, 0) == 20.0);
  BOOST_TEST(result(1, 0) == 22.0);
  BOOST_TEST(result(0, 1) == 7.0);
  BOOST_TEST(result(1, 1) == 8.0);
}

void perform3DTestConsistentMappingVector(Mapping &mapping)
{
  int dimensions     = 3;
  int dataDimensions = dimensions;

  using Eigen::Vector3d;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", dataDimensions, 0_dataID);
  int           inDataID = inData->getID();
  inMesh->createVertex(Vector3d(0.0, 0.0, 0.0));
  inMesh->createVertex(Vector3d(1.0, 0.0, 0.0));
  inMesh->createVertex(Vector3d(1.0, 1.0, 0.0));
  inMesh->createVertex(Vector3d(0.0, 1.0, 0.0));

  inMesh->createVertex(Vector3d(0.0, 0.0, 1.0));
  inMesh->createVertex(Vector3d(1.0, 0.0, 1.0));
  inMesh->createVertex(Vector3d(1.0, 1.0, 1.0));
  inMesh->createVertex(Vector3d(0.0, 1.0, 1.0));

  inMesh->createVertex(Vector3d(2.0, 0.0, 0.0));
  inMesh->createVertex(Vector3d(3.0, 0.0, 0.0));
  inMesh->createVertex(Vector3d(3.0, 1.0, 0.0));
  inMesh->createVertex(Vector3d(2.0, 1.0, 0.0));

  inMesh->createVertex(Vector3d(2.0, 0.0, 1.0));
  inMesh->createVertex(Vector3d(3.0, 0.0, 1.0));
  inMesh->createVertex(Vector3d(3.0, 1.0, 1.0));
  inMesh->createVertex(Vector3d(2.0, 1.0, 1.0));

  inMesh->createVertex(Vector3d(4.0, 0.0, 0.0));
  inMesh->createVertex(Vector3d(5.0, 0.0, 0.0));
  inMesh->createVertex(Vector3d(5.0, 1.0, 0.0));
  inMesh->createVertex(Vector3d(4.0, 1.0, 0.0));

  inMesh->createVertex(Vector3d(4.0, 0.0, 1.0));
  inMesh->createVertex(Vector3d(5.0, 0.0, 1.0));
  inMesh->createVertex(Vector3d(5.0, 1.0, 1.0));
  inMesh->createVertex(Vector3d(4.0, 1.0, 1.0));

  inMesh->allocateDataValues();
  addGlobalIndex(inMesh);

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", dataDimensions, 1_dataID);
  int           outDataID = outData->getID();
  mesh::Vertex &vertex    = outMesh->createVertex(Vector3d(0, 0, 0));
  mesh::Vertex &vertex1   = outMesh->createVertex(Vector3d(3.5, 0.5, 0.5));
  mesh::Vertex &vertex2   = outMesh->createVertex(Vector3d(2.5, 0.5, 1.0));
  outMesh->createVertex(Vector3d(0, 0, 0.5));
  outMesh->createVertex(Vector3d(5, 1, 0.5));
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  for (unsigned int i = 0; i < inData->values().size() / dataDimensions; ++i) {
    inData->values()(i * dataDimensions)     = 7;
    inData->values()(i * dataDimensions + 1) = 38;
    inData->values()(i * dataDimensions + 2) = 19;
  }

  // Check for the independent interpolation of each component when applying a constant field
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values()(0) == 7.0);
  BOOST_TEST(outData->values()(1) == 38.0);
  BOOST_TEST(outData->values()(2) == 19.0);
  BOOST_TEST(outData->values()(9) == 7.0);
  BOOST_TEST(outData->values()(10) == 38.0);
  BOOST_TEST(outData->values()(11) == 19.0);

  // Check for the consistency between individual components
  for (unsigned int i = 0; i < inData->values().size() / dataDimensions; ++i) {
    inData->values()(i * dataDimensions)     = std::pow(i * dataDimensions, 2);
    inData->values()(i * dataDimensions + 1) = 2 * std::pow(i * dataDimensions, 2);
    inData->values()(i * dataDimensions + 2) = 3 * std::pow(i * dataDimensions, 2);
  }

  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values()(12) == 3673.7308684337013);
  BOOST_TEST(outData->values()(13) == 2 * outData->values()(12));
  BOOST_TEST(outData->values()(14) == 3 * outData->values()(12));

  // The remaining parts should already be covered by the other 3D/2D tests
}

void perform3DTestConservativeMapping(Mapping &mapping)
{
  int dimensions = 3;
  using Eigen::Vector3d;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", 1, 0_dataID);
  const DataID  inDataID = inData->getID();
  mesh::Vertex &vertex   = inMesh->createVertex(Vector3d(0, 0, 0));
  mesh::Vertex &vertex1  = inMesh->createVertex(Vector3d(3.5, 0.5, 0.5));
  mesh::Vertex &vertex2  = inMesh->createVertex(Vector3d(2.5, 0.5, 1.0));
  inMesh->createVertex(Vector3d(1, 1, 1));
  inMesh->createVertex(Vector3d(5, 1, 0.5));
  inMesh->allocateDataValues();
  addGlobalIndex(inMesh);

  // Create mesh to map from
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", 1, 1_dataID);
  const DataID  outDataID = outData->getID();
  outMesh->createVertex(Vector3d(0.0, 0.0, 0.0));
  outMesh->createVertex(Vector3d(1.0, 0.0, 0.0));
  outMesh->createVertex(Vector3d(1.0, 1.0, 0.0));
  outMesh->createVertex(Vector3d(0.0, 1.0, 0.0));

  outMesh->createVertex(Vector3d(0.0, 0.0, 1.0));
  outMesh->createVertex(Vector3d(1.0, 0.0, 1.0));
  outMesh->createVertex(Vector3d(1.0, 1.0, 1.0));
  outMesh->createVertex(Vector3d(0.0, 1.0, 1.0));

  outMesh->createVertex(Vector3d(2.0, 0.0, 0.0));
  outMesh->createVertex(Vector3d(3.0, 0.0, 0.0));
  outMesh->createVertex(Vector3d(3.0, 1.0, 0.0));
  outMesh->createVertex(Vector3d(2.0, 1.0, 0.0));

  outMesh->createVertex(Vector3d(2.0, 0.0, 1.0));
  outMesh->createVertex(Vector3d(3.0, 0.0, 1.0));
  outMesh->createVertex(Vector3d(3.0, 1.0, 1.0));
  outMesh->createVertex(Vector3d(2.0, 1.0, 1.0));

  outMesh->createVertex(Vector3d(4.0, 0.0, 0.0));
  outMesh->createVertex(Vector3d(5.0, 0.0, 0.0));
  outMesh->createVertex(Vector3d(5.0, 1.0, 0.0));
  outMesh->createVertex(Vector3d(4.0, 1.0, 0.0));

  outMesh->createVertex(Vector3d(4.0, 0.0, 1.0));
  outMesh->createVertex(Vector3d(5.0, 0.0, 1.0));
  outMesh->createVertex(Vector3d(5.0, 1.0, 1.0));
  outMesh->createVertex(Vector3d(4.0, 1.0, 1.0));

  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  auto &values = inData->values();
  // Some arbitrary input values in order to check conservation
  values << 1.0, 2.0, 2.0, 1.0, 3.0;

  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values().sum() == inData->values().sum());

  // Check for the exact reproduction of individual values
  values << 12.0, 5.0, 7.0, 8.0, 9.0;

  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  const double expectedSum = 41;
  BOOST_TEST(outData->values().sum() == expectedSum);
  BOOST_TEST(outData->values()(6) == 8);
  BOOST_TEST(outData->values()(0) == 12);

  // Check for consistency (applying a heavy load to the front layer z = 1)
  values << 0.0, 50.0, 107.0, 108.0, 48.0;

  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values().sum() == inData->values().sum());
  BOOST_TEST(outData->values()(0) == 0.0);
  BOOST_TEST(outData->values()(6) > outData->values()(1));
  BOOST_TEST(outData->values()(13) > outData->values()(9));

  // Check for symmetry when applying a central load (not guaranteed ?)
  values << 0.0, 100.0, 0.0, 0.0, 0.0;

  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values().sum() == inData->values().sum());
  BOOST_TEST(outData->values()(9) == outData->values()(13));
  BOOST_TEST(outData->values()(9) == outData->values()(10));
  BOOST_TEST(outData->values()(16) == outData->values()(19));
  // TODO: There is no symmetry between the (3, x, x) and (4, x , x) layer
  // as the domain length is 5 and we have a partition center at x = 2.9 and x = 3.75
  // BOOST_TEST(outData->values()(16) == outData->values()(9));
  BOOST_TEST(outData->values()(20) == outData->values()(23));
  BOOST_TEST(outData->values()(16) == outData->values()(20));
  BOOST_TEST(outData->values()(9) == 15.470584170385226);
}

void perform3DTestJustInTimeMappingConservative(Mapping &mapping)
{
  const int dimensions     = 3;
  const int dataComponents = 1;

  // the dummy from mesh
  mesh::PtrMesh inMesh = mesh::MeshConfiguration::getJustInTimeMappingMesh(dimensions);

  // Create mesh to map from
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData outData = outMesh->createData("OutData", dataComponents, 1_dataID);

  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(1.0, 1.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));

  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 1.0));
  outMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 1.0));
  outMesh->createVertex(Eigen::Vector3d(1.0, 1.0, 1.0));
  outMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 1.0));

  outMesh->createVertex(Eigen::Vector3d(2.0, 0.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(3.0, 0.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(3.0, 1.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(2.0, 1.0, 0.0));

  outMesh->createVertex(Eigen::Vector3d(2.0, 0.0, 1.0));
  outMesh->createVertex(Eigen::Vector3d(3.0, 0.0, 1.0));
  outMesh->createVertex(Eigen::Vector3d(3.0, 1.0, 1.0));
  outMesh->createVertex(Eigen::Vector3d(2.0, 1.0, 1.0));

  outMesh->createVertex(Eigen::Vector3d(4.0, 0.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(5.0, 0.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(5.0, 1.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(4.0, 1.0, 0.0));

  outMesh->createVertex(Eigen::Vector3d(4.0, 0.0, 1.0));
  outMesh->createVertex(Eigen::Vector3d(5.0, 0.0, 1.0));
  outMesh->createVertex(Eigen::Vector3d(5.0, 1.0, 1.0));
  outMesh->createVertex(Eigen::Vector3d(4.0, 1.0, 1.0));

  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  mapping.computeMapping();
  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Now, we can setup the MappingDataCache with the structures
  // from computeMapping before starting the mappings
  mapping::impl::MappingDataCache cache(dataComponents);
  mapping.initializeMappingDataCache(cache);

  // Test infrastructure/ the user input
  Eigen::MatrixXd             coords(dimensions, 1);
  Eigen::MatrixXd             inData(dataComponents, 1);
  Eigen::Map<Eigen::MatrixXd> outValues(outData->values().data(), dataComponents, outMesh->vertices().size());

  // before writing first, we have to reset the cache to contain zeros
  cache.resetData();
  outValues.setZero();

  coords.col(0) = Eigen::RowVector3d(5.0, 0.0, 1.0);
  inData(0, 0)  = 4.3;
  // note that the outvalues are only filled in the completeJustInTimeMapping
  mapping.mapConservativeAt(coords, inData, cache, outValues);

  // Once all values are written, we complete the mapping (expensive part)
  mapping.completeJustInTimeMapping(cache, outValues);
  BOOST_TEST(outValues.sum() == inData.sum());
  Eigen::VectorXd expectedValues(24);
  expectedValues.setZero();
  expectedValues(21) = 4.3;
  BOOST_TEST(outData->values() == expectedValues, boost::test_tools::per_element());

  // clear data again for the next go
  cache.resetData();
  outValues.setZero();

  // if we call the function twice for the same coordinate, we get twice the contribution
  mapping.mapConservativeAt(coords, inData, cache, outValues);
  mapping.mapConservativeAt(coords, inData, cache, outValues);
  mapping.completeJustInTimeMapping(cache, outValues);
  BOOST_TEST(outValues.sum() == 2 * inData.sum());
  BOOST_TEST(outData->values() == 2 * expectedValues, boost::test_tools::per_element());

  // clear data again for the next go
  cache.resetData();
  outValues.setZero();

  double expectedSum = 0;
  mapping.mapConservativeAt(coords, inData, cache, outValues);
  expectedSum += inData.sum();

  // We can also change the shape of the user input
  coords.resize(dimensions, 3);
  inData.resize(dataComponents, 3);
  coords.col(0) = Eigen::RowVector3d(3.5, 0.5, 0.5);
  coords.col(1) = Eigen::RowVector3d(4.5, 1.0, 1.0);
  coords.col(2) = Eigen::RowVector3d(0.1, 0.0, 1.0);
  inData(0, 0)  = 4.3;
  inData(0, 1)  = 17.3;
  inData(0, 2)  = 5.8;
  mapping.mapConservativeAt(coords, inData, cache, outValues);
  expectedSum += inData.sum();

  mapping.completeJustInTimeMapping(cache, outValues);
  BOOST_TEST(outValues.sum() == expectedSum);
  BOOST_TEST(outValues(0, 6) == 0.063700286667, boost::test_tools::tolerance(1e-7));
  BOOST_TEST(outValues(0, 23) == 9.2744883594, boost::test_tools::tolerance(1e-7));

  // next go
  cache.resetData();
  outValues.setZero();
  expectedSum = 0;

  // We can use different shapes for different calls
  coords.col(0) = Eigen::RowVector3d(0, 0, 0);
  coords.col(1) = Eigen::RowVector3d(3.5, 0.5, 0.5);
  coords.col(2) = Eigen::RowVector3d(2.5, 0.5, 1.0);
  inData(0, 0)  = 0;
  inData(0, 1)  = 50;
  inData(0, 2)  = 107;
  mapping.mapConservativeAt(coords, inData, cache, outValues);
  expectedSum += inData.sum();

  coords.resize(dimensions, 2);
  inData.resize(dataComponents, 2);
  coords.col(0) = Eigen::RowVector3d(1, 1, 1);
  coords.col(1) = Eigen::RowVector3d(5, 1, 0.5);
  inData(0, 0)  = 108;
  inData(0, 1)  = 48;
  mapping.mapConservativeAt(coords, inData, cache, outValues);
  expectedSum += inData.sum();

  mapping.completeJustInTimeMapping(cache, outValues);
  BOOST_TEST(outValues.sum() == expectedSum);

  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values()(0) == 0.0);
  BOOST_TEST(outData->values()(6) > outData->values()(1));
  BOOST_TEST(outData->values()(13) > outData->values()(9));

  // Now we check for re-computing the mapping after changing the mesh
  mapping.clear();
  outMesh->index().clear();
  outMesh->createVertex(Eigen::Vector3d(6, 0, 0));
  outMesh->allocateDataValues();

  BOOST_TEST(mapping.hasComputedMapping() == false);
  mapping.computeMapping();
  BOOST_TEST(mapping.hasComputedMapping() == true);

  mapping.initializeMappingDataCache(cache);

  Eigen::Map<Eigen::MatrixXd> outValues2(outData->values().data(), dataComponents, outMesh->vertices().size());

  // before writing first, we have to reset the cache to contain zeros
  cache.resetData();
  outValues2.setZero();
  expectedSum = 0;

  coords.resize(dimensions, 1);
  inData.resize(dataComponents, 1);
  coords.col(0) = Eigen::RowVector3d(0.5, 0.5, 0.5);
  inData(0, 0)  = 4.3;
  mapping.mapConservativeAt(coords, inData, cache, outValues2);
  expectedSum += inData.sum();
  coords.col(0) = Eigen::RowVector3d(6.0, 0, 0);
  inData(0, 0)  = 42;
  mapping.mapConservativeAt(coords, inData, cache, outValues2);
  expectedSum += inData.sum();

  mapping.completeJustInTimeMapping(cache, outValues2);
  BOOST_TEST(outValues2.sum() == expectedSum);
  BOOST_TEST(outValues2(0, 24) == 42, boost::test_tools::tolerance(1e-12));
  BOOST_TEST(outValues2(0, 0) == 0.505106848, boost::test_tools::tolerance(1e-7));
}

void perform2DTestJustInTimeMappingConservative(Mapping &mapping)
{
  const int dimensions     = 2;
  const int dataComponents = 2;

  // the dummy from mesh
  mesh::PtrMesh inMesh = mesh::MeshConfiguration::getJustInTimeMappingMesh(dimensions);

  // Create mesh to map from
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData outData = outMesh->createData("OutData", dataComponents, 1_dataID);

  outMesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  outMesh->createVertex(Eigen::Vector2d(1.0, 0.0));
  outMesh->createVertex(Eigen::Vector2d(1.0, 1.0));
  outMesh->createVertex(Eigen::Vector2d(0.0, 1.0));

  outMesh->createVertex(Eigen::Vector2d(2.0, 0.0));
  outMesh->createVertex(Eigen::Vector2d(3.0, 0.0));
  outMesh->createVertex(Eigen::Vector2d(3.0, 1.0));
  outMesh->createVertex(Eigen::Vector2d(2.0, 1.0));

  outMesh->createVertex(Eigen::Vector2d(4.0, 0.0));
  outMesh->createVertex(Eigen::Vector2d(5.0, 0.0));
  outMesh->createVertex(Eigen::Vector2d(5.0, 1.0));
  outMesh->createVertex(Eigen::Vector2d(4.0, 1.0));

  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  mapping.computeMapping();
  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Now, we can setup the MappingDataCache with the structures
  // from computeMapping before starting the mappings
  mapping::impl::MappingDataCache cache(dataComponents);
  mapping.initializeMappingDataCache(cache);

  // Test infrastructure/ the user input
  Eigen::MatrixXd             coords(dimensions, 1);
  Eigen::MatrixXd             inData(dataComponents, 1);
  Eigen::Map<Eigen::MatrixXd> outValues(outData->values().data(), dataComponents, outMesh->vertices().size());

  // before writing first, we have to reset the cache to contain zeros
  cache.resetData();
  outValues.setZero();

  coords.col(0) = Eigen::RowVector2d(5.0, 0.0);
  inData.col(0) = Eigen::RowVector2d(7.3, 14.6);

  // note that the outvalues are only filled in the completeJustInTimeMapping
  mapping.mapConservativeAt(coords, inData, cache, outValues);

  // Once all values are written, we complete the mapping (expensive part)
  mapping.completeJustInTimeMapping(cache, outValues);
  BOOST_TEST(outValues.sum() == inData.sum());
  Eigen::MatrixXd expectedValues(dataComponents, 12);
  expectedValues.setZero();
  expectedValues.col(9) = Eigen::RowVector2d(7.3, 14.6);
  BOOST_TEST(outValues.row(0) == expectedValues.row(0), boost::test_tools::per_element());
  BOOST_TEST(outValues.row(1) == expectedValues.row(1), boost::test_tools::per_element());
  BOOST_TEST(outValues.row(0) == 0.5 * outValues.row(1), boost::test_tools::per_element());

  // clear data again for the next go
  cache.resetData();
  outValues.setZero();

  // // if we call the function twice for the same coordinate, we get twice the contribution
  mapping.mapConservativeAt(coords, inData, cache, outValues);
  mapping.mapConservativeAt(coords, inData, cache, outValues);
  mapping.completeJustInTimeMapping(cache, outValues);
  BOOST_TEST(outValues.row(0) == 2 * expectedValues.row(0), boost::test_tools::per_element());
  BOOST_TEST(outValues.row(1) == 2 * expectedValues.row(1), boost::test_tools::per_element());
  BOOST_TEST(outValues.row(0) == 0.5 * outValues.row(1), boost::test_tools::per_element());

  // // clear data again for the next go
  cache.resetData();
  outValues.setZero();

  Eigen::Vector2d expectedSum;
  expectedSum.setZero();
  mapping.mapConservativeAt(coords, inData, cache, outValues);
  expectedSum += inData.rowwise().sum();

  // We can also change the shape of the user input
  coords.resize(dimensions, 3);
  inData.resize(dataComponents, 3);
  coords.col(0) = Eigen::RowVector2d(3.5, 0.5);
  coords.col(1) = Eigen::RowVector2d(4.5, 1.0);
  coords.col(2) = Eigen::RowVector2d(0.1, 0.0);

  inData.col(0) = Eigen::RowVector2d(4.3, 8.6);
  inData.col(1) = Eigen::RowVector2d(17.3, 34.6);
  inData.col(2) = Eigen::RowVector2d(5.8, 11.6);

  mapping.mapConservativeAt(coords, inData, cache, outValues);
  expectedSum += inData.rowwise().sum();

  mapping.completeJustInTimeMapping(cache, outValues);
  // We don't reach full conservative'ness with the RBF config
  BOOST_TEST(outValues.rowwise().sum()(0) == expectedSum(0), boost::test_tools::tolerance(0.6));
  BOOST_TEST(outValues.rowwise().sum()(1) == expectedSum(1), boost::test_tools::tolerance(0.6));
  BOOST_TEST(outValues.row(0) == 0.5 * outValues.row(1), boost::test_tools::per_element());
  BOOST_TEST(outValues(0, 10) == 8.76997588, boost::test_tools::tolerance(1e-7));

  // Now we check for re-computing the mapping after changing the mesh
  mapping.clear();
  outMesh->index().clear();
  outMesh->createVertex(Eigen::Vector2d(6.0, 0.0));
  outMesh->allocateDataValues();

  BOOST_TEST(mapping.hasComputedMapping() == false);
  mapping.computeMapping();
  BOOST_TEST(mapping.hasComputedMapping() == true);

  mapping.initializeMappingDataCache(cache);

  Eigen::Map<Eigen::MatrixXd> outValues2(outData->values().data(), dataComponents, outMesh->vertices().size());

  // before writing first, we have to reset the cache to contain zeros
  cache.resetData();
  outValues2.setZero();
  expectedSum.setZero();

  coords.resize(dimensions, 1);
  inData.resize(dataComponents, 1);
  coords.col(0) = Eigen::RowVector2d(0.5, 0.5);
  inData.col(0) = Eigen::RowVector2d(4.3, 43);
  mapping.mapConservativeAt(coords, inData, cache, outValues2);
  expectedSum += inData.rowwise().sum();
  coords.col(0) = Eigen::RowVector2d(6.0, 0);
  inData.col(0) = Eigen::RowVector2d(5.1, 51);
  mapping.mapConservativeAt(coords, inData, cache, outValues2);
  expectedSum += inData.rowwise().sum();

  mapping.completeJustInTimeMapping(cache, outValues2);
  BOOST_TEST(outValues2.rowwise().sum()(0) == expectedSum(0), boost::test_tools::tolerance(0.6));
  BOOST_TEST(outValues2.rowwise().sum()(1) == expectedSum(1), boost::test_tools::tolerance(0.6));
  BOOST_TEST(outValues2.col(12) == Eigen::RowVector2d(5.1, 51), boost::test_tools::per_element());
  BOOST_TEST(outValues2(0, 0) == 1.13162327, boost::test_tools::tolerance(1e-7));
}

void perform3DTestConservativeMappingVector(Mapping &mapping)
{
  const int dimensions    = 3;
  const int dataDimension = dimensions;
  using Eigen::Vector3d;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", dataDimension, 0_dataID);
  const DataID  inDataID = inData->getID();
  mesh::Vertex &vertex   = inMesh->createVertex(Vector3d(0, 0, 0));
  mesh::Vertex &vertex1  = inMesh->createVertex(Vector3d(3.5, 0.5, 0.5));
  mesh::Vertex &vertex2  = inMesh->createVertex(Vector3d(2.5, 0.5, 1.0));
  inMesh->createVertex(Vector3d(1, 1, 1));
  inMesh->createVertex(Vector3d(5, 1, 0.5));
  inMesh->allocateDataValues();
  addGlobalIndex(inMesh);

  // Create mesh to map from
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", dataDimension, 1_dataID);
  const DataID  outDataID = outData->getID();
  outMesh->createVertex(Vector3d(0.0, 0.0, 0.0));
  outMesh->createVertex(Vector3d(1.0, 0.0, 0.0));
  outMesh->createVertex(Vector3d(1.0, 1.0, 0.0));
  outMesh->createVertex(Vector3d(0.0, 1.0, 0.0));

  outMesh->createVertex(Vector3d(0.0, 0.0, 1.0));
  outMesh->createVertex(Vector3d(1.0, 0.0, 1.0));
  outMesh->createVertex(Vector3d(1.0, 1.0, 1.0));
  outMesh->createVertex(Vector3d(0.0, 1.0, 1.0));

  outMesh->createVertex(Vector3d(2.0, 0.0, 0.0));
  outMesh->createVertex(Vector3d(3.0, 0.0, 0.0));
  outMesh->createVertex(Vector3d(3.0, 1.0, 0.0));
  outMesh->createVertex(Vector3d(2.0, 1.0, 0.0));

  outMesh->createVertex(Vector3d(2.0, 0.0, 1.0));
  outMesh->createVertex(Vector3d(3.0, 0.0, 1.0));
  outMesh->createVertex(Vector3d(3.0, 1.0, 1.0));
  outMesh->createVertex(Vector3d(2.0, 1.0, 1.0));

  outMesh->createVertex(Vector3d(4.0, 0.0, 0.0));
  outMesh->createVertex(Vector3d(5.0, 0.0, 0.0));
  outMesh->createVertex(Vector3d(5.0, 1.0, 0.0));
  outMesh->createVertex(Vector3d(4.0, 1.0, 0.0));

  outMesh->createVertex(Vector3d(4.0, 0.0, 1.0));
  outMesh->createVertex(Vector3d(5.0, 0.0, 1.0));
  outMesh->createVertex(Vector3d(5.0, 1.0, 1.0));
  outMesh->createVertex(Vector3d(4.0, 1.0, 1.0));

  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Check for the correct sum in the first component
  for (unsigned int i = 0; i < inData->values().size() / dataDimension; ++i) {
    inData->values()(i * dataDimension)     = 7;
    inData->values()(i * dataDimension + 1) = 0;
    inData->values()(i * dataDimension + 2) = 0;
  }

  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values().sum() == inData->values().sum());
  BOOST_TEST(outData->values().sum() == 35);
  for (unsigned int i = 0; i < outData->values().size() / dataDimension; ++i) {
    BOOST_TEST(outData->values()(i * dataDimension + 1) == 0);
    BOOST_TEST(outData->values()(i * dataDimension + 2) == 0);
  }

  // Check for the correct sum in the second component
  for (unsigned int i = 0; i < inData->values().size() / dataDimension; ++i) {
    inData->values()(i * dataDimension)     = 0;
    inData->values()(i * dataDimension + 1) = 27;
    inData->values()(i * dataDimension + 2) = 0;
  }
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values().sum() == inData->values().sum());
  BOOST_TEST(outData->values().sum() == 135);

  for (unsigned int i = 0; i < outData->values().size() / dataDimension; ++i) {
    BOOST_TEST(outData->values()(i * dataDimension) == 0);
    BOOST_TEST(outData->values()(i * dataDimension + 2) == 0);
  }

  // Check for the correct sum in the third component
  for (unsigned int i = 0; i < inData->values().size() / dataDimension; ++i) {
    inData->values()(i * dataDimension)     = 0;
    inData->values()(i * dataDimension + 1) = 0;
    inData->values()(i * dataDimension + 2) = 3;
  }
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values().sum() == inData->values().sum());
  BOOST_TEST(outData->values().sum() == 15);

  for (unsigned int i = 0; i < outData->values().size() / dataDimension; ++i) {
    BOOST_TEST(outData->values()(i * dataDimension) == 0);
    BOOST_TEST(outData->values()(i * dataDimension + 1) == 0);
  }

  // Check for the correct relation between copmonents
  for (unsigned int i = 0; i < inData->values().size() / dataDimension; ++i) {
    inData->values()(i * dataDimension)     = std::pow(i * dataDimension, 3);
    inData->values()(i * dataDimension + 1) = 5 * std::pow(i * dataDimension, 3);
    inData->values()(i * dataDimension + 2) = 10 * std::pow(i * dataDimension, 3);
  }
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values().sum() == inData->values().sum());

  for (unsigned int i = 0; i < outData->values().size() / dataDimension; ++i) {
    // Some result values are close to zero and we cannot compare the relation between these values
    if (outData->values()(i * dataDimension) > 1e-10) {
      BOOST_TEST(outData->values()(i * dataDimension + 1) == 5 * outData->values()(i * dataDimension));
      BOOST_TEST(outData->values()(i * dataDimension + 2) == 10 * outData->values()(i * dataDimension));
    }
  }
  BOOST_TEST(outData->values()(55) == 4087.8100404079933);

  // The remaining parts should already be covered by the other 3D/2D tests
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(PartitionOfUnityMappingTests)
{
  PRECICE_TEST();
  mapping::CompactPolynomialC0                          function(3);
  mapping::PartitionOfUnityMapping<CompactPolynomialC0> consistentMap2D(Mapping::CONSISTENT, 2, function, Polynomial::SEPARATE, 5, 0.4, false);
  perform2DTestConsistentMapping(consistentMap2D);
  mapping::PartitionOfUnityMapping<CompactPolynomialC0> consistentMap2DVector(Mapping::CONSISTENT, 2, function, Polynomial::SEPARATE, 5, 0.4, false);
  perform2DTestConsistentMappingVector(consistentMap2DVector);
  mapping::PartitionOfUnityMapping<CompactPolynomialC6> consistentMap2DDeadAxis(Mapping::CONSISTENT, 2, mapping::CompactPolynomialC6(6), Polynomial::SEPARATE, 5, 0.4, false);
  performTestConsistentMapDeadAxis(consistentMap2DDeadAxis, 2);
  mapping::PartitionOfUnityMapping<CompactPolynomialC6> consistentMap3DDeadAxis(Mapping::CONSISTENT, 3, mapping::CompactPolynomialC6(8), Polynomial::SEPARATE, 5, 0.265, false);
  performTestConsistentMapDeadAxis(consistentMap3DDeadAxis, 3);
  mapping::PartitionOfUnityMapping<CompactPolynomialC0> conservativeMap2D(Mapping::CONSERVATIVE, 2, function, Polynomial::SEPARATE, 5, 0.4, false);
  perform2DTestConservativeMapping(conservativeMap2D);
  mapping::PartitionOfUnityMapping<CompactPolynomialC0> conservativeMap2DVector(Mapping::CONSERVATIVE, 2, function, Polynomial::SEPARATE, 5, 0.4, false);
  perform2DTestConservativeMappingVector(conservativeMap2DVector);
  mapping::PartitionOfUnityMapping<CompactPolynomialC0> consistentMap3D(Mapping::CONSISTENT, 3, function, Polynomial::SEPARATE, 5, 0.265, false);
  perform3DTestConsistentMapping(consistentMap3D);
  mapping::PartitionOfUnityMapping<CompactPolynomialC0> consistentMap3DVector(Mapping::CONSISTENT, 3, function, Polynomial::SEPARATE, 5, 0.265, false);
  perform3DTestConsistentMappingVector(consistentMap3DVector);
  mapping::PartitionOfUnityMapping<CompactPolynomialC0> conservativeMap3D(Mapping::CONSERVATIVE, 3, function, Polynomial::SEPARATE, 5, 0.265, false);
  perform3DTestConservativeMapping(conservativeMap3D);
  mapping::PartitionOfUnityMapping<CompactPolynomialC0> conservativeMap3DVector(Mapping::CONSERVATIVE, 3, function, Polynomial::SEPARATE, 5, 0.265, false);
  perform3DTestConservativeMappingVector(conservativeMap3DVector);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(JustInTimeMapping)
{
  PRECICE_TEST();
  // using scalar data
  mapping::CompactPolynomialC0                          c0function(3);
  mapping::PartitionOfUnityMapping<CompactPolynomialC0> polynomial3Dconsistent(Mapping::CONSISTENT, 3, c0function, Polynomial::SEPARATE, 5, 0.265, false);
  perform3DTestJustInTimeMappingWithPolynomial(polynomial3Dconsistent);
  // using vector data
  mapping::CompactPolynomialC4                          c4function(3);
  mapping::PartitionOfUnityMapping<CompactPolynomialC4> noPolynomial2Dconsistent(Mapping::CONSISTENT, 2, c4function, Polynomial::OFF, 5, 0.265, false);
  perform2DTestJustInTimeMappingNoPolynomial(noPolynomial2Dconsistent);

  // using scalar data
  mapping::CompactPolynomialC2                          c2function(3);
  mapping::PartitionOfUnityMapping<CompactPolynomialC2> polynomial3Dconservative(Mapping::CONSERVATIVE, 3, c2function, Polynomial::SEPARATE, 5, 0.265, false);
  perform3DTestJustInTimeMappingConservative(polynomial3Dconservative);
  // using vector data
  mapping::CompactPolynomialC6                          c6function(10);
  mapping::PartitionOfUnityMapping<CompactPolynomialC6> noPolynomial2Dconservative(Mapping::CONSERVATIVE, 2, c6function, Polynomial::OFF, 5, 0.265, false);
  perform2DTestJustInTimeMappingConservative(noPolynomial2Dconservative);
}

// Test for small meshes, where the number of requested vertices per cluster is bigger than the global
PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TestSingleClusterPartitionOfUnity)
{
  PRECICE_TEST();
  mapping::CompactPolynomialC0                          function(3);
  mapping::PartitionOfUnityMapping<CompactPolynomialC0> mapping(Mapping::CONSISTENT, 3, function, Polynomial::SEPARATE, 50, 0.4, false);

  int dimensions = 3;
  using Eigen::Vector3d;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", 1, 0_dataID);
  int           inDataID = inData->getID();
  inMesh->createVertex(Vector3d(0.0, 0.0, 0.0));
  inMesh->createVertex(Vector3d(1.0, 0.0, 0.0));
  inMesh->createVertex(Vector3d(1.0, 1.0, 0.0));
  inMesh->createVertex(Vector3d(0.0, 1.0, 0.0));

  inMesh->createVertex(Vector3d(0.0, 0.0, 1.0));
  inMesh->createVertex(Vector3d(1.0, 0.0, 1.0));
  inMesh->createVertex(Vector3d(1.0, 1.0, 1.0));
  inMesh->createVertex(Vector3d(0.0, 1.0, 1.0));

  inMesh->createVertex(Vector3d(2.0, 0.0, 0.0));
  inMesh->createVertex(Vector3d(3.0, 0.0, 0.0));
  inMesh->createVertex(Vector3d(3.0, 1.0, 0.0));
  inMesh->createVertex(Vector3d(2.0, 1.0, 0.0));

  inMesh->createVertex(Vector3d(2.0, 0.0, 1.0));
  inMesh->createVertex(Vector3d(3.0, 0.0, 1.0));
  inMesh->createVertex(Vector3d(3.0, 1.0, 1.0));
  inMesh->createVertex(Vector3d(2.0, 1.0, 1.0));

  inMesh->createVertex(Vector3d(4.0, 0.0, 0.0));
  inMesh->createVertex(Vector3d(5.0, 0.0, 0.0));
  inMesh->createVertex(Vector3d(5.0, 1.0, 0.0));
  inMesh->createVertex(Vector3d(4.0, 1.0, 0.0));

  inMesh->createVertex(Vector3d(4.0, 0.0, 1.0));
  inMesh->createVertex(Vector3d(5.0, 0.0, 1.0));
  inMesh->createVertex(Vector3d(5.0, 1.0, 1.0));
  inMesh->createVertex(Vector3d(4.0, 1.0, 1.0));

  inMesh->allocateDataValues();
  addGlobalIndex(inMesh);

  auto &values = inData->values();
  // Set the values in the parallel (z) plane 3*values
  values << 1.0, 2.0, 2.0, 1.0, 3.0, 6.0, 6.0, 3.0, 3.0, 4.0, 4.0, 3.0, 9.0, 12.0, 12.0, 9.0, 5.0, 6.0, 6.0, 5.0, 15.0, 18.0, 18.0, 15.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", 1, 1_dataID);
  int           outDataID = outData->getID();
  mesh::Vertex &vertex    = outMesh->createVertex(Vector3d(0, 0, 0));
  mesh::Vertex &vertex1   = outMesh->createVertex(Vector3d(3.5, 0.5, 0.5));
  mesh::Vertex &vertex2   = outMesh->createVertex(Vector3d(2.5, 0.5, 1.0));
  outMesh->createVertex(Vector3d(0, 0, 0.5));
  outMesh->createVertex(Vector3d(5, 1, 0.5));
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  vertex.setCoords(Vector3d(0.0, 0.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  double value  = outData->values()(0);
  double value1 = outData->values()(1);
  double value2 = outData->values()(2);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.0);
  BOOST_TEST(value1 == 9.0);
  BOOST_TEST(value2 == 10.5);

  vertex.setCoords(Vector3d(0.0, 0.5, 0.5));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Vector3d(1.0, 0.0, 0.0));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Vector3d(1.0, 0.5, 0.5));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 4.0);

  vertex.setCoords(Vector3d(0.5, 0.5, 1.0));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value > 4.6);

  vertex.setCoords(Vector3d(0.5, 1.0, 0.0));
  mapping.clear();
  outData->values().setZero();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value < 1.4);
}

BOOST_AUTO_TEST_SUITE_END() // Serial

BOOST_AUTO_TEST_SUITE(Parallel)

void addGlobalIndex(mesh::PtrMesh &mesh, int offset = 0)
{
  for (mesh::Vertex &v : mesh->vertices()) {
    v.setGlobalIndex(v.getID() + offset);
  }
}

/// Holds rank, owner, position and value of a single vertex
struct VertexSpecification {
  int                 rank;
  int                 owner;
  std::vector<double> position;
  std::vector<double> value;
};

/*
MeshSpecification format:
{ {rank, owner rank, {x, y, z}, {v}}, ... }

also see struct VertexSpecification.

- -1 on rank means all ranks
- -1 on owner rank means no rank
- x, y, z is position of vertex, z is optional, 2D mesh will be created then
- v is the value of the respective vertex. Only 1D supported at this time.

ReferenceSpecification format:
{ {rank, {v}, ... }
- -1 on rank means all ranks
- v is the expected value of n-th vertex on that particular rank
*/
using MeshSpecification = std::vector<VertexSpecification>;

/// Contains which values are expected on which rank: rank -> vector of data.
using ReferenceSpecification = std::vector<std::pair<int, std::vector<double>>>;

void getDistributedMesh(const TestContext &      context,
                        MeshSpecification const &vertices,
                        mesh::PtrMesh &          mesh,
                        mesh::PtrData &          data,
                        int                      globalIndexOffset = 0,
                        bool                     meshIsSmaller     = false)
{
  Eigen::VectorXd d;

  int i = 0;
  for (auto &vertex : vertices) {
    if (vertex.rank == context.rank or vertex.rank == -1) {
      if (vertex.position.size() == 3) // 3-dimensional
        mesh->createVertex(Eigen::Vector3d(vertex.position.data()));
      else if (vertex.position.size() == 2) // 2-dimensional
        mesh->createVertex(Eigen::Vector2d(vertex.position.data()));

      int valueDimension = vertex.value.size();

      if (vertex.owner == context.rank)
        mesh->vertices().back().setOwner(true);
      else
        mesh->vertices().back().setOwner(false);

      d.conservativeResize(i * valueDimension + valueDimension);
      // Get data in every dimension
      for (int dim = 0; dim < valueDimension; ++dim) {
        d(i * valueDimension + dim) = vertex.value.at(dim);
      }
      i++;
    }
  }
  addGlobalIndex(mesh, globalIndexOffset);
  mesh->allocateDataValues();
  // All tests use eight vertices
  // if (meshIsSmaller) {
  //   mesh->setGlobalNumberOfVertices(7);
  // } else {
  //   mesh->setGlobalNumberOfVertices(8);
  // }
  data->values() = d;
}

void testDistributed(const TestContext &    context,
                     Mapping &              mapping,
                     MeshSpecification      inMeshSpec,
                     MeshSpecification      outMeshSpec,
                     ReferenceSpecification referenceSpec,
                     int                    inGlobalIndexOffset = 0,
                     bool                   meshIsSmaller       = false)
{
  int meshDimension  = inMeshSpec.at(0).position.size();
  int valueDimension = inMeshSpec.at(0).value.size();

  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", meshDimension, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", valueDimension, 0_dataID);
  int           inDataID = inData->getID();

  getDistributedMesh(context, inMeshSpec, inMesh, inData, inGlobalIndexOffset);

  mesh::PtrMesh outMesh(new mesh::Mesh("outMesh", meshDimension, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", valueDimension, 1_dataID);
  int           outDataID = outData->getID();

  getDistributedMesh(context, outMeshSpec, outMesh, outData, 0, meshIsSmaller);

  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  mapping.computeMapping();
  BOOST_TEST(mapping.hasComputedMapping() == true);
  mapping.map(inDataID, outDataID);

  int index = 0;
  for (auto &referenceVertex : referenceSpec) {
    if (referenceVertex.first == context.rank or referenceVertex.first == -1) {
      for (int dim = 0; dim < valueDimension; ++dim) {
        BOOST_TEST_INFO("Index of vertex: " << index << " - Dimension: " << dim);
        BOOST_TEST(outData->values()(index * valueDimension + dim) == referenceVertex.second.at(dim));
      }
      ++index;
    }
  }
  BOOST_TEST(outData->values().size() == index * valueDimension);
}

PRECICE_TEST_SETUP(""_on(4_ranks).setupIntraComm())
BOOST_AUTO_TEST_CASE(DistributedConsistent2D)
{
  PRECICE_TEST();
  std::vector<int> globalIndexOffsets = {0, 0, 0, 0};

  MeshSpecification in{// Consistent mapping: The inMesh is communicated
                       {-1, 0, {0, 0}, {1}},
                       {-1, 0, {0, 1}, {2}},
                       {-1, 1, {1, 0}, {3}},
                       {-1, 1, {1, 1}, {4}},
                       {-1, 2, {2, 0}, {5}},
                       {-1, 2, {2, 1}, {6}},
                       {-1, 3, {3, 0}, {7}},
                       {-1, 3, {3, 1}, {8}}};
  MeshSpecification out{// The outMesh is local, distributed among all ranks
                        {0, -1, {0, 0}, {0}},
                        {0, -1, {0, 1}, {0}},
                        {1, -1, {1, 0}, {0}},
                        {1, -1, {1, 1}, {0}},
                        {2, -1, {2, 0}, {0}},
                        {2, -1, {2, 1}, {0}},
                        {3, -1, {3, 0}, {0}},
                        {3, -1, {3, 1}, {0}}};

  ReferenceSpecification ref{// Tests for {0, 1} on the first rank, {1, 2} on the second, ...
                             {0, {1}},
                             {0, {2}},
                             {1, {3}},
                             {1, {4}},
                             {2, {5}},
                             {2, {6}},
                             {3, {7}},
                             {3, {8}}};

  mapping::PartitionOfUnityMapping<CompactPolynomialC6> consistentMap2D(Mapping::CONSISTENT, 2, CompactPolynomialC6(3.), Polynomial::SEPARATE, 5, 0.3, false);
  testDistributed(context, consistentMap2D, in, out, ref, globalIndexOffsets.at(context.rank));
}

// Same as above, but including empty ranks not participating
PRECICE_TEST_SETUP(""_on(4_ranks).setupIntraComm())
BOOST_AUTO_TEST_CASE(DistributedConsistent2DEmptyOut)
{
  PRECICE_TEST();
  std::vector<int> globalIndexOffsets = {0, 0, 0, 0};

  MeshSpecification in{// Consistent mapping: The inMesh is communicated
                       {-1, 0, {0, 0}, {1}},
                       {-1, 0, {0, 1}, {2}},
                       {-1, 1, {1, 0}, {3}},
                       {-1, 1, {1, 1}, {4}},
                       {-1, 2, {2, 0}, {5}},
                       {-1, 2, {2, 1}, {6}},
                       {-1, 3, {3, 0}, {7}},
                       {-1, 3, {3, 1}, {8}}};
  MeshSpecification out{// The outMesh is local, distributed among all ranks
                        {0, 0, {0, 0}, {0}},
                        {0, 0, {0, 1}, {0}},
                        {1, 1, {1, 0}, {0}},
                        {1, 1, {1, 1}, {0}},
                        {2, 2, {2, 0}, {0}},
                        {2, 2, {2, 1}, {0}},
                        {2, 2, {3, 0}, {0}},
                        {2, 2, {3, 1}, {0}}};

  ReferenceSpecification ref{// Tests for {0, 1} on the first rank, {1, 2} on the second, ...
                             {0, {1}},
                             {0, {2}},
                             {1, {3}},
                             {1, {4}},
                             {2, {5}},
                             {2, {6}},
                             {2, {7}},
                             {2, {8}}};

  mapping::PartitionOfUnityMapping<CompactPolynomialC6> consistentMap2D(Mapping::CONSISTENT, 2, CompactPolynomialC6(3.), Polynomial::SEPARATE, 5, 0.3, false);
  testDistributed(context, consistentMap2D, in, out, ref, globalIndexOffsets.at(context.rank));
}

// Same as above, but including empty ranks not participating
PRECICE_TEST_SETUP(""_on(4_ranks).setupIntraComm())
BOOST_AUTO_TEST_CASE(DistributedConsistent2DEmptyRank)
{
  PRECICE_TEST();
  std::vector<int> globalIndexOffsets = {0, 0, 0, 0};

  MeshSpecification in{// Consistent mapping: The inMesh is communicated
                       {-1, 0, {0, 0}, {1}},
                       {-1, 0, {0, 1}, {2}},
                       {-1, 1, {1, 0}, {3}},
                       {-1, 1, {1, 1}, {4}},
                       {-1, 2, {2, 0}, {5}},
                       {-1, 2, {2, 1}, {6}}};

  MeshSpecification out{// The outMesh is local, distributed among all ranks
                        {0, 0, {0, 0}, {0}},
                        {0, 0, {0, 1}, {0}},
                        {1, 1, {1, 0}, {0}},
                        {1, 1, {1, 1}, {0}},
                        {2, 2, {2, 0}, {0}},
                        {2, 2, {2, 1}, {0}}};

  ReferenceSpecification ref{// Tests for {0, 1} on the first rank, {1, 2} on the second, ...
                             {0, {1}},
                             {0, {2}},
                             {1, {3}},
                             {1, {4}},
                             {2, {5}},
                             {2, {6}}};

  mapping::PartitionOfUnityMapping<CompactPolynomialC6> consistentMap2D(Mapping::CONSISTENT, 2, CompactPolynomialC6(3.), Polynomial::SEPARATE, 5, 0.3, false);
  testDistributed(context, consistentMap2D, in, out, ref, globalIndexOffsets.at(context.rank));
}

PRECICE_TEST_SETUP(""_on(4_ranks).setupIntraComm())
BOOST_AUTO_TEST_CASE(DistributedConservative2D)
{
  PRECICE_TEST();
  std::vector<int> globalIndexOffsets = {0, 0, 0, 0};

  MeshSpecification in{// Conservative mapping: The inMesh is local
                       {0, -1, {0, 0}, {1}},
                       {0, -1, {0, 1}, {2}},
                       {1, -1, {1, 0}, {3}},
                       {1, -1, {1, 1}, {4}},
                       {2, -1, {2, 0}, {5}},
                       {2, -1, {2, 1}, {6}},
                       {3, -1, {3, 0}, {7}},
                       {3, -1, {3, 1}, {8}}};
  MeshSpecification out{// The outMesh is remote, distributed among all ranks
                        {0, 0, {0, 0}, {0}},
                        {0, 0, {0, 1}, {0}},
                        {1, 1, {1, 0}, {0}},
                        {1, 1, {1, 1}, {0}},
                        {2, 2, {2, 0}, {0}},
                        {2, 2, {2, 1}, {0}},
                        {3, 3, {3, 0}, {0}},
                        {3, 3, {3, 1}, {0}}};

  ReferenceSpecification ref{// Tests for {0, 1} on the first rank, {1, 2} on the second, ...
                             {0, {1}},
                             {0, {2}},
                             {1, {3}},
                             {1, {4}},
                             {2, {5}},
                             {2, {6}},
                             {3, {7}},
                             {3, {8}}};

  mapping::PartitionOfUnityMapping<CompactPolynomialC6> conservativeMap2D(Mapping::CONSERVATIVE, 2, CompactPolynomialC6(3.), Polynomial::SEPARATE, 5, 0.3, false);
  testDistributed(context, conservativeMap2D, in, out, ref, globalIndexOffsets.at(context.rank));
}

// Same as above, but checking with empty ranks
PRECICE_TEST_SETUP(""_on(4_ranks).setupIntraComm())
BOOST_AUTO_TEST_CASE(DistributedConservative2DEmptyRank)
{
  PRECICE_TEST();
  std::vector<int> globalIndexOffsets = {0, 2, 4, 6};

  MeshSpecification in{// Conservative mapping: The inMesh is local
                       {0, -1, {0, 0}, {1}},
                       {0, -1, {0, 1}, {2}},
                       {1, -1, {1, 0}, {3}},
                       {1, -1, {1, 1}, {4}},
                       {2, -1, {2, 0}, {5}},
                       {2, -1, {2, 1}, {6}},
                       {2, -1, {3, 0}, {7}},
                       {2, -1, {3, 1}, {8}}};
  MeshSpecification out{// The outMesh is remote, distributed among all ranks
                        {-1, 0, {0, 0}, {0}},
                        {-1, 0, {0, 1}, {0}},
                        {-1, 1, {1, 0}, {0}},
                        {-1, 1, {1, 1}, {0}},
                        {-1, 2, {2, 0}, {0}},
                        {-1, 2, {2, 1}, {0}},
                        {-1, 3, {3, 0}, {0}},
                        {-1, 3, {3, 1}, {0}}};

  ReferenceSpecification ref{// Tests for {0, 1} on the first rank, {1, 2} on the second, ...
                             {0, {1}},
                             {0, {2}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {1, {0}},
                             {1, {0}},
                             {1, {3}},
                             {1, {4}},
                             {1, {0}},
                             {1, {0}},
                             {1, {0}},
                             {1, {0}},
                             {2, {0}},
                             {2, {0}},
                             {2, {0}},
                             {2, {0}},
                             {2, {5}},
                             {2, {6}},
                             {2, {7}},
                             {2, {8}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}}};

  mapping::PartitionOfUnityMapping<CompactPolynomialC6> conservativeMap2D(Mapping::CONSERVATIVE, 2, CompactPolynomialC6(3.), Polynomial::SEPARATE, 5, 0.3, false);
  testDistributed(context, conservativeMap2D, in, out, ref, globalIndexOffsets.at(context.rank));
}

// Same as above, but checking the primary rank for all output vertices
PRECICE_TEST_SETUP(""_on(2_ranks).setupIntraComm())
BOOST_AUTO_TEST_CASE(DistributedConservative2DTwoRanks)
{
  PRECICE_TEST();
  std::vector<int> globalIndexOffsets = {0, 2};

  MeshSpecification in{// Conservative mapping: The inMesh is local
                       {0, -1, {0, 0}, {1}},
                       {0, -1, {0, 1}, {2}},
                       {0, -1, {1, 0}, {3}},
                       {0, -1, {1, 1}, {4}},
                       {1, -1, {2, 0}, {5}},
                       {1, -1, {2, 1}, {6}},
                       {1, -1, {3, 0}, {7}},
                       {1, -1, {3, 1}, {8}}};
  MeshSpecification out{// The outMesh is remote, distributed among all ranks
                        {-1, 0, {0, 0}, {0}},
                        {-1, 0, {0, 1}, {0}},
                        {-1, 0, {1, 0}, {0}},
                        {-1, 0, {1, 1}, {0}},
                        {-1, 0, {2, 0}, {0}},
                        {-1, 0, {2, 1}, {0}},
                        {-1, 0, {3, 0}, {0}},
                        {-1, 0, {3, 1}, {0}}};

  ReferenceSpecification ref{// Tests for {0, 1} on the first rank, {1, 2} on the second, ...
                             {0, {1}},
                             {0, {2}},
                             {0, {3}},
                             {0, {4}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {1, {0}},
                             {1, {0}},
                             {1, {0}},
                             {1, {0}},
                             {1, {5}},
                             {1, {6}},
                             {1, {7}},
                             {1, {8}}};

  mapping::PartitionOfUnityMapping<CompactPolynomialC6> conservativeMap2D(Mapping::CONSERVATIVE, 2, CompactPolynomialC6(3.), Polynomial::SEPARATE, 5, 0.3, false);
  testDistributed(context, conservativeMap2D, in, out, ref, globalIndexOffsets.at(context.rank));
}

void testTagging(const TestContext &context,
                 MeshSpecification  inMeshSpec,
                 MeshSpecification  outMeshSpec,
                 MeshSpecification  shouldTagFirstRound,
                 MeshSpecification  shouldTagSecondRound,
                 bool               consistent)
{
  int meshDimension  = inMeshSpec.at(0).position.size();
  int valueDimension = inMeshSpec.at(0).value.size();

  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", meshDimension, testing::nextMeshID()));
  mesh::PtrData inData = inMesh->createData("InData", valueDimension, 0_dataID);
  getDistributedMesh(context, inMeshSpec, inMesh, inData);

  mesh::PtrMesh outMesh(new mesh::Mesh("outMesh", meshDimension, testing::nextMeshID()));
  mesh::PtrData outData = outMesh->createData("OutData", valueDimension, 1_dataID);
  getDistributedMesh(context, outMeshSpec, outMesh, outData);

  Mapping::Constraint                                   constr = consistent ? Mapping::CONSISTENT : Mapping::CONSERVATIVE;
  mapping::PartitionOfUnityMapping<CompactPolynomialC4> mapping(constr, 2, CompactPolynomialC4(2), Polynomial::SEPARATE, 2, 0.3, false);
  inMesh->computeBoundingBox();
  outMesh->computeBoundingBox();

  mapping.setMeshes(inMesh, outMesh);
  mapping.tagMeshFirstRound();

  for (const auto &v : inMesh->vertices()) {
    auto pos   = std::find_if(shouldTagFirstRound.begin(), shouldTagFirstRound.end(),
                            [meshDimension, &v](const VertexSpecification &spec) {
                              return std::equal(spec.position.data(), spec.position.data() + meshDimension, v.getCoords().data());
                            });
    bool found = pos != shouldTagFirstRound.end();
    BOOST_TEST(found >= v.isTagged(),
               "FirstRound: Vertex " << v << " is tagged, but should not be.");
    BOOST_TEST(found <= v.isTagged(),
               "FirstRound: Vertex " << v << " is not tagged, but should be.");
  }

  mapping.tagMeshSecondRound();

  for (const auto &v : inMesh->vertices()) {
    auto posFirst    = std::find_if(shouldTagFirstRound.begin(), shouldTagFirstRound.end(),
                                 [meshDimension, &v](const VertexSpecification &spec) {
                                   return std::equal(spec.position.data(), spec.position.data() + meshDimension, v.getCoords().data());
                                 });
    bool foundFirst  = posFirst != shouldTagFirstRound.end();
    auto posSecond   = std::find_if(shouldTagSecondRound.begin(), shouldTagSecondRound.end(),
                                  [meshDimension, &v](const VertexSpecification &spec) {
                                    return std::equal(spec.position.data(), spec.position.data() + meshDimension, v.getCoords().data());
                                  });
    bool foundSecond = posSecond != shouldTagSecondRound.end();
    BOOST_TEST(foundFirst <= v.isTagged(), "SecondRound: Vertex " << v
                                                                  << " is not tagged, but should be from the first round.");
    BOOST_TEST(foundSecond <= v.isTagged(), "SecondRound: Vertex " << v
                                                                   << " is not tagged, but should be.");
    BOOST_TEST((foundSecond or foundFirst) >= v.isTagged(), "SecondRound: Vertex " << v
                                                                                   << " is tagged, but should not be.");
  }
}

PRECICE_TEST_SETUP(""_on(4_ranks).setupIntraComm())
BOOST_AUTO_TEST_CASE(testTagFirstRound)
{
  PRECICE_TEST();

  MeshSpecification outMeshSpec = {
      {0, -1, {0, 0}, {0}}};
  MeshSpecification inMeshSpec = {
      {0, -1, {-1, 0}, {1}}, // inside
      {0, -1, {-3, 0}, {1}}, // outside
      {0, 0, {1, 0}, {1}},   // inside, owner
      {0, -1, {3, 0}, {1}},  // outside
      {0, -1, {0, -1}, {1}}, // inside
      {0, -1, {0, -3}, {1}}, // outside
      {0, -1, {0, 1}, {1}},  // inside
      {0, -1, {0, 3}, {1}}   // outside
  };
  MeshSpecification shouldTagFirstRound = {
      {0, -1, {-1, 0}, {1}},
      {0, -1, {1, 0}, {1}},
      {0, -1, {0, -1}, {1}},
      {0, -1, {0, 1}, {1}}};
  // No tagging happening in round two
  MeshSpecification shouldTagSecondRound;

  testTagging(context, inMeshSpec, outMeshSpec, shouldTagFirstRound, shouldTagSecondRound, true);
  // For conservative just swap meshes.
  testTagging(context, outMeshSpec, inMeshSpec, shouldTagFirstRound, shouldTagSecondRound, false);
}

BOOST_AUTO_TEST_SUITE_END() // Parallel

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

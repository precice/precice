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
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Utils.hpp"
#include "mesh/Vertex.hpp"
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
  BOOST_TEST(mapping.hasComputedMapping() == false);
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Vector2d(0.0, 1.0));
  mapping.clear();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Vector2d(1.0, 0.0));
  mapping.clear();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Vector2d(1.0, 0.5));
  mapping.clear();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Vector2d(1.0, 1.0));
  mapping.clear();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Vector2d(0.5, 0.0));
  mapping.clear();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.5);

  vertex.setCoords(Vector2d(0.5, 0.5));
  mapping.clear();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.5);

  vertex.setCoords(Vector2d(0.5, 1.0));
  mapping.clear();
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
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values()(0) == 1.0);
  BOOST_TEST(outData->values()(1) == 2.0);

  vertex.setCoords(Vector2d(0.0, 1.0));
  mapping.clear();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values()(0) == 1.0);
  BOOST_TEST(outData->values()(1) == 2.0);

  vertex.setCoords(Vector2d(1.0, 0.0));
  mapping.clear();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values()(0) == 2.0);
  BOOST_TEST(outData->values()(1) == 4.0);

  vertex.setCoords(Vector2d(1.0, 0.5));
  mapping.clear();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values()(0) == 2.0);
  BOOST_TEST(outData->values()(1) == 4.0);

  vertex.setCoords(Vector2d(1.0, 1.0));
  mapping.clear();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values()(0) == 2.0);
  BOOST_TEST(outData->values()(1) == 4.0);

  vertex.setCoords(Vector2d(0.5, 0.0));
  mapping.clear();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values()(0) == 1.5);
  BOOST_TEST(outData->values()(1) == 3.0);

  vertex.setCoords(Vector2d(0.5, 0.5));
  mapping.clear();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outData->values()(0) == 1.5);
  BOOST_TEST(outData->values()(1) == 3.0);

  vertex.setCoords(Vector2d(0.5, 1.0));
  mapping.clear();
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
    BOOST_TEST(value == 19935.268150244759);
    BOOST_TEST(value1 == 19935.266962237896);
    BOOST_TEST(value2 == 19935.276990769267);
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
  inMesh->setGlobalNumberOfVertices(inMesh->vertices().size());

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
  outMesh->setGlobalNumberOfVertices(outMesh->vertices().size());

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
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(outData->values().sum() == inData->values().sum());
  BOOST_TEST(outData->values()(5) == outData->values()(8));

  // Test the conservation property if we have everywhere non-zero input data
  values << 5.0, 10.0, 7.0, 3.0, 4.0;

  mapping.clear();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);

  BOOST_TEST(outData->values().sum() == inData->values().sum());
  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Test the conservation property if we have everywhere non-zero input data
  values << 3.0, 4.0, 5.0, 7.0, 9.0;

  mapping.clear();
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
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Vector3d(0.0, 1.0, 1.0));
  mapping.clear();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 3.0);

  vertex.setCoords(Vector3d(1.0, 0.0, 0.0));
  mapping.clear();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Vector3d(1.0, 0.5, 0.5));
  mapping.clear();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 4.0);

  vertex.setCoords(Vector3d(1.0, 1.0, 1.0));
  mapping.clear();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 6.0);

  vertex.setCoords(Vector3d(0.5, 0.0, 0.5));
  mapping.clear();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 3.0);

  vertex.setCoords(Vector3d(0.5, 0.5, 1.0));
  mapping.clear();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 4.5);

  vertex.setCoords(Vector3d(0.5, 1.0, 0.0));
  mapping.clear();
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.5);
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

BOOST_AUTO_TEST_CASE(PartitionOfUnityMappingTests)
{
  PRECICE_TEST(1_rank);
  mapping::CompactPolynomialC0                          function(3);
  std::array<bool, 3>                                   deadAxis({{false, false, false}});
  mapping::PartitionOfUnityMapping<CompactPolynomialC0> consistentMap2D(Mapping::CONSISTENT, 2, function, deadAxis, Polynomial::SEPARATE, 5, 0.4, false);
  perform2DTestConsistentMapping(consistentMap2D);
  mapping::PartitionOfUnityMapping<CompactPolynomialC0> consistentMap2DVector(Mapping::CONSISTENT, 2, function, deadAxis, Polynomial::SEPARATE, 5, 0.4, false);
  perform2DTestConsistentMappingVector(consistentMap2DVector);
  mapping::PartitionOfUnityMapping<CompactPolynomialC6> consistentMap2DDeadAxis(Mapping::CONSISTENT, 2, mapping::CompactPolynomialC6(6), deadAxis, Polynomial::SEPARATE, 5, 0.4, false);
  performTestConsistentMapDeadAxis(consistentMap2DDeadAxis, 2);
  mapping::PartitionOfUnityMapping<CompactPolynomialC6> consistentMap3DDeadAxis(Mapping::CONSISTENT, 3, mapping::CompactPolynomialC6(8), deadAxis, Polynomial::SEPARATE, 5, 0.4, false);
  performTestConsistentMapDeadAxis(consistentMap3DDeadAxis, 3);
  mapping::PartitionOfUnityMapping<CompactPolynomialC0> conservativeMap2D(Mapping::CONSERVATIVE, 2, function, deadAxis, Polynomial::SEPARATE, 5, 0.4, false);
  perform2DTestConservativeMapping(conservativeMap2D);
  mapping::PartitionOfUnityMapping<CompactPolynomialC0> conservativeMap2DVector(Mapping::CONSERVATIVE, 2, function, deadAxis, Polynomial::SEPARATE, 5, 0.4, false);
  perform2DTestConservativeMappingVector(conservativeMap2DVector);
  mapping::PartitionOfUnityMapping<CompactPolynomialC0> consistentMap3D(Mapping::CONSISTENT, 3, function, deadAxis, Polynomial::SEPARATE, 5, 0.4, false);
  perform3DTestConsistentMapping(consistentMap3D);
  mapping::PartitionOfUnityMapping<CompactPolynomialC0> consistentMap3DVector(Mapping::CONSISTENT, 3, function, deadAxis, Polynomial::SEPARATE, 5, 0.4, false);
  perform3DTestConsistentMappingVector(consistentMap3DVector);
  mapping::PartitionOfUnityMapping<CompactPolynomialC0> conservativeMap3D(Mapping::CONSERVATIVE, 3, function, deadAxis, Polynomial::SEPARATE, 5, 0.4, false);
  perform3DTestConservativeMapping(conservativeMap3D);
  mapping::PartitionOfUnityMapping<CompactPolynomialC0> conservativeMap3DVector(Mapping::CONSERVATIVE, 3, function, deadAxis, Polynomial::SEPARATE, 5, 0.4, false);
  perform3DTestConservativeMappingVector(conservativeMap3DVector);
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

BOOST_AUTO_TEST_CASE(DistributedConsistent2D)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
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

  mapping::PartitionOfUnityMapping<CompactPolynomialC6> consistentMap2D(Mapping::CONSISTENT, 2, CompactPolynomialC6(3.), {{false, false, false}}, Polynomial::SEPARATE, 5, 0.3, false);
  testDistributed(context, consistentMap2D, in, out, ref, globalIndexOffsets.at(context.rank));
}

// Same as above, but including empty ranks not participating
BOOST_AUTO_TEST_CASE(DistributedConsistent2DEmptyOut)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
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

  mapping::PartitionOfUnityMapping<CompactPolynomialC6> consistentMap2D(Mapping::CONSISTENT, 2, CompactPolynomialC6(3.), {{false, false, false}}, Polynomial::SEPARATE, 5, 0.3, false);
  testDistributed(context, consistentMap2D, in, out, ref, globalIndexOffsets.at(context.rank));
}

// Same as above, but including empty ranks not participating
BOOST_AUTO_TEST_CASE(DistributedConsistent2DEmptyRank)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
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

  mapping::PartitionOfUnityMapping<CompactPolynomialC6> consistentMap2D(Mapping::CONSISTENT, 2, CompactPolynomialC6(3.), {{false, false, false}}, Polynomial::SEPARATE, 5, 0.3, false);
  testDistributed(context, consistentMap2D, in, out, ref, globalIndexOffsets.at(context.rank));
}

BOOST_AUTO_TEST_CASE(DistributedConservative2D)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
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

  mapping::PartitionOfUnityMapping<CompactPolynomialC6> conservativeMap2D(Mapping::CONSERVATIVE, 2, CompactPolynomialC6(3.), {{false, false, false}}, Polynomial::SEPARATE, 5, 0.3, false);
  testDistributed(context, conservativeMap2D, in, out, ref, globalIndexOffsets.at(context.rank));
}

// Same as above, but checking with empty ranks
BOOST_AUTO_TEST_CASE(DistributedConservative2DEmptyRank)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
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

  mapping::PartitionOfUnityMapping<CompactPolynomialC6> conservativeMap2D(Mapping::CONSERVATIVE, 2, CompactPolynomialC6(3.), {{false, false, false}}, Polynomial::SEPARATE, 5, 0.3, false);
  testDistributed(context, conservativeMap2D, in, out, ref, globalIndexOffsets.at(context.rank));
}

// Same as above, but checking the primary rank for all output vertices
BOOST_AUTO_TEST_CASE(DistributedConservative2DTwoRanks)
{
  PRECICE_TEST(""_on(2_ranks).setupIntraComm());
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

  mapping::PartitionOfUnityMapping<CompactPolynomialC6> conservativeMap2D(Mapping::CONSERVATIVE, 2, CompactPolynomialC6(3.), {{false, false, false}}, Polynomial::SEPARATE, 5, 0.3, false);
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
  mapping::PartitionOfUnityMapping<CompactPolynomialC4> mapping(constr, 2, CompactPolynomialC4(2), {{false, false, false}}, Polynomial::SEPARATE, 2, 0.3, false);
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

BOOST_AUTO_TEST_CASE(testTagFirstRound)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm())

  MeshSpecification outMeshSpec = {
      {0, -1, {0, 0}, {0}}};
  MeshSpecification inMeshSpec = {
      {0, -1, {-1, 0}, {1}}, //inside
      {0, -1, {-3, 0}, {1}}, //outside
      {0, 0, {1, 0}, {1}},   //inside, owner
      {0, -1, {3, 0}, {1}},  //outside
      {0, -1, {0, -1}, {1}}, //inside
      {0, -1, {0, -3}, {1}}, //outside
      {0, -1, {0, 1}, {1}},  //inside
      {0, -1, {0, 3}, {1}}   //outside
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

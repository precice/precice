#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>
#include "logging/Logger.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/RadialBasisFctMapping.hpp"
#include "mapping/config/MappingConfiguration.hpp"
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

// We use internal linkage in order to avoid a separate TU for these functions
namespace {

void addGlobalIndex(mesh::PtrMesh &mesh, int offset = 0)
{
  for (mesh::Vertex &v : mesh->vertices()) {
    v.setGlobalIndex(v.getID() + offset);
  }
}

void testSerialScaledConsistent(mesh::PtrMesh inMesh, mesh::PtrMesh outMesh, mesh::PtrData inData, mesh::PtrData outData)
{
  auto inputIntegral  = mesh::integrateSurface(inMesh, inData);
  auto outputIntegral = mesh::integrateSurface(outMesh, outData);

  for (int dim = 0; dim < inputIntegral.size(); ++dim) {
    BOOST_TEST(inputIntegral(dim) == outputIntegral(dim));
  }
}

// Definition of the helper functions such that we can test different solver setups
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
  inMesh->allocateDataValues();
  addGlobalIndex(inMesh);
  inMesh->setGlobalNumberOfVertices(inMesh->vertices().size());

  auto &values = inData->values();
  values << 1.0, 2.0, 2.0, 1.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", 1, 1_dataID);
  int           outDataID = outData->getID();
  mesh::Vertex &vertex    = outMesh->createVertex(Vector2d(0, 0));
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);
  outMesh->setGlobalNumberOfVertices(outMesh->vertices().size());

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  vertex.setCoords(Vector2d(0.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  double value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Vector2d(0.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Vector2d(0.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Vector2d(1.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Vector2d(1.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Vector2d(1.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Vector2d(0.5, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.5);

  vertex.setCoords(Vector2d(0.5, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.5);

  vertex.setCoords(Vector2d(0.5, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.5);
}

void perform2DTestConsistentMappingVector(Mapping &mapping)
{
  int dimensions = 2;
  using Eigen::Vector2d;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", 2, 0_dataID);
  int           inDataID = inData->getID();
  inMesh->createVertex(Vector2d(0.0, 0.0));
  inMesh->createVertex(Vector2d(1.0, 0.0));
  inMesh->createVertex(Vector2d(1.0, 1.0));
  inMesh->createVertex(Vector2d(0.0, 1.0));
  inMesh->allocateDataValues();
  addGlobalIndex(inMesh);
  inMesh->setGlobalNumberOfVertices(inMesh->vertices().size());

  auto &values = inData->values();
  values << 1.0, 4.0, 2.0, 5.0, 2.0, 5.0, 1.0, 4.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", 2, 1_dataID);
  int           outDataID = outData->getID();
  mesh::Vertex &vertex    = outMesh->createVertex(Vector2d(0, 0));
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);
  outMesh->setGlobalNumberOfVertices(outMesh->vertices().size());

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  vertex.setCoords(Vector2d(0.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  double value1 = outData->values()(0);
  double value2 = outData->values()(1);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value1 == 1.0);
  BOOST_TEST(value2 == 4.0);

  vertex.setCoords(Vector2d(0.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value1 = outData->values()(0);
  value2 = outData->values()(1);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value1 == 1.0);
  BOOST_TEST(value2 == 4.0);

  vertex.setCoords(Vector2d(0.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value1 = outData->values()(0);
  value2 = outData->values()(1);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value1 == 1.0);
  BOOST_TEST(value2 == 4.0);

  vertex.setCoords(Vector2d(1.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value1 = outData->values()(0);
  value2 = outData->values()(1);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value1 == 2.0);
  BOOST_TEST(value2 == 5.0);

  vertex.setCoords(Vector2d(1.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value1 = outData->values()(0);
  value2 = outData->values()(1);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value1 == 2.0);
  BOOST_TEST(value2 == 5.0);

  vertex.setCoords(Vector2d(1.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value1 = outData->values()(0);
  value2 = outData->values()(1);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value1 == 2.0);
  BOOST_TEST(value2 == 5.0);

  vertex.setCoords(Vector2d(0.5, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value1 = outData->values()(0);
  value2 = outData->values()(1);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value1 == 1.5);
  BOOST_TEST(value2 == 4.5);

  vertex.setCoords(Vector2d(0.5, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value1 = outData->values()(0);
  value2 = outData->values()(1);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value1 == 1.5);
  BOOST_TEST(value2 == 4.5);

  vertex.setCoords(Vector2d(0.5, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value1 = outData->values()(0);
  value2 = outData->values()(1);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value1 == 1.5);
  BOOST_TEST(value2 == 4.5);
}

void perform3DTestConsistentMapping(Mapping &mapping)
{
  int dimensions = 3;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", 1, 0_dataID);
  int           inDataID = inData->getID();
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(1.0, 1.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 1.0));
  inMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 1.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 1.0));
  inMesh->createVertex(Eigen::Vector3d(1.0, 1.0, 1.0));
  inMesh->allocateDataValues();
  addGlobalIndex(inMesh);
  inMesh->setGlobalNumberOfVertices(inMesh->vertices().size());

  auto &values = inData->values();
  values << 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", 1, 1_dataID);
  int           outDataID = outData->getID();
  mesh::Vertex &vertex    = outMesh->createVertex(Eigen::Vector3d::Zero());
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);
  outMesh->setGlobalNumberOfVertices(outMesh->vertices().size());

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  vertex.setCoords(Eigen::Vector3d(0.0, 0.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  double value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Eigen::Vector3d(0.0, 0.5, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Eigen::Vector3d(0.5, 0.5, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Eigen::Vector3d(1.0, 0.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Eigen::Vector3d(1.0, 1.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Eigen::Vector3d(0.0, 0.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Eigen::Vector3d(1.0, 0.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Eigen::Vector3d(1.0, 1.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Eigen::Vector3d(0.5, 0.5, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Eigen::Vector3d(0.0, 0.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value, 1.5);

  vertex.setCoords(Eigen::Vector3d(1.0, 0.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.5);

  vertex.setCoords(Eigen::Vector3d(0.0, 1.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.5);

  vertex.setCoords(Eigen::Vector3d(1.0, 1.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.5);

  vertex.setCoords(Eigen::Vector3d(0.5, 0.5, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.5);
}

void perform2DTestScaledConsistentMapping(Mapping &mapping)
{
  int dimensions = 2;
  using Eigen::Vector2d;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", 1, 0_dataID);
  int           inDataID = inData->getID();
  auto &        inV1     = inMesh->createVertex(Vector2d(0.0, 0.0));
  auto &        inV2     = inMesh->createVertex(Vector2d(1.0, 0.0));
  auto &        inV3     = inMesh->createVertex(Vector2d(1.0, 1.0));
  auto &        inV4     = inMesh->createVertex(Vector2d(0.0, 1.0));

  inMesh->createEdge(inV1, inV2);
  inMesh->createEdge(inV2, inV3);
  inMesh->createEdge(inV3, inV4);
  inMesh->createEdge(inV1, inV4);

  inMesh->allocateDataValues();
  addGlobalIndex(inMesh);
  inMesh->setGlobalNumberOfVertices(inMesh->vertices().size());

  auto &inValues = inData->values();
  inValues << 1.0, 2.0, 2.0, 1.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", 1, 1_dataID);
  int           outDataID = outData->getID();
  auto &        outV1     = outMesh->createVertex(Vector2d(0.0, 0.0));
  auto &        outV2     = outMesh->createVertex(Vector2d(0.0, 1.0));
  auto &        outV3     = outMesh->createVertex(Vector2d(1.1, 1.1));
  auto &        outV4     = outMesh->createVertex(Vector2d(0.1, 1.1));
  outMesh->createEdge(outV1, outV2);
  outMesh->createEdge(outV2, outV3);
  outMesh->createEdge(outV3, outV4);
  outMesh->createEdge(outV1, outV4);
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);
  outMesh->setGlobalNumberOfVertices(outMesh->vertices().size());

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  mapping.computeMapping();
  mapping.map(inDataID, outDataID);

  testSerialScaledConsistent(inMesh, outMesh, inData, outData);
}

void perform3DTestScaledConsistentMapping(Mapping &mapping)
{
  int dimensions = 3;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", 1, 0_dataID);
  int           inDataID = inData->getID();
  auto &        inV1     = inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  auto &        inV2     = inMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  auto &        inV3     = inMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.5));
  auto &        inV4     = inMesh->createVertex(Eigen::Vector3d(2.0, 0.0, 0.0));
  auto &        inV5     = inMesh->createVertex(Eigen::Vector3d(0.0, 2.0, 0.0));
  auto &        inV6     = inMesh->createVertex(Eigen::Vector3d(0.0, 2.0, 1.0));
  auto &        inE1     = inMesh->createEdge(inV1, inV2);
  auto &        inE2     = inMesh->createEdge(inV2, inV3);
  auto &        inE3     = inMesh->createEdge(inV1, inV3);
  auto &        inE4     = inMesh->createEdge(inV4, inV5);
  auto &        inE5     = inMesh->createEdge(inV5, inV6);
  auto &        inE6     = inMesh->createEdge(inV4, inV6);
  inMesh->createTriangle(inE1, inE2, inE3);
  inMesh->createTriangle(inE4, inE5, inE6);

  inMesh->allocateDataValues();
  addGlobalIndex(inMesh);
  inMesh->setGlobalNumberOfVertices(inMesh->vertices().size());

  auto &inValues = inData->values();
  inValues << 1.0, 2.0, 4.0, 6.0, 8.0, 9.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", 1, 1_dataID);
  int           outDataID = outData->getID();
  auto &        outV1     = outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  auto &        outV2     = outMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  auto &        outV3     = outMesh->createVertex(Eigen::Vector3d(0.0, 1.1, 0.6));
  auto &        outE1     = outMesh->createEdge(outV1, outV2);
  auto &        outE2     = outMesh->createEdge(outV2, outV3);
  auto &        outE3     = outMesh->createEdge(outV1, outV3);
  outMesh->createTriangle(outE1, outE2, outE3);

  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);
  outMesh->setGlobalNumberOfVertices(outMesh->vertices().size());

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);
  mapping.computeMapping();
  BOOST_TEST(mapping.hasComputedMapping() == true);
  mapping.map(inDataID, outDataID);

  testSerialScaledConsistent(inMesh, outMesh, inData, outData);
}

void perform2DTestConservativeMapping(Mapping &mapping)
{
  const int    dimensions = 2;
  const double tolerance  = 1e-6;
  using Eigen::Vector2d;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", 1, 0_dataID);
  int           inDataID = inData->getID();
  mesh::Vertex &vertex0  = inMesh->createVertex(Vector2d(0, 0));
  mesh::Vertex &vertex1  = inMesh->createVertex(Vector2d(0, 0));
  inMesh->allocateDataValues();
  inData->values() << 1.0, 2.0;
  addGlobalIndex(inMesh);
  inMesh->setGlobalNumberOfVertices(inMesh->vertices().size());

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", 1, 1_dataID);
  int           outDataID = outData->getID();
  outMesh->createVertex(Vector2d(0.0, 0.0));
  outMesh->createVertex(Vector2d(1.0, 0.0));
  outMesh->createVertex(Vector2d(1.0, 1.0));
  outMesh->createVertex(Vector2d(0.0, 1.0));
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);
  outMesh->setGlobalNumberOfVertices(outMesh->vertices().size());

  auto &values = outData->values();

  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  vertex0.setCoords(Vector2d(0.5, 0.0));
  vertex1.setCoords(Vector2d(0.5, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(testing::equals(values, Eigen::Vector4d(0.5, 0.5, 1.0, 1.0), tolerance));

  vertex0.setCoords(Vector2d(0.0, 0.5));
  vertex1.setCoords(Vector2d(1.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(testing::equals(values, Eigen::Vector4d(0.5, 1.0, 1.0, 0.5), tolerance));

  vertex0.setCoords(Vector2d(0.0, 1.0));
  vertex1.setCoords(Vector2d(1.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(testing::equals(values, Eigen::Vector4d(0.0, 2.0, 0.0, 1.0), tolerance));

  vertex0.setCoords(Vector2d(0.0, 0.0));
  vertex1.setCoords(Vector2d(1.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(testing::equals(values, Eigen::Vector4d(1.0, 0.0, 2.0, 0.0), tolerance));

  vertex0.setCoords(Vector2d(0.4, 0.5));
  vertex1.setCoords(Vector2d(0.6, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(values.sum() == 3.0);
}

void perform2DTestConservativeMappingVector(Mapping &mapping)
{
  const int    dimensions = 2;
  const double tolerance  = 1e-6;
  using Eigen::Vector2d;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", 2, 0_dataID);
  int           inDataID = inData->getID();
  mesh::Vertex &vertex0  = inMesh->createVertex(Vector2d(0, 0));
  mesh::Vertex &vertex1  = inMesh->createVertex(Vector2d(0, 0));
  inMesh->allocateDataValues();
  inData->values() << 1.0, 4.0, 2.0, 5.0;
  addGlobalIndex(inMesh);
  inMesh->setGlobalNumberOfVertices(inMesh->vertices().size());

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", 2, 1_dataID);
  int           outDataID = outData->getID();
  outMesh->createVertex(Vector2d(0.0, 0.0));
  outMesh->createVertex(Vector2d(1.0, 0.0));
  outMesh->createVertex(Vector2d(1.0, 1.0));
  outMesh->createVertex(Vector2d(0.0, 1.0));
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);
  outMesh->setGlobalNumberOfVertices(outMesh->vertices().size());

  auto &values = outData->values();

  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  vertex0.setCoords(Vector2d(0.5, 0.0));
  vertex1.setCoords(Vector2d(0.5, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  Eigen::VectorXd refValues(8);
  refValues << 0.5, 2, 0.5, 2, 1.0, 2.5, 1.0, 2.5;
  BOOST_TEST(testing::equals(values, refValues, tolerance)); // fails

  vertex0.setCoords(Vector2d(0.0, 0.5));
  vertex1.setCoords(Vector2d(1.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  refValues << 0.5, 2, 1.0, 2.5, 1.0, 2.5, 0.5, 2;
  BOOST_TEST(testing::equals(values, refValues, tolerance)); // fails

  vertex0.setCoords(Vector2d(0.0, 1.0));
  vertex1.setCoords(Vector2d(1.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  refValues << 0.0, 0.0, 2.0, 5.0, 0.0, 0.0, 1.0, 4.0;
  BOOST_TEST(testing::equals(values, refValues, tolerance)); // fails

  vertex0.setCoords(Vector2d(0.0, 0.0));
  vertex1.setCoords(Vector2d(1.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  refValues << 1.0, 4.0, 0.0, 0.0, 2.0, 5.0, 0.0, 0.0;
  BOOST_TEST(testing::equals(values, refValues, tolerance)); // fails

  vertex0.setCoords(Vector2d(0.4, 0.5));
  vertex1.setCoords(Vector2d(0.6, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(values.sum() == 12.0); // fails
}

void perform3DTestConservativeMapping(Mapping &mapping)
{
  using Eigen::Vector3d;
  int dimensions = 3;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", 1, 0_dataID);
  int           inDataID = inData->getID();
  mesh::Vertex &vertex0  = inMesh->createVertex(Vector3d(0, 0, 0));
  mesh::Vertex &vertex1  = inMesh->createVertex(Vector3d(0, 0, 0));
  inMesh->allocateDataValues();
  inData->values() << 1.0, 2.0;
  addGlobalIndex(inMesh);
  inMesh->setGlobalNumberOfVertices(inMesh->vertices().size());

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", 1, 1_dataID);
  int           outDataID = outData->getID();
  outMesh->createVertex(Vector3d(0.0, 0.0, 0.0));
  outMesh->createVertex(Vector3d(1.0, 0.0, 0.0));
  outMesh->createVertex(Vector3d(1.0, 1.0, 0.0));
  outMesh->createVertex(Vector3d(0.0, 1.0, 0.0));
  outMesh->createVertex(Vector3d(0.0, 0.0, 1.0));
  outMesh->createVertex(Vector3d(1.0, 0.0, 1.0));
  outMesh->createVertex(Vector3d(1.0, 1.0, 1.0));
  outMesh->createVertex(Vector3d(0.0, 1.0, 1.0));
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);
  outMesh->setGlobalNumberOfVertices(outMesh->vertices().size());

  auto & values      = outData->values();
  double expectedSum = inData->values().sum();

  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  vertex0.setCoords(Vector3d(0.5, 0.0, 0.0));
  vertex1.setCoords(Vector3d(0.5, 1.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping());
  BOOST_TEST(values.sum() == expectedSum);
}
} // namespace

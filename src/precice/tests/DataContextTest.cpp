#include <Eigen/Core>
#include <string>
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "precice/impl/ReadDataContext.hpp"
#include "precice/impl/WriteDataContext.hpp"
#include "testing/DataContextFixture.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::impl;

BOOST_AUTO_TEST_SUITE(PreciceTests)

BOOST_AUTO_TEST_SUITE(DataContextTests)

BOOST_AUTO_TEST_CASE(testDataContextWriteMapping)
{
  PRECICE_TEST(1_rank);

  testing::DataContextFixture fixture;

  // Create mesh object for from mesh
  int           dimensions  = 3;
  mesh::PtrMesh ptrFromMesh = std::make_shared<mesh::Mesh>("ParticipantMesh", dimensions, testing::nextMeshID());
  mesh::PtrData ptrFromData = ptrFromMesh->createData("MappedData", dimensions, 0_dataID);

  ptrFromMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  ptrFromMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  ptrFromMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 1.0));

  // Create mesh object for from mesh
  mesh::PtrMesh ptrToMesh = std::make_shared<mesh::Mesh>("OtherMesh", dimensions, testing::nextMeshID());
  mesh::PtrData ptrToData = ptrToMesh->createData("MappedData", dimensions, 1_dataID);

  ptrToMesh->createVertex(Eigen::Vector3d(0.0, 0.1, 0.0));
  ptrToMesh->createVertex(Eigen::Vector3d(1.0, 0.1, 0.0));
  ptrToMesh->createVertex(Eigen::Vector3d(0.0, 0.1, 1.0));

  MeshContext toMeshContext;
  toMeshContext.mesh = ptrToMesh;

  WriteDataContext dataContext(ptrFromData, ptrFromMesh);

  MappingContext mappingContext;
  mappingContext.fromMeshID = ptrFromMesh->getID();
  mappingContext.toMeshID   = ptrToMesh->getID();

  BOOST_TEST(ptrToData->getID() != ptrFromData->getID());
  BOOST_TEST(ptrToMesh->getID() != ptrFromMesh->getID());

  BOOST_TEST(!fixture.hasMapping(dataContext));
  BOOST_TEST(fixture.getProvidedDataID(dataContext) == ptrFromData->getID());
  BOOST_TEST(dataContext.getMeshID() == ptrFromMesh->getID());

  dataContext.appendMappingConfiguration(mappingContext, toMeshContext);

  // mapping is configured. Write mapping, therefore _providedData == _fromData
  BOOST_TEST(fixture.hasMapping(dataContext));
  BOOST_TEST(fixture.getFromDataID(dataContext, 0) == ptrFromData->getID());
  BOOST_TEST(fixture.getToDataID(dataContext, 0) == ptrToData->getID());
  BOOST_TEST(fixture.getProvidedDataID(dataContext) != ptrToData->getID());
  BOOST_TEST(fixture.getProvidedDataID(dataContext) == ptrFromData->getID());
  BOOST_TEST(dataContext.getMeshID() != ptrToMesh->getID());
  BOOST_TEST(dataContext.getMeshID() == ptrFromMesh->getID());
  BOOST_TEST(fixture.hasWriteMapping(dataContext));
  BOOST_TEST(!fixture.hasReadMapping(dataContext));
  BOOST_TEST(fixture.mappingContexts(dataContext)[0].fromMeshID == mappingContext.fromMeshID);
  BOOST_TEST(fixture.mappingContexts(dataContext)[0].toMeshID == mappingContext.toMeshID);
  BOOST_TEST(fixture.mappingContexts(dataContext)[0].mapping == mappingContext.mapping);
}

BOOST_AUTO_TEST_CASE(testDataContextMultipleWriteMapping)
{
  PRECICE_TEST(1_rank);

  testing::DataContextFixture fixture;

  // Create mesh object for from mesh
  int           dimensions  = 3;
  mesh::PtrMesh ptrFromMesh = std::make_shared<mesh::Mesh>("ParticipantMesh", dimensions, testing::nextMeshID());
  mesh::PtrData ptrFromData = ptrFromMesh->createData("MappedData", dimensions, 0_dataID);

  ptrFromMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  ptrFromMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  ptrFromMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 1.0));

  // Create mesh object for to mesh
  mesh::PtrMesh ptrToMesh = std::make_shared<mesh::Mesh>("OtherMesh", dimensions, testing::nextMeshID());
  mesh::PtrData ptrToData = ptrToMesh->createData("MappedData", dimensions, 1_dataID);

  ptrToMesh->createVertex(Eigen::Vector3d(0.0, 0.1, 0.0));
  ptrToMesh->createVertex(Eigen::Vector3d(1.0, 0.1, 0.0));
  ptrToMesh->createVertex(Eigen::Vector3d(0.0, 0.1, 1.0));

  MeshContext toMeshContext1;
  toMeshContext1.mesh = ptrToMesh;

  WriteDataContext dataContext(ptrFromData, ptrFromMesh);

  MappingContext mappingContext;
  mappingContext.fromMeshID = ptrFromMesh->getID();
  mappingContext.toMeshID   = ptrToMesh->getID();

  BOOST_TEST(ptrToData->getID() != ptrFromData->getID());
  BOOST_TEST(ptrToMesh->getID() != ptrFromMesh->getID());

  BOOST_TEST(!fixture.hasMapping(dataContext));
  BOOST_TEST(fixture.getProvidedDataID(dataContext) == ptrFromData->getID());
  BOOST_TEST(dataContext.getMeshID() == ptrFromMesh->getID());

  // Add the first mapping we configured for this context
  dataContext.appendMappingConfiguration(mappingContext, toMeshContext1);

  // Add a second mapping targeting a different to mesh
  // Create the object for to mesh and the data
  mesh::PtrMesh ptrToMesh2 = std::make_shared<mesh::Mesh>("SecondOtherMesh", dimensions, testing::nextMeshID());
  mesh::PtrData ptrToData2 = ptrToMesh2->createData("MappedData", dimensions, 2_dataID);

  ptrToMesh2->createVertex(Eigen::Vector3d(0.0, 1.1, 0.0));
  ptrToMesh2->createVertex(Eigen::Vector3d(2.0, 1.1, 0.0));
  ptrToMesh2->createVertex(Eigen::Vector3d(0.0, 2.1, 4.0));

  MeshContext toMeshContext2;
  toMeshContext2.mesh = ptrToMesh2;

  MappingContext mappingContext2;
  mappingContext2.fromMeshID = ptrFromMesh->getID();
  mappingContext2.toMeshID   = ptrToMesh2->getID();

  // the mapping configuration
  dataContext.appendMappingConfiguration(mappingContext2, toMeshContext2);

  // First, we repeat the checks from above in order to check that nothing changed
  BOOST_TEST(fixture.hasMapping(dataContext));
  BOOST_TEST(fixture.getProvidedDataID(dataContext) != ptrToData->getID());
  BOOST_TEST(fixture.getProvidedDataID(dataContext) == ptrFromData->getID());
  BOOST_TEST(dataContext.getMeshID() != ptrToMesh->getID());
  BOOST_TEST(dataContext.getMeshID() == ptrFromMesh->getID());
  BOOST_TEST(fixture.hasWriteMapping(dataContext));
  BOOST_TEST(!fixture.hasReadMapping(dataContext));
  // Test dedicated content of the first mapping configuration
  BOOST_TEST(fixture.getFromDataID(dataContext, 0) == ptrFromData->getID());
  BOOST_TEST(fixture.getToDataID(dataContext, 0) == ptrToData->getID());
  BOOST_TEST(fixture.mappingContexts(dataContext)[0].fromMeshID == mappingContext.fromMeshID);
  BOOST_TEST(fixture.mappingContexts(dataContext)[0].toMeshID == mappingContext.toMeshID);
  BOOST_TEST(fixture.mappingContexts(dataContext)[0].mapping == mappingContext.mapping);

  // Now, test the newly added mapping
  BOOST_TEST(fixture.getProvidedDataID(dataContext) != ptrToData2->getID());
  BOOST_TEST(dataContext.getMeshID() != ptrToMesh2->getID());
  BOOST_TEST(fixture.hasWriteMapping(dataContext));
  BOOST_TEST(!fixture.hasReadMapping(dataContext));
  // Test dedicated content of the first mapping configuration
  BOOST_TEST(fixture.getFromDataID(dataContext, 1) == ptrFromData->getID());
  BOOST_TEST(fixture.getToDataID(dataContext, 1) == ptrToData2->getID());
  BOOST_TEST(fixture.mappingContexts(dataContext)[1].fromMeshID == mappingContext2.fromMeshID);
  BOOST_TEST(fixture.mappingContexts(dataContext)[1].toMeshID == mappingContext2.toMeshID);
  BOOST_TEST(fixture.mappingContexts(dataContext)[1].mapping == mappingContext2.mapping);
}

BOOST_AUTO_TEST_CASE(testDataContextReadMapping)
{
  PRECICE_TEST(1_rank);

  testing::DataContextFixture fixture;

  // Create mesh object
  int           dimensions = 3;
  mesh::PtrMesh ptrToMesh  = std::make_shared<mesh::Mesh>("ParticipantMesh", dimensions, testing::nextMeshID());
  mesh::PtrData ptrToData  = ptrToMesh->createData("MappedData", dimensions, 0_dataID);

  ptrToMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  ptrToMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  ptrToMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 1.0));

  // Create mesh object for from mesh
  mesh::PtrMesh ptrFromMesh = std::make_shared<mesh::Mesh>("OtherMesh", dimensions, testing::nextMeshID());
  mesh::PtrData ptrFromData = ptrFromMesh->createData("MappedData", dimensions, 1_dataID);

  ptrFromMesh->createVertex(Eigen::Vector3d(0.0, 0.1, 0.0));
  ptrFromMesh->createVertex(Eigen::Vector3d(1.0, 0.1, 0.0));
  ptrFromMesh->createVertex(Eigen::Vector3d(0.0, 0.1, 1.0));

  MeshContext fromMeshContext;
  fromMeshContext.mesh = ptrFromMesh;

  ReadDataContext dataContext(ptrToData, ptrToMesh);

  MappingContext mappingContext;
  mappingContext.fromMeshID = ptrFromMesh->getID();
  mappingContext.toMeshID   = ptrToMesh->getID();

  BOOST_TEST(ptrToData->getID() != ptrFromData->getID());
  BOOST_TEST(ptrToMesh->getID() != ptrFromMesh->getID());

  BOOST_TEST(!fixture.hasMapping(dataContext));
  BOOST_TEST(fixture.getProvidedDataID(dataContext) == ptrToData->getID());
  BOOST_TEST(dataContext.getMeshID() == ptrToMesh->getID());

  // Add the mapping we configured for this context
  // For read data contexts, there is only one context allowed
  dataContext.appendMappingConfiguration(mappingContext, fromMeshContext);

  // mapping is configured. Write mapping, therefore _providedData == _toData
  BOOST_TEST(fixture.hasMapping(dataContext));
  BOOST_TEST(fixture.getFromDataID(dataContext, 0) == ptrFromData->getID());
  BOOST_TEST(fixture.getToDataID(dataContext, 0) == ptrToData->getID());
  BOOST_TEST(fixture.getProvidedDataID(dataContext) == ptrToData->getID());
  BOOST_TEST(fixture.getProvidedDataID(dataContext) != ptrFromData->getID());
  BOOST_TEST(dataContext.getMeshID() == ptrToMesh->getID());
  BOOST_TEST(dataContext.getMeshID() != ptrFromMesh->getID());
  BOOST_TEST(!fixture.hasWriteMapping(dataContext));
  BOOST_TEST(fixture.hasReadMapping(dataContext));
  BOOST_TEST(fixture.mappingContexts(dataContext)[0].fromMeshID == mappingContext.fromMeshID);
  BOOST_TEST(fixture.mappingContexts(dataContext)[0].toMeshID == mappingContext.toMeshID);
  BOOST_TEST(fixture.mappingContexts(dataContext)[0].mapping == mappingContext.mapping);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

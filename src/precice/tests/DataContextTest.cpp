#include <Eigen/Core>
#include <string>
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "precice/impl/ReadDataContext.hpp"
#include "precice/impl/WriteDataContext.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::impl;

BOOST_AUTO_TEST_SUITE(PreciceTests)

BOOST_AUTO_TEST_SUITE(DataContextTests)

BOOST_AUTO_TEST_CASE(testDataContextWriteMapping)
{
  PRECICE_TEST(1_rank);

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

  MeshContext toMeshContext(dimensions);
  toMeshContext.mesh = ptrToMesh;

  WriteDataContext dataContext(ptrFromData, ptrFromMesh);

  MappingContext mappingContext;
  mappingContext.fromMeshID = ptrFromMesh->getID();
  mappingContext.toMeshID   = ptrToMesh->getID();

  BOOST_TEST(ptrToData->getID() != ptrFromData->getID());
  BOOST_TEST(ptrToMesh->getID() != ptrFromMesh->getID());

  BOOST_TEST(!dataContext.hasMapping());
  BOOST_TEST(dataContext.getProvidedDataID() == ptrFromData->getID());
  BOOST_TEST(dataContext.getMeshID() == ptrFromMesh->getID());

  dataContext.configureMapping(mappingContext, toMeshContext);

  // mapping is configured. Write mapping, therefore _providedData == _fromData
  BOOST_TEST(dataContext.hasMapping());
  BOOST_TEST(dataContext.getFromDataID() == ptrFromData->getID());
  BOOST_TEST(dataContext.getToDataID() == ptrToData->getID());
  BOOST_TEST(dataContext.getProvidedDataID() != ptrToData->getID());
  BOOST_TEST(dataContext.getProvidedDataID() == ptrFromData->getID());
  BOOST_TEST(dataContext.getMeshID() != ptrToMesh->getID());
  BOOST_TEST(dataContext.getMeshID() == ptrFromMesh->getID());
  BOOST_TEST(dataContext.hasWriteMapping());
  BOOST_TEST(!dataContext.hasReadMapping());
  BOOST_TEST(dataContext.mappingContext().fromMeshID == mappingContext.fromMeshID);
  BOOST_TEST(dataContext.mappingContext().toMeshID == mappingContext.toMeshID);
  BOOST_TEST(dataContext.mappingContext().hasMappedData == mappingContext.hasMappedData);
  BOOST_TEST(dataContext.mappingContext().mapping == mappingContext.mapping);
  BOOST_TEST(dataContext.mappingContext().timing == mappingContext.timing);
}

BOOST_AUTO_TEST_CASE(testDataContextReadMapping)
{
  PRECICE_TEST(1_rank);

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

  MeshContext fromMeshContext(dimensions);
  fromMeshContext.mesh = ptrFromMesh;

  ReadDataContext dataContext(ptrToData, ptrToMesh);

  MappingContext mappingContext;
  mappingContext.fromMeshID = ptrFromMesh->getID();
  mappingContext.toMeshID   = ptrToMesh->getID();

  BOOST_TEST(ptrToData->getID() != ptrFromData->getID());
  BOOST_TEST(ptrToMesh->getID() != ptrFromMesh->getID());

  BOOST_TEST(!dataContext.hasMapping());
  BOOST_TEST(dataContext.getProvidedDataID() == ptrToData->getID());
  BOOST_TEST(dataContext.getMeshID() == ptrToMesh->getID());

  dataContext.configureMapping(mappingContext, fromMeshContext);

  // mapping is configured. Write mapping, therefore _providedData == _toData
  BOOST_TEST(dataContext.hasMapping());
  BOOST_TEST(dataContext.getFromDataID() == ptrFromData->getID());
  BOOST_TEST(dataContext.getToDataID() == ptrToData->getID());
  BOOST_TEST(dataContext.getProvidedDataID() == ptrToData->getID());
  BOOST_TEST(dataContext.getProvidedDataID() != ptrFromData->getID());
  BOOST_TEST(dataContext.getMeshID() == ptrToMesh->getID());
  BOOST_TEST(dataContext.getMeshID() != ptrFromMesh->getID());
  BOOST_TEST(!dataContext.hasWriteMapping());
  BOOST_TEST(dataContext.hasReadMapping());
  BOOST_TEST(dataContext.mappingContext().fromMeshID == mappingContext.fromMeshID);
  BOOST_TEST(dataContext.mappingContext().toMeshID == mappingContext.toMeshID);
  BOOST_TEST(dataContext.mappingContext().hasMappedData == mappingContext.hasMappedData);
  BOOST_TEST(dataContext.mappingContext().mapping == mappingContext.mapping);
  BOOST_TEST(dataContext.mappingContext().timing == mappingContext.timing);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

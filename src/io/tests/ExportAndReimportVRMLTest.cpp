// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ExportAndReimportVRMLTest.hpp"
#include "io/ExportVRML.hpp"
#include "io/ImportVRML.hpp"
#include "io/ExportVTK.hpp"
#include "mesh/Data.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/PropertyContainer.hpp"
#include "geometry/Cuboid.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"
#include "utils/Dimensions.hpp"
#include "tarch/la/WrappedVector.h"
#include <map>
#include <string>

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::io::tests::ExportAndReimportVRMLTest)

namespace precice {
namespace io {
namespace tests {


tarch::logging::Log ExportAndReimportVRMLTest::
   _log ( "precice::io::tests::ExportAndReimportVRMLTest" );

ExportAndReimportVRMLTest:: ExportAndReimportVRMLTest()
:
  TestCase ( "io::tests::ExportAndReimportVRMLTest" )
{}

void ExportAndReimportVRMLTest:: run()
{
# ifndef PRECICE_NO_SPIRIT2
  PRECICE_MASTER_ONLY {
    testMethod ( testInternallyCreatedGeometry );
    //testMethod ( testReimportDriftRatchet );
  }
# endif // not PRECICE_NO_SPIRIT2
}

void ExportAndReimportVRMLTest:: testInternallyCreatedGeometry()
{
  preciceTrace("testInternallyCreatedGeometry()");
  mesh::Mesh::resetGeometryIDsGlobally();
  // Create geometry
  bool flipNormals = false;
  int dim = 2;
  mesh::Mesh mesh("TestCuboid", dim, flipNormals);
  utils::DynVector offset(2, 0.0);
  double dx = 1.0;
  utils::DynVector length(2, 1.0);
  geometry::Cuboid cuboid(offset, dx, length);
  // Set geoemtry sub-ids
  std::string nameSubID0("side-0");
  mesh.setSubID(nameSubID0);

  std::string nameSubID1 = "side-1";
  mesh.setSubID(nameSubID1);

  std::string nameSubID2 = "side-2";
  mesh.setSubID(nameSubID2);
  mesh::PtrData data = mesh.createData("TestData", 2);
  cuboid.create(mesh);

  // Query mesh information
  size_t vertexCount = mesh.vertices().size();
  size_t edgeCount = mesh.edges().size();
  utils::Vector2D coordsVertexOne(mesh.vertices()[0].getCoords());
  utils::Vector2D coordsVertexN(mesh.vertices()[vertexCount-1].getCoords());
  utils::Vector2D dataOne ( 2.0 );
  utils::Vector2D dataN ( 4.0 );
  int vertex0ID = mesh.vertices()[0].getID();
  int vertexNID = mesh.vertices()[vertexCount-1].getID();
  using tarch::la::slice;
  slice<2>(data->values(), vertex0ID*2) = dataOne;
  slice<2>(data->values(), vertexNID*2) = dataN;
  std::string filename("io-ExportandReimportVRMLTest-testInternallyCreatedGeometry.wrl");

  // Export geometry
  bool exportNormals = false;
  io::ExportVRML exportMesh(exportNormals);
  exportMesh.doExportCheckpoint(filename, mesh);

  // Reimport geometry
  mesh::Mesh reimportedMesh("TestCuboid", 2, false);
  mesh::PtrData reimportedData =
      reimportedMesh.createData(data->getName(), data->getDimensions());
  reimportedMesh.setSubID(nameSubID0);
  reimportedMesh.setSubID(nameSubID1);
  reimportedMesh.setSubID(nameSubID2);

  std::string location = "";
  io::ImportVRML importMesh(location);
  importMesh.doImportCheckpoint(filename, reimportedMesh, true);

  // Validate mesh information
  validateEquals(reimportedMesh.vertices().size(), vertexCount);
  validateEquals(reimportedMesh.edges().size(), edgeCount);
  validate(tarch::la::equals(reimportedMesh.vertices()[0].getCoords(),
                             coordsVertexOne));
  validate(tarch::la::equals(
           reimportedMesh.vertices()[vertexCount-1].getCoords(), coordsVertexN));
  utils::DynVector& values = reimportedData->values();
  int vertexID = reimportedMesh.vertices()[0].getID();
  validateNumericalEquals(values[vertexID], dataOne[0]);
  validateNumericalEquals(values[vertexID], dataOne[1]);
  utils::Vector2D readDataOne;
  utils::Vector2D readDataN;
  for (size_t i=0; i < 2; i++){
    readDataOne[i] = data->values()[vertex0ID * 2 + i];
    readDataN[i] = data->values()[vertexNID * 2 + i];
  }
  validate(tarch::la::equals(readDataOne, dataOne));
  validate(tarch::la::equals(readDataN, dataN));

  validateEquals(reimportedMesh.propertyContainers().size(), 3);
  std::vector<int> properties;
  int id = mesh::PropertyContainer::INDEX_GEOMETRY_ID;
  reimportedMesh.vertices()[0].getProperties(id, properties);
  validateEqualsWithMessage(properties.size(), 3, properties);
  properties.clear();

  reimportedMesh.vertices()[1].getProperties(id, properties);
  validateEqualsWithMessage(properties.size(), 3, properties);
  properties.clear();

  reimportedMesh.vertices()[2].getProperties(id, properties);
  validateEqualsWithMessage(properties.size(), 2, properties);
  properties.clear();

  reimportedMesh.vertices()[3].getProperties(id, properties);
  validateEqualsWithMessage(properties.size(), 2, properties);
  properties.clear();

  reimportedMesh.edges()[0].getProperties(id, properties);
  validateEqualsWithMessage(properties.size(), 2, properties);
  properties.clear();

  reimportedMesh.edges()[1].getProperties(id, properties);
  validateEqualsWithMessage(properties.size(), 1, properties);
  properties.clear();

  reimportedMesh.edges()[2].getProperties(id, properties);
  validateEqualsWithMessage(properties.size(), 2, properties);
  properties.clear();

  reimportedMesh.edges()[3].getProperties(id, properties);
  validateEqualsWithMessage(properties.size(), 2, properties);
  properties.clear();

  mesh::Mesh reimportedMesh1("TestCuboid", 2, false);
  importMesh.doImport(filename, reimportedMesh1);
  validateEquals(reimportedMesh1.data().size(), 0);
  validateEquals(reimportedMesh1.propertyContainers().size(), 0);
}

void ExportAndReimportVRMLTest:: testReimportDriftRatchet()
{
  preciceTrace("testReimportDriftRatchet()");
  mesh::Mesh mesh("test-cuboid", 3, false);
  io::ImportVRML importMesh("");
  importMesh.doImport("io-ExportVRMLTest-testExportCuboid-3d.wrl", mesh);
  mesh.computeState();
  ExportVTK exportVTK(true);
  exportVTK.doExport("io-ExportAndReimportVRMLTest-testReimportDriftRatchet.vtk", mesh);
}


}}} // namespace precice, io, tests

// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ScaleActionTest.hpp"
#include "action/ScaleByAreaAction.hpp"
#include "action/ScaleByDtAction.hpp"
#include "action/config/ActionConfiguration.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "mesh/Data.hpp"
#include "utils/Parallel.hpp"
#include "utils/Dimensions.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::action::tests::ScaleActionTest)

namespace precice {
namespace action {
namespace tests {

tarch::logging::Log ScaleActionTest::
  _log ("precice::action::tests::ScaleActionTest");

ScaleActionTest:: ScaleActionTest()
:
   TestCase ("precice::action::tests::ScaleActionTest")
{}

void ScaleActionTest:: run()
{
  PRECICE_MASTER_ONLY {
    testMethod(testDivideByArea);
    testMethod(testScaleByComputedTimestepLength);
    testMethod(testScaleByComputedTimestepPartLength);
    testMethod(testConfiguration);
  }
}

void ScaleActionTest:: testDivideByArea()
{
   preciceTrace("testDivideByArea()");
   using namespace mesh;
   PtrMesh mesh(new Mesh("Mesh", 2, true));
   PtrData data = mesh->createData("test-data", 1);
   int dataID = data->getID();
   Vertex& v0 = mesh->createVertex(utils::Vector2D(0.0, 0.0));
   Vertex& v1 = mesh->createVertex(utils::Vector2D(1.0, 0.0));
   Vertex& v2 = mesh->createVertex(utils::Vector2D(1.0, 1.0));
   mesh->createEdge(v0, v1);
   mesh->createEdge(v1, v2);
   mesh->computeState();
   mesh->allocateDataValues();
   utils::DynVector& values = data->values();
   assignList(values) = 2.0, 3.0, 4.0;

   validateNumericalEquals(values[0], 2.0);
   validateNumericalEquals(values[1], 3.0);
   validateNumericalEquals(values[2], 4.0);

   // Scale properties
   action::ScaleByAreaAction scale(
     action::ScaleByAreaAction::ALWAYS_PRIOR, dataID, mesh,
     action::ScaleByAreaAction::SCALING_DIVIDE_BY_AREA);
//   scale.loadMeshContext ( geoContext );
   scale.performAction(0.0, 0.0, 0.0, 0.0);

   // Validate expected results
   validateNumericalEquals(values[0], 4.0);
   validateNumericalEquals(values[1], 3.0);
   validateNumericalEquals(values[2], 8.0);
}

void ScaleActionTest:: testScaleByComputedTimestepLength()
{
  preciceTrace("testScaleByComputedTimestepLength()");
  using namespace mesh;
  PtrMesh mesh(new Mesh("Mesh", 3, true));
  PtrData sourceData = mesh->createData("SourceData", 1);
  PtrData targetData = mesh->createData("TargetData", 1);
  int sourceDataID = sourceData->getID();
  int targetDataID = targetData->getID();
  mesh->createVertex(utils::Vector3D(0.0));
  mesh->createVertex(utils::Vector3D(1.0));
  mesh->createVertex(utils::Vector3D(2.0));
//  mesh->computeState ();
  mesh->allocateDataValues();
  utils::DynVector& sourceValues = sourceData->values();
  utils::DynVector& targetValues = targetData->values();
  assignList(sourceValues) = 2.0, 3.0, 4.0;
  assign(targetValues) = 0.0;

  action::ScaleByDtAction scale(
      action::ScaleByDtAction::ALWAYS_PRIOR, sourceDataID, targetDataID, mesh,
      action::ScaleByDtAction::SCALING_BY_COMPUTED_DT_RATIO);

  scale.performAction(0.0, 0.0, 0.0, 1.0);
  validateNumericalEquals(sourceValues[0], 2.0);
  validateNumericalEquals(sourceValues[1], 3.0);
  validateNumericalEquals(sourceValues[2], 4.0);
  validateNumericalEquals(targetValues[0], 0.0);
  validateNumericalEquals(targetValues[1], 0.0);
  validateNumericalEquals(targetValues[2], 0.0);

  scale.performAction(0.0, 0.5, 0.5, 1.0);
  validateNumericalEquals(sourceValues[0], 2.0);
  validateNumericalEquals(sourceValues[1], 3.0);
  validateNumericalEquals(sourceValues[2], 4.0);
  validateNumericalEquals(targetValues[0], 1.0);
  validateNumericalEquals(targetValues[1], 1.5);
  validateNumericalEquals(targetValues[2], 2.0);

  scale.performAction(0.0, 0.25, 0.75, 1.0);
  validateNumericalEquals(sourceValues[0], 2.0);
  validateNumericalEquals(sourceValues[1], 3.0);
  validateNumericalEquals(sourceValues[2], 4.0);
  validateNumericalEquals(targetValues[0], 0.5);
  validateNumericalEquals(targetValues[1], 0.75);
  validateNumericalEquals(targetValues[2], 1.0);

  scale.performAction(0.0, 0.25, 1.0, 1.0);
  validateNumericalEquals(sourceValues[0], 2.0);
  validateNumericalEquals(sourceValues[1], 3.0);
  validateNumericalEquals(sourceValues[2], 4.0);
  validateNumericalEquals(targetValues[0], 0.5);
  validateNumericalEquals(targetValues[1], 0.75);
  validateNumericalEquals(targetValues[2], 1.0);
}

void ScaleActionTest:: testScaleByComputedTimestepPartLength()
{
  preciceTrace("testScaleByComputedTimestepPartLength()");
  using namespace mesh;
  PtrMesh mesh(new Mesh("Mesh", 3, true));
  PtrData sourceData = mesh->createData("SourceData", 1);
  PtrData targetData = mesh->createData("TargetData", 1);
  int sourceDataID = sourceData->getID();
  int targetDataID = targetData->getID();
  mesh->createVertex(utils::Vector3D(0.0));
  mesh->createVertex(utils::Vector3D(1.0));
  mesh->createVertex(utils::Vector3D(2.0));
  mesh->allocateDataValues();
  utils::DynVector& sourceValues = sourceData->values();
  utils::DynVector& targetValues = targetData->values();
  assignList(sourceValues) = 2.0, 3.0, 4.0;
  assign(targetValues) = 0.0;

  action::ScaleByDtAction scale(
      action::ScaleByDtAction::ALWAYS_PRIOR, sourceDataID, targetDataID, mesh,
      action::ScaleByDtAction::SCALING_BY_COMPUTED_DT_PART_RATIO);

  scale.performAction(0.0, 0.0, 0.0, 1.0);
  validateNumericalEquals(sourceValues[0], 2.0);
  validateNumericalEquals(sourceValues[1], 3.0);
  validateNumericalEquals(sourceValues[2], 4.0);
  validateNumericalEquals(targetValues[0], 0.0);
  validateNumericalEquals(targetValues[1], 0.0);
  validateNumericalEquals(targetValues[2], 0.0);

  scale.performAction(0.0, 0.5, 0.5, 1.0);
  validateNumericalEquals(sourceValues[0], 2.0);
  validateNumericalEquals(sourceValues[1], 3.0);
  validateNumericalEquals(sourceValues[2], 4.0);
  validateNumericalEquals(targetValues[0], 1.0);
  validateNumericalEquals(targetValues[1], 1.5);
  validateNumericalEquals(targetValues[2], 2.0);

  scale.performAction(0.0, 0.5, 1.0, 1.0);
  validateNumericalEquals(sourceValues[0], 2.0);
  validateNumericalEquals(sourceValues[1], 3.0);
  validateNumericalEquals(sourceValues[2], 4.0);
  validateNumericalEquals(targetValues[0], 2.0);
  validateNumericalEquals(targetValues[1], 3.0);
  validateNumericalEquals(targetValues[2], 4.0);
}

void ScaleActionTest:: testConfiguration()
{
  preciceTrace("testConfiguration()");
  {
    preciceDebug("Test 1");
    std::string filename = utils::Globals::getPathToSources() +
                           "/action/tests/ScaleActionTest-testConfiguration-1.xml";
    utils::XMLTag tag = utils::getRootTag();
    mesh::PtrDataConfiguration dataConfig(new mesh::DataConfiguration(tag));
    dataConfig->setDimensions(3);
    mesh::PtrMeshConfiguration meshConfig(new mesh::MeshConfiguration(tag,dataConfig));
    meshConfig->setDimensions(3);
    action::ActionConfiguration config(tag, meshConfig);
    utils::configure(tag, filename);
    validateEquals(config.actions().size(), 1);
    action::PtrAction action = config.actions().front();
    validate(action.get() != NULL);
  }
  {
    preciceDebug("Test 2");
    std::string filename = utils::Globals::getPathToSources() +
                           "/action/tests/ScaleActionTest-testConfiguration-2.xml";
    utils::XMLTag tag = utils::getRootTag();
    mesh::PtrDataConfiguration dataConfig(new mesh::DataConfiguration(tag));
    dataConfig->setDimensions(3);
    mesh::PtrMeshConfiguration meshConfig(new mesh::MeshConfiguration(tag,dataConfig));
    meshConfig->setDimensions(3);
    action::ActionConfiguration config(tag, meshConfig);
    utils::configure(tag, filename);
    validateEquals(config.actions().size(), 1);
    action::PtrAction action = config.actions().front();
    validate(action.get() != NULL);
  }
  {
    preciceDebug("Test 3");
    std::string filename = utils::Globals::getPathToSources() +
                           "/action/tests/ScaleActionTest-testConfiguration-3.xml";
    utils::XMLTag tag = utils::getRootTag();
    mesh::PtrDataConfiguration dataConfig(new mesh::DataConfiguration(tag));
    dataConfig->setDimensions(3);
    mesh::PtrMeshConfiguration meshConfig(new mesh::MeshConfiguration(tag,dataConfig));
    meshConfig->setDimensions(3);
    action::ActionConfiguration config(tag, meshConfig);
    utils::configure(tag, filename);
    validateEquals(config.actions().size(), 1);
    action::PtrAction action = config.actions().front();
    validate(action.get() != NULL);
  }
  {
    preciceDebug("Test 4");
    std::string filename = utils::Globals::getPathToSources() +
                           "/action/tests/ScaleActionTest-testConfiguration-4.xml";
    utils::XMLTag tag = utils::getRootTag();
    mesh::PtrDataConfiguration dataConfig(new mesh::DataConfiguration(tag));
    dataConfig->setDimensions(3);
    mesh::PtrMeshConfiguration meshConfig(new mesh::MeshConfiguration(tag,dataConfig));
    meshConfig->setDimensions(3);
    action::ActionConfiguration config(tag, meshConfig);
    utils::configure(tag, filename);
    validateEquals(config.actions().size(), 1);
    action::PtrAction action = config.actions().front();
    validate(action.get() != NULL);
  }
}

}}} // namespace precice, action, tests

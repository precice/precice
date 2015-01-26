// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "PythonActionTest.hpp"
#include "action/PythonAction.hpp"
#include "utils/Globals.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Data.hpp"

#ifndef PRECICE_NO_PYTHON
#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::action::tests::PythonActionTest)
#endif // not PRECICE_NO_PYTHON

namespace precice {
namespace action {
namespace tests {

tarch::logging::Log PythonActionTest::
  _log("precice::action::tests::PythonActionTest");

PythonActionTest:: PythonActionTest()
:
  TestCase("precice::action::tests::PythonActionTest")
{}

void PythonActionTest:: run()
{
  PRECICE_MASTER_ONLY{
    testMethod(testAllMethods);
    testMethod(testOmitMethods);
  }
}

void PythonActionTest:: testAllMethods()
{
  preciceTrace("testAllMethods()");
  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 3, false));
  mesh->createVertex(utils::Vector3D(1.0));
  mesh->createVertex(utils::Vector3D(2.0));
  mesh->createVertex(utils::Vector3D(3.0));
  int targetID = mesh->createData("TargetData",1)->getID();
  int sourceID = mesh->createData("SourceData",1)->getID();
  mesh->allocateDataValues();
  std::string path = utils::Globals::getPathToSources() + "/action/tests/";
  PythonAction action(PythonAction::ALWAYS_PRIOR, path, "TestAllAction", mesh,
                      targetID, sourceID);
  assignList(mesh->data(sourceID)->values()) = 0.1, 0.2, 0.3;
  assign(mesh->data(targetID)->values()) = 0.0;
  action.performAction(0.0, 0.0, 0.0, 0.0);
  tarch::la::Vector<3,double> result(2.1, 3.2, 4.3);
  validateWithMessage(tarch::la::equals(mesh->data(targetID)->values(), result),
                      mesh->data(targetID)->values());
  assign(mesh->data(sourceID)->values()) = 0.0;
  assignList(result) = 1.0, 2.0, 3.0;
  action.performAction(0.0, 0.0, 0.0, 0.0);
  validateWithMessage(tarch::la::equals(mesh->data(targetID)->values(), result),
                      mesh->data(targetID)->values());
}

void PythonActionTest:: testOmitMethods()
{
  preciceTrace("testOmitMethods()");
  std::string path = utils::Globals::getPathToSources() + "/action/tests/";
  {
    preciceDebug("Test 1");
    mesh::PtrMesh mesh;
    PythonAction action(PythonAction::ALWAYS_PRIOR, path, "TestOmitAction1",
                        mesh, -1, -1);
    action.performAction(0.0, 0.0, 0.0, 0.0);
  }
  {
    preciceDebug("Test 2");
    mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 3, false));
    mesh->createVertex(utils::Vector3D(0.0));
    mesh::PtrData data = mesh->createData("TargetData", 1);
    mesh->allocateDataValues();
    PythonAction action(PythonAction::ALWAYS_PRIOR, path, "TestOmitAction2",
                        mesh, data->getID(), -1);
    action.performAction(0.0, 0.0, 0.0, 0.0);
  }
  {
    preciceDebug("Test 3");
    mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 3, false));
    mesh->createVertex(utils::Vector3D(0.0));
    mesh::PtrData data = mesh->createData("SourceData", 1);
    mesh->allocateDataValues();
    PythonAction action(PythonAction::ALWAYS_PRIOR, path, "TestOmitAction3",
                        mesh, -1, data->getID());
    action.performAction(0.0, 0.0, 0.0, 0.0);
  }
}

}}} // namespace precice, action, tests

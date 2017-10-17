#include "PythonActionTest.hpp"
#include "action/PythonAction.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Data.hpp"
#include "math/math.hpp"
#include "utils/Globals.hpp"

#ifndef PRECICE_NO_PYTHON
#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::action::tests::PythonActionTest)
#endif // not PRECICE_NO_PYTHON

namespace precice {
namespace action {
namespace tests {

logging::Logger PythonActionTest::_log("action::tests::PythonActionTest");

PythonActionTest:: PythonActionTest()
:
  TestCase("action::tests::PythonActionTest")
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
  TRACE();
  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 3, false));
  mesh->createVertex(Eigen::Vector3d::Constant(1.0));
  mesh->createVertex(Eigen::Vector3d::Constant(2.0));
  mesh->createVertex(Eigen::Vector3d::Constant(3.0));
  int targetID = mesh->createData("TargetData",1)->getID();
  int sourceID = mesh->createData("SourceData",1)->getID();
  mesh->allocateDataValues();
  std::string path = utils::getPathToSources() + "/action/tests/";
  PythonAction action(PythonAction::ALWAYS_PRIOR, path, "TestAllAction", mesh,
                      targetID, sourceID);
  mesh->data(sourceID)->values() << 0.1, 0.2, 0.3;
  mesh->data(targetID)->values() = Eigen::VectorXd::Zero(mesh->data(targetID)->values().size());
  action.performAction(0.0, 0.0, 0.0, 0.0);
  Eigen::Vector3d result(2.1, 3.2, 4.3);
  validateWithMessage(math::equals(mesh->data(targetID)->values(), result), mesh->data(targetID)->values());
  mesh->data(sourceID)->values() = Eigen::VectorXd::Zero(mesh->data(sourceID)->values().size());
  result << 1.0, 2.0, 3.0;
  action.performAction(0.0, 0.0, 0.0, 0.0);
  validateWithMessage(math::equals(mesh->data(targetID)->values(), result), mesh->data(targetID)->values());
}

void PythonActionTest:: testOmitMethods()
{
  TRACE();
  std::string path = utils::getPathToSources() + "/action/tests/";
  {
    DEBUG("Test 1");
    mesh::PtrMesh mesh;
    PythonAction action(PythonAction::ALWAYS_PRIOR, path, "TestOmitAction1",
                        mesh, -1, -1);
    action.performAction(0.0, 0.0, 0.0, 0.0);
  }
  {
    DEBUG("Test 2");
    mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 3, false));
    mesh->createVertex(Eigen::Vector3d::Zero());
    mesh::PtrData data = mesh->createData("TargetData", 1);
    mesh->allocateDataValues();
    PythonAction action(PythonAction::ALWAYS_PRIOR, path, "TestOmitAction2",
                        mesh, data->getID(), -1);
    action.performAction(0.0, 0.0, 0.0, 0.0);
  }
  {
    DEBUG("Test 3");
    mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 3, false));
    mesh->createVertex(Eigen::Vector3d::Zero());
    mesh::PtrData data = mesh->createData("SourceData", 1);
    mesh->allocateDataValues();
    PythonAction action(PythonAction::ALWAYS_PRIOR, path, "TestOmitAction3",
                        mesh, -1, data->getID());
    action.performAction(0.0, 0.0, 0.0, 0.0);
  }
}

}}} // namespace precice, action, tests

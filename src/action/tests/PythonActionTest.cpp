#ifndef PRECICE_NO_PYTHON

#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include <string>
#include "action/Action.hpp"
#include "action/PythonAction.hpp"
#include "logging/Logger.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::action;

BOOST_AUTO_TEST_SUITE(ActionTests)
BOOST_AUTO_TEST_SUITE(Python)

BOOST_AUTO_TEST_CASE(AllMethods)
{
  PRECICE_TEST(1_rank);
  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 3, false, testing::nextMeshID()));
  mesh->createVertex(Eigen::Vector3d::Constant(1.0));
  mesh->createVertex(Eigen::Vector3d::Constant(2.0));
  mesh->createVertex(Eigen::Vector3d::Constant(3.0));
  int targetID = mesh->createData("TargetData", 1)->getID();
  int sourceID = mesh->createData("SourceData", 1)->getID();
  mesh->allocateDataValues();
  std::string  path = testing::getPathToSources() + "/action/tests/";
  PythonAction action(PythonAction::WRITE_MAPPING_PRIOR, path, "TestAllAction", mesh, targetID, sourceID);
  mesh->data(sourceID)->values() << 0.1, 0.2, 0.3;
  mesh->data(targetID)->values() = Eigen::VectorXd::Zero(mesh->data(targetID)->values().size());
  action.performAction(0.0, 0.0, 0.0, 0.0);
  Eigen::Vector3d result(2.1, 3.2, 4.3);
  BOOST_TEST(testing::equals(mesh->data(targetID)->values(), result));
  mesh->data(sourceID)->values() = Eigen::VectorXd::Zero(mesh->data(sourceID)->values().size());
  result << 1.0, 2.0, 3.0;
  action.performAction(0.0, 0.0, 0.0, 0.0);
  BOOST_TEST(testing::equals(mesh->data(targetID)->values(), result));
}

BOOST_AUTO_TEST_CASE(OmitMethods)
{
  PRECICE_TEST(1_rank);
  std::string path = testing::getPathToSources() + "/action/tests/";
  {
    mesh::PtrMesh mesh;
    PythonAction  action(PythonAction::WRITE_MAPPING_PRIOR, path, "TestOmitAction1", mesh, -1, -1);
    action.performAction(0.0, 0.0, 0.0, 0.0);
  }
  {
    mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 3, false, testing::nextMeshID()));
    mesh->createVertex(Eigen::Vector3d::Zero());
    mesh::PtrData data = mesh->createData("TargetData", 1);
    mesh->allocateDataValues();
    PythonAction action(PythonAction::WRITE_MAPPING_PRIOR, path, "TestOmitAction2", mesh, data->getID(), -1);
    action.performAction(0.0, 0.0, 0.0, 0.0);
  }
  {
    mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 3, false, testing::nextMeshID()));
    mesh->createVertex(Eigen::Vector3d::Zero());
    mesh::PtrData data = mesh->createData("SourceData", 1);
    mesh->allocateDataValues();
    PythonAction action(PythonAction::WRITE_MAPPING_PRIOR, path, "TestOmitAction3", mesh, -1, data->getID());
    action.performAction(0.0, 0.0, 0.0, 0.0);
  }
}

BOOST_AUTO_TEST_SUITE_END() // Python
BOOST_AUTO_TEST_SUITE_END() // ActionTest

#endif // not PRECICE_NO_PYTHON

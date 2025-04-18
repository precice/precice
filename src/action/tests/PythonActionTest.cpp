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

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(PerformActionWithGlobalIterationsCounter)
{
  PRECICE_TEST();
  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 3, testing::nextMeshID()));
  mesh->createVertex(Eigen::Vector3d::Constant(1.0));
  mesh->createVertex(Eigen::Vector3d::Constant(2.0));
  mesh->createVertex(Eigen::Vector3d::Constant(3.0));
  int targetID = mesh->createData("TargetData", 1, 0_dataID)->getID();
  int sourceID = mesh->createData("SourceData", 1, 1_dataID)->getID();
  mesh->allocateDataValues();
  std::string  path = testing::getPathToSources() + "/action/tests/";
  PythonAction action(PythonAction::WRITE_MAPPING_POST, path, "TestAction", mesh, targetID, sourceID);

  mesh->data(sourceID)->emplaceSampleAtTime(1, {0.1, 0.2, 0.3});
  mesh->data(targetID)->emplaceSampleAtTime(1, {0.0, 0.0, 0.0});

  action.performAction();

  Eigen::VectorXd result(3);
  result << 1.1, 1.2, 1.3;
  BOOST_TEST(testing::equals(mesh->data(targetID)->stamples().back().sample.values, result));
  mesh->data(sourceID)->emplaceSampleAtTime(1, {0.0, 0.0, 0.0});

  action.performAction();

  result << 2.0, 2.0, 2.0;
  BOOST_TEST(testing::equals(mesh->data(targetID)->stamples().back().sample.values, result));
}

BOOST_AUTO_TEST_SUITE_END() // Python
BOOST_AUTO_TEST_SUITE_END() // ActionTest

#endif // not PRECICE_NO_PYTHON

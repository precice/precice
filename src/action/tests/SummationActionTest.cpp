#include <Eigen/Core>
#include <algorithm>
#include <list>
#include <memory>
#include <string>
#include <vector>
#include "action/Action.hpp"
#include "action/SharedPointer.hpp"
#include "action/SummationAction.hpp"
#include "action/config/ActionConfiguration.hpp"
#include "logging/Logger.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "xml/XMLTag.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(ActionTests)
BOOST_AUTO_TEST_SUITE(Summation)

BOOST_AUTO_TEST_CASE(SummationOneDimensional)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;
  PtrMesh          mesh(new Mesh("Mesh", 3, testing::nextMeshID()));
  int              dimensions  = 1;
  PtrData          sourceData1 = mesh->createData("SourceData1", dimensions, 0_dataID);
  PtrData          sourceData2 = mesh->createData("SourceData2", dimensions, 1_dataID);
  PtrData          sourceData3 = mesh->createData("SourceData3", dimensions, 2_dataID);
  PtrData          targetData  = mesh->createData("TargetData", dimensions, 3_dataID);
  std::vector<int> sourceDataIDs{sourceData1->getID(), sourceData2->getID(), sourceData3->getID()};
  int              targetDataID = targetData->getID();
  mesh->createVertex(Eigen::Vector3d::Constant(0.0));
  mesh->createVertex(Eigen::Vector3d::Constant(1.0));
  mesh->createVertex(Eigen::Vector3d::Constant(2.0));

  mesh->allocateDataValues();

  Eigen::VectorXd v1(3), v2(3), v3(3);
  v1 << 2.0, 3.0, 4.0;
  v2 << 1.0, 2.0, 3.0;
  v3 << 2.0, 3.0, 4.0;
  sourceData1->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{dimensions, v1});
  sourceData2->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{dimensions, v2});
  sourceData3->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{dimensions, v3});
  // targetData->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{dimensions,Eigen::VectorXd::Zero(targetValues.size())});

  action::SummationAction sum(
      action::SummationAction::WRITE_MAPPING_POST, sourceDataIDs, targetDataID, mesh);

  sum.performAction(0.0);

  const auto &sourceValues1 = sourceData1->values();
  const auto &sourceValues2 = sourceData2->values();
  const auto &sourceValues3 = sourceData3->values();
  const auto &targetValues  = targetData->values();
  BOOST_TEST(sourceValues1(0) == 2.0);
  BOOST_TEST(sourceValues1(1) == 3.0);
  BOOST_TEST(sourceValues1(2) == 4.0);
  BOOST_TEST(sourceValues2(0) == 1.0);
  BOOST_TEST(sourceValues2(1) == 2.0);
  BOOST_TEST(sourceValues2(2) == 3.0);
  BOOST_TEST(sourceValues3(0) == 2.0);
  BOOST_TEST(sourceValues3(1) == 3.0);
  BOOST_TEST(sourceValues3(2) == 4.0);
  BOOST_TEST(targetValues(0) == 5.0);
  BOOST_TEST(targetValues(1) == 8.0);
  BOOST_TEST(targetValues(2) == 11.0);
}

BOOST_AUTO_TEST_CASE(SummationThreeDimensional)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;
  int              dimensions = 3;
  PtrMesh          mesh(new Mesh("Mesh", dimensions, testing::nextMeshID()));
  PtrData          sourceData1 = mesh->createData("SourceData1", dimensions, 0_dataID);
  PtrData          sourceData2 = mesh->createData("SourceData2", dimensions, 1_dataID);
  PtrData          targetData  = mesh->createData("TargetData", dimensions, 2_dataID);
  std::vector<int> sourceDataIDs{sourceData1->getID(), sourceData2->getID()};
  int              targetDataID = targetData->getID();
  mesh->createVertex(Eigen::Vector3d::Constant(0.0));
  mesh->createVertex(Eigen::Vector3d::Constant(1.0));
  mesh->createVertex(Eigen::Vector3d::Constant(2.0));
  mesh->allocateDataValues();

  Eigen::VectorXd v1(9), v2(9);
  v1 << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
  v2 << 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0;

  sourceData1->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{dimensions, v1});
  sourceData2->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{dimensions, v2});
  // targetData->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{dimensions,Eigen::VectorXd::Zero(targetValues.size())})

  action::SummationAction sum(
      action::SummationAction::WRITE_MAPPING_POST, sourceDataIDs, targetDataID, mesh);

  sum.performAction(0.0);

  const auto &sourceValues1 = sourceData1->values();
  const auto &sourceValues2 = sourceData2->values();
  const auto &targetValues  = targetData->values();
  BOOST_TEST(sourceValues1(0) == 1.0);
  BOOST_TEST(sourceValues2(0) == 2.0);
  BOOST_TEST(targetValues(0) == 3.0);

  BOOST_TEST(sourceValues1(1) == 2.0);
  BOOST_TEST(sourceValues2(1) == 3.0);
  BOOST_TEST(targetValues(1) == 5.0);

  BOOST_TEST(sourceValues1(2) == 3.0);
  BOOST_TEST(sourceValues2(2) == 4.0);
  BOOST_TEST(targetValues(2) == 7.0);

  BOOST_TEST(sourceValues1(3) == 4.0);
  BOOST_TEST(sourceValues2(3) == 5.0);
  BOOST_TEST(targetValues(3) == 9.0);

  BOOST_TEST(sourceValues1(4) == 5.0);
  BOOST_TEST(sourceValues2(4) == 6.0);
  BOOST_TEST(targetValues(4) == 11.0);

  BOOST_TEST(sourceValues1(5) == 6.0);
  BOOST_TEST(sourceValues2(5) == 7.0);
  BOOST_TEST(targetValues(5) == 13.0);

  BOOST_TEST(sourceValues1(6) == 7.0);
  BOOST_TEST(sourceValues2(6) == 8.0);
  BOOST_TEST(targetValues(6) == 15.0);

  BOOST_TEST(sourceValues1(7) == 8.0);
  BOOST_TEST(sourceValues2(7) == 9.0);
  BOOST_TEST(targetValues(7) == 17.0);

  BOOST_TEST(sourceValues1(8) == 9.0);
  BOOST_TEST(sourceValues2(8) == 10.0);
  BOOST_TEST(targetValues(8) == 19.0);
}

BOOST_AUTO_TEST_CASE(SummationThreeDimensionalSubcycling)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;
  int              dimensions = 3;
  PtrMesh          mesh(new Mesh("Mesh", dimensions, testing::nextMeshID()));
  PtrData          sourceData1 = mesh->createData("SourceData1", dimensions, 0_dataID);
  PtrData          sourceData2 = mesh->createData("SourceData2", dimensions, 1_dataID);
  PtrData          targetData  = mesh->createData("TargetData", dimensions, 2_dataID);
  std::vector<int> sourceDataIDs{sourceData1->getID(), sourceData2->getID()};
  int              targetDataID = targetData->getID();
  mesh->createVertex(Eigen::Vector3d::Constant(0.0));
  mesh->createVertex(Eigen::Vector3d::Constant(1.0));
  mesh->createVertex(Eigen::Vector3d::Constant(2.0));
  mesh->allocateDataValues();

  Eigen::VectorXd v1_05(9), v2_05(9), v1_1(9), v2_1(9);
  v1_05 << 11.0, 22.0, 33.0, 44.0, 55.0, 66.0, 77.0, 88.0, 99.0;
  v2_05 << 12.0, 23.0, 34.0, 45.0, 56.0, 67.0, 78.0, 89.0, 910.0;
  v1_1 << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
  v2_1 << 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0;

  sourceData1->setSampleAtTime(time::Storage::WINDOW_END * 0.5, time::Sample{dimensions, v1_05});
  sourceData2->setSampleAtTime(time::Storage::WINDOW_END * 0.5, time::Sample{dimensions, v2_05});
  sourceData1->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{dimensions, v1_1});
  sourceData2->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{dimensions, v2_1});

  action::SummationAction sum(
      action::SummationAction::WRITE_MAPPING_POST, sourceDataIDs, targetDataID, mesh);

  sum.performAction(0.0);

  // By default data buffers after action hold samples from WINDOW_END

  const auto &sourceValues1 = sourceData1->values();
  const auto &sourceValues2 = sourceData2->values();
  const auto &targetValues  = targetData->values();
  BOOST_TEST(sourceValues1(0) == 1.0);
  BOOST_TEST(sourceValues2(0) == 2.0);
  BOOST_TEST(targetValues(0) == 3.0);

  BOOST_TEST(sourceValues1(1) == 2.0);
  BOOST_TEST(sourceValues2(1) == 3.0);
  BOOST_TEST(targetValues(1) == 5.0);

  BOOST_TEST(sourceValues1(2) == 3.0);
  BOOST_TEST(sourceValues2(2) == 4.0);
  BOOST_TEST(targetValues(2) == 7.0);

  BOOST_TEST(sourceValues1(3) == 4.0);
  BOOST_TEST(sourceValues2(3) == 5.0);
  BOOST_TEST(targetValues(3) == 9.0);

  BOOST_TEST(sourceValues1(4) == 5.0);
  BOOST_TEST(sourceValues2(4) == 6.0);
  BOOST_TEST(targetValues(4) == 11.0);

  BOOST_TEST(sourceValues1(5) == 6.0);
  BOOST_TEST(sourceValues2(5) == 7.0);
  BOOST_TEST(targetValues(5) == 13.0);

  BOOST_TEST(sourceValues1(6) == 7.0);
  BOOST_TEST(sourceValues2(6) == 8.0);
  BOOST_TEST(targetValues(6) == 15.0);

  BOOST_TEST(sourceValues1(7) == 8.0);
  BOOST_TEST(sourceValues2(7) == 9.0);
  BOOST_TEST(targetValues(7) == 17.0);

  BOOST_TEST(sourceValues1(8) == 9.0);
  BOOST_TEST(sourceValues2(8) == 10.0);
  BOOST_TEST(targetValues(8) == 19.0);

  // Load and check data from 0.5 * time::Storage::WINDOW_END

  auto &loadedStample1 = sourceData1->stamples().front();
  BOOST_TEST(loadedStample1.timestamp == time::Storage::WINDOW_END * 0.5);
  sourceData1->values() = loadedStample1.sample.values;

  auto &loadedStample2 = sourceData2->stamples().front();
  BOOST_TEST(loadedStample2.timestamp == time::Storage::WINDOW_END * 0.5);
  sourceData2->values() = loadedStample2.sample.values;

  auto &loadedStample3 = targetData->stamples().front();
  BOOST_TEST(loadedStample3.timestamp == time::Storage::WINDOW_END * 0.5);
  targetData->values() = loadedStample3.sample.values;

  BOOST_TEST(sourceValues1(0) == 11.0);
  BOOST_TEST(sourceValues2(0) == 12.0);
  BOOST_TEST(targetValues(0) == 23.0);

  BOOST_TEST(sourceValues1(1) == 22.0);
  BOOST_TEST(sourceValues2(1) == 23.0);
  BOOST_TEST(targetValues(1) == 45.0);

  BOOST_TEST(sourceValues1(2) == 33.0);
  BOOST_TEST(sourceValues2(2) == 34.0);
  BOOST_TEST(targetValues(2) == 67.0);

  BOOST_TEST(sourceValues1(3) == 44.0);
  BOOST_TEST(sourceValues2(3) == 45.0);
  BOOST_TEST(targetValues(3) == 89.0);

  BOOST_TEST(sourceValues1(4) == 55.0);
  BOOST_TEST(sourceValues2(4) == 56.0);
  BOOST_TEST(targetValues(4) == 111.0);

  BOOST_TEST(sourceValues1(5) == 66.0);
  BOOST_TEST(sourceValues2(5) == 67.0);
  BOOST_TEST(targetValues(5) == 133.0);

  BOOST_TEST(sourceValues1(6) == 77.0);
  BOOST_TEST(sourceValues2(6) == 78.0);
  BOOST_TEST(targetValues(6) == 155.0);

  BOOST_TEST(sourceValues1(7) == 88.0);
  BOOST_TEST(sourceValues2(7) == 89.0);
  BOOST_TEST(targetValues(7) == 177.0);

  BOOST_TEST(sourceValues1(8) == 99.0);
  BOOST_TEST(sourceValues2(8) == 910.0);
  BOOST_TEST(targetValues(8) == 1009.0);
}

BOOST_AUTO_TEST_CASE(Configuration)
{
  PRECICE_TEST(1_rank);
  std::string                filename = testing::getPathToSources() + "/action/tests/SummationActionTest-testConfiguration-1.xml";
  xml::XMLTag                tag      = xml::getRootTag();
  mesh::PtrDataConfiguration dataConfig(new mesh::DataConfiguration(tag));
  dataConfig->setDimensions(3);
  mesh::PtrMeshConfiguration meshConfig(new mesh::MeshConfiguration(tag, dataConfig));
  meshConfig->setDimensions(3);
  action::ActionConfiguration config(tag, meshConfig);
  xml::configure(tag, xml::ConfigurationContext{}, filename);
  BOOST_TEST(config.actions().size() == 1);
  auto &action = config.actions().front();
  BOOST_TEST(static_cast<bool>(action));
}

BOOST_AUTO_TEST_SUITE_END() // Summation
BOOST_AUTO_TEST_SUITE_END() // ActionTest

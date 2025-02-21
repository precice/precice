#include "com/SerializedStamples.hpp"
#include "com/tests/helper.hpp"
#include "cplscheme/CouplingData.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/EigenHelperFunctions.hpp"

using namespace precice;
using namespace precice::com;
using precice::testing::makeCouplingData;

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(SerializedStamples)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(SerializeValues)
{
  PRECICE_TEST();
  std::vector<int> vertexOffsets{4, 8, 8, 10};

  const int meshDimensions = 3;
  const int dataDimensions = 1;
  const int nValues        = 4;
  const int nTimeSteps     = 3;

  mesh::PtrMesh dummyMesh(new mesh::Mesh("DummyMesh", 3, testing::nextMeshID()));
  dummyMesh->setVertexOffsets(vertexOffsets);
  dummyMesh->createVertex(Eigen::Vector3d{0, 0, 0});
  dummyMesh->createVertex(Eigen::Vector3d{1, 0, 0});
  dummyMesh->createVertex(Eigen::Vector3d{2, 0, 0});
  dummyMesh->createVertex(Eigen::Vector3d{3, 0, 0});

  mesh::PtrData fromData(new mesh::Data("from", -1, dataDimensions));
  mesh::PtrData toData(new mesh::Data("to", -1, dataDimensions));

  Eigen::VectorXd insert0(nValues);
  insert0 << 1.0, 2.0, 3.0, 4.0;
  Eigen::VectorXd insert05(nValues);
  insert05 << 10.0, 20.0, 30.0, 40.0;
  Eigen::VectorXd insert1(nValues);
  insert1 << 100.0, 200.0, 300.0, 400.0;

  cplscheme::PtrCouplingData fromDataPtr = makeCouplingData(fromData, dummyMesh);
  cplscheme::PtrCouplingData toDataPtr   = makeCouplingData(toData, dummyMesh);

  fromDataPtr->setSampleAtTime(0, time::Sample{dataDimensions, insert0});
  fromDataPtr->setSampleAtTime(0.5, time::Sample{dataDimensions, insert05});
  fromDataPtr->setSampleAtTime(1, time::Sample{dataDimensions, insert1});

  const auto serialized = serialize::SerializedStamples::serialize(fromDataPtr);

  Eigen::VectorXd expectedSerialized(nTimeSteps * nValues);
  expectedSerialized << 1.0, 10.0, 100.0, 2.0, 20.0, 200.0, 3.0, 30.0, 300.0, 4.0, 40.0, 400.0;

  for (int i = 0; i < nTimeSteps * nValues; i++) {
    BOOST_TEST(testing::equals(serialized.values()(i), expectedSerialized(i)));
  }
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(DeserializeValues)
{
  PRECICE_TEST();
  std::vector<int> vertexOffsets{4, 8, 8, 10};

  const int nValues        = 4;
  const int dataDimensions = 1;
  const int nTimeSteps     = 3;

  Eigen::VectorXd serializedValues(nTimeSteps * nValues);
  serializedValues << 1.0, 10.0, 100.0, 2.0, 20.0, 200.0, 3.0, 30.0, 300.0, 4.0, 40.0, 400.0;

  mesh::PtrMesh dummyMesh(new mesh::Mesh("DummyMesh", 3, testing::nextMeshID()));
  dummyMesh->setVertexOffsets(vertexOffsets);
  dummyMesh->createVertex(Eigen::Vector3d{0, 0, 0});
  dummyMesh->createVertex(Eigen::Vector3d{1, 0, 0});
  dummyMesh->createVertex(Eigen::Vector3d{2, 0, 0});
  dummyMesh->createVertex(Eigen::Vector3d{3, 0, 0});

  mesh::PtrData              toData(new mesh::Data("to", -1, dataDimensions));
  cplscheme::PtrCouplingData toDataPtr = makeCouplingData(toData, dummyMesh);

  Eigen::VectorXd initValues(nValues);
  initValues.setZero();

  toDataPtr->setSampleAtTime(0, time::Sample(toDataPtr->getDimensions(), initValues));

  Eigen::VectorXd timeStamps(nTimeSteps);
  timeStamps << 0, 0.5, 1;

  auto serialized = serialize::SerializedStamples::empty(timeStamps, toDataPtr);

  serialized.values() = serializedValues;

  serialized.deserializeInto(timeStamps, toDataPtr);

  std::vector<Eigen::VectorXd> expectedValues;

  Eigen::VectorXd insert0(nValues);
  insert0 << 1.0, 2.0, 3.0, 4.0;
  Eigen::VectorXd insert05(nValues);
  insert05 << 10.0, 20.0, 30.0, 40.0;
  Eigen::VectorXd insert1(nValues);
  insert1 << 100.0, 200.0, 300.0, 400.0;

  expectedValues.push_back(insert0);
  expectedValues.push_back(insert05);
  expectedValues.push_back(insert1);

  BOOST_TEST(toDataPtr->timeStepsStorage().stamples().size() == nTimeSteps);
  for (int t = 0; t < nTimeSteps; t++) {
    BOOST_TEST(testing::equals(toDataPtr->timeStepsStorage().stamples()[t].timestamp, timeStamps[t]));
    for (int i = 0; i < nValues; i++) {
      BOOST_TEST(testing::equals(toDataPtr->timeStepsStorage().stamples()[t].sample.values(i), expectedValues[t](i)));
    }
  }
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(SerializeValuesAndGradients)
{
  PRECICE_TEST();
  std::vector<int> vertexOffsets{4, 8, 8, 10};

  const int meshDimensions = 3;
  const int dataDimensions = 1;
  const int nValues        = 4;
  const int nTimeSteps     = 3;

  mesh::PtrMesh dummyMesh(new mesh::Mesh("DummyMesh", meshDimensions, testing::nextMeshID()));
  dummyMesh->setVertexOffsets(vertexOffsets);
  dummyMesh->createVertex(Eigen::Vector3d{0, 0, 0});
  dummyMesh->createVertex(Eigen::Vector3d{1, 0, 0});
  dummyMesh->createVertex(Eigen::Vector3d{2, 0, 0});
  dummyMesh->createVertex(Eigen::Vector3d{3, 0, 0});

  mesh::PtrData fromData(new mesh::Data("from", -1, dataDimensions));
  fromData->requireDataGradient();

  Eigen::VectorXd insert0(nValues);
  insert0 << 1.0, 2.0, 3.0, 4.0;
  Eigen::VectorXd insert05(nValues);
  insert05 << 10.0, 20.0, 30.0, 40.0;
  Eigen::VectorXd insert1(nValues);
  insert1 << 100.0, 200.0, 300.0, 400.0;
  Eigen::MatrixXd insertGradients0(meshDimensions, nValues * dataDimensions);
  insertGradients0 << 10.0, 20.0, 30.0, 40.0,
      11.0, 21.0, 31.0, 41.0,
      12.0, 22.0, 32.0, 42.0;
  Eigen::MatrixXd insertGradients05(meshDimensions, nValues * dataDimensions);
  insertGradients05 << 100.0, 200.0, 300.0, 400.0,
      101.0, 201.0, 301.0, 401.0,
      102.0, 202.0, 302.0, 402.0;
  Eigen::MatrixXd insertGradients1(meshDimensions, nValues * dataDimensions);
  insertGradients1 << 1000.0, 2000.0, 3000.0, 4000.0,
      1001.0, 2001.0, 3001.0, 4001.0,
      1002.0, 2002.0, 3002.0, 4002.0;

  cplscheme::PtrCouplingData fromDataPtr = makeCouplingData(fromData, dummyMesh);

  fromDataPtr->setSampleAtTime(0, time::Sample{dataDimensions, insert0, insertGradients0});
  fromDataPtr->setSampleAtTime(0.5, time::Sample{dataDimensions, insert05, insertGradients05});
  fromDataPtr->setSampleAtTime(1, time::Sample{dataDimensions, insert1, insertGradients1});

  const auto serialized = serialize::SerializedStamples::serialize(fromDataPtr);

  Eigen::VectorXd expectedSerializedValues(nTimeSteps * nValues);
  expectedSerializedValues << 1.0, 10.0, 100.0, 2.0, 20.0, 200.0, 3.0, 30.0, 300.0, 4.0, 40.0, 400.0;

  Eigen::VectorXd expectedSerializedGradients(nTimeSteps * nValues * meshDimensions);
  expectedSerializedGradients << 10.0, 100.0, 1000.0, 11.0, 101.0, 1001.0, 12.0, 102.0, 1002.0,
      20.0, 200.0, 2000.0, 21.0, 201.0, 2001.0, 22.0, 202.0, 2002.0,
      30.0, 300.0, 3000.0, 31.0, 301.0, 3001.0, 32.0, 302.0, 3002.0,
      40.0, 400.0, 4000.0, 41.0, 401.0, 4001.0, 42.0, 402.0, 4002.0;

  for (int i = 0; i < nTimeSteps * nValues; i++) {
    BOOST_TEST(testing::equals(serialized.values()(i), expectedSerializedValues(i)));
  }

  for (int i = 0; i < nTimeSteps * meshDimensions * nValues; i++) {
    BOOST_TEST(testing::equals(serialized.gradients()(i), expectedSerializedGradients(i)));
  }
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(DeserializeValuesAndGradients)
{
  PRECICE_TEST();
  std::vector<int> vertexOffsets{4, 8, 8, 10};

  const int meshDimensions = 3;
  const int nValues        = 4;
  const int dataDimensions = 1;
  const int nTimeSteps     = 3;

  Eigen::VectorXd serializedValues(nTimeSteps * nValues);
  serializedValues << 1.0, 10.0, 100.0, 2.0, 20.0, 200.0, 3.0, 30.0, 300.0, 4.0, 40.0, 400.0;

  Eigen::VectorXd serializedGradients(nTimeSteps * nValues * meshDimensions);
  serializedGradients << 10.0, 100.0, 1000.0, 11.0, 101.0, 1001.0, 12.0, 102.0, 1002.0,
      20.0, 200.0, 2000.0, 21.0, 201.0, 2001.0, 22.0, 202.0, 2002.0,
      30.0, 300.0, 3000.0, 31.0, 301.0, 3001.0, 32.0, 302.0, 3002.0,
      40.0, 400.0, 4000.0, 41.0, 401.0, 4001.0, 42.0, 402.0, 4002.0;

  mesh::PtrMesh dummyMesh(new mesh::Mesh("DummyMesh", 3, testing::nextMeshID()));
  dummyMesh->setVertexOffsets(vertexOffsets);
  dummyMesh->createVertex(Eigen::Vector3d{0, 0, 0});
  dummyMesh->createVertex(Eigen::Vector3d{1, 0, 0});
  dummyMesh->createVertex(Eigen::Vector3d{2, 0, 0});
  dummyMesh->createVertex(Eigen::Vector3d{3, 0, 0});

  mesh::PtrData toData(new mesh::Data("to", -1, dataDimensions));
  toData->requireDataGradient();

  cplscheme::PtrCouplingData toDataPtr = makeCouplingData(toData, dummyMesh);

  Eigen::VectorXd initValues(nValues);
  initValues.setZero();

  Eigen::MatrixXd initGradients(meshDimensions, nValues);
  initGradients.setZero();

  toDataPtr->setSampleAtTime(0, time::Sample(toDataPtr->getDimensions(), initValues, initGradients));

  Eigen::VectorXd timeStamps(nTimeSteps);
  timeStamps << 0, 0.5, 1;

  auto serialized = serialize::SerializedStamples::empty(timeStamps, toDataPtr);

  serialized.values()    = serializedValues;
  serialized.gradients() = serializedGradients;

  serialized.deserializeInto(timeStamps, toDataPtr);

  std::vector<Eigen::VectorXd> expectedValues;

  Eigen::VectorXd insert0(nValues);
  insert0 << 1.0, 2.0, 3.0, 4.0;
  Eigen::VectorXd insert05(nValues);
  insert05 << 10.0, 20.0, 30.0, 40.0;
  Eigen::VectorXd insert1(nValues);
  insert1 << 100.0, 200.0, 300.0, 400.0;

  expectedValues.push_back(insert0);
  expectedValues.push_back(insert05);
  expectedValues.push_back(insert1);

  std::vector<Eigen::MatrixXd> expectedGradients;

  Eigen::MatrixXd insertGradients0(meshDimensions, dataDimensions * nValues);
  insertGradients0 << 10.0, 20.0, 30.0, 40.0,
      11.0, 21.0, 31.0, 41.0,
      12.0, 22.0, 32.0, 42.0;
  Eigen::MatrixXd insertGradients05(meshDimensions, nValues * dataDimensions);
  insertGradients05 << 100.0, 200.0, 300.0, 400.0,
      101.0, 201.0, 301.0, 401.0,
      102.0, 202.0, 302.0, 402.0;
  Eigen::MatrixXd insertGradients1(meshDimensions, nValues * dataDimensions);
  insertGradients1 << 1000.0, 2000.0, 3000.0, 4000.0,
      1001.0, 2001.0, 3001.0, 4001.0,
      1002.0, 2002.0, 3002.0, 4002.0;

  expectedGradients.push_back(insertGradients0);
  expectedGradients.push_back(insertGradients05);
  expectedGradients.push_back(insertGradients1);

  BOOST_TEST(toDataPtr->timeStepsStorage().stamples().size() == nTimeSteps);
  for (int t = 0; t < nTimeSteps; t++) {
    BOOST_TEST(testing::equals(toDataPtr->timeStepsStorage().stamples()[t].timestamp, timeStamps[t]));
    for (int i = 0; i < nValues; i++) {
      BOOST_TEST(testing::equals(toDataPtr->timeStepsStorage().stamples()[t].sample.values(i), expectedValues[t](i)));
    }
    for (int i = 0; i < meshDimensions * nValues; i++) {
      BOOST_TEST(testing::equals(toDataPtr->timeStepsStorage().stamples()[t].sample.gradients(i), expectedGradients[t](i)));
    }
  }
}

BOOST_AUTO_TEST_SUITE_END() // SerializedStamples
BOOST_AUTO_TEST_SUITE_END() // CommunicationTests

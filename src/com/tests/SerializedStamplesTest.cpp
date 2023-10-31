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

BOOST_AUTO_TEST_CASE(SerializeValues)
{
  std::vector<int> vertexOffsets{4, 8, 8, 10};

  const int meshDimensions = 3;
  const int dataDimensions = 1;
  const int nValues        = 4;
  const int nTimeSteps     = 3;

  mesh::PtrMesh dummyMesh(new mesh::Mesh("DummyMesh", 3, testing::nextMeshID()));
  dummyMesh->setVertexOffsets(vertexOffsets);

  mesh::PtrData fromData(new mesh::Data("from", -1, dataDimensions));
  mesh::PtrData toData(new mesh::Data("to", -1, dataDimensions));

  Eigen::VectorXd insert0(nValues);
  insert0 << 1.0, 2.0, 3.0, 4.0;
  Eigen::VectorXd insert05(nValues);
  insert05 << 10.0, 20.0, 30.0, 40.0;
  Eigen::VectorXd insert1(nValues);
  insert1 << 100.0, 200.0, 300.0, 400.0;

  cplscheme::PtrCouplingData fromDataPtr = makeCouplingData(fromData, dummyMesh, true);
  cplscheme::PtrCouplingData toDataPtr   = makeCouplingData(toData, dummyMesh, true);

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

BOOST_AUTO_TEST_CASE(DeserializeValues)
{
  std::vector<int> vertexOffsets{4, 8, 8, 10};

  const int meshDimensions = 3;
  const int nValues        = 4;
  const int nTimeSteps     = 3;

  Eigen::VectorXd serializedValues(nTimeSteps * nValues);
  serializedValues << 1.0, 10.0, 100.0, 2.0, 20.0, 200.0, 3.0, 30.0, 300.0, 4.0, 40.0, 400.0;

  mesh::PtrMesh dummyMesh(new mesh::Mesh("DummyMesh", 3, testing::nextMeshID()));
  dummyMesh->setVertexOffsets(vertexOffsets);

  mesh::PtrData              toData(new mesh::Data("to", -1, 1));
  cplscheme::PtrCouplingData toDataPtr = makeCouplingData(toData, dummyMesh, true);
  toDataPtr->sample().values           = Eigen::VectorXd(nValues);
  toDataPtr->sample().setZero();

  toDataPtr->setSampleAtTime(0, toDataPtr->sample());

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

BOOST_AUTO_TEST_CASE(SerializeValuesAndGradients)
{
  std::vector<int> vertexOffsets{4, 8, 8, 10};

  const int meshDimensions = 3;
  const int dataDimensions = 1;
  const int nValues        = 4;
  const int nTimeSteps     = 3;

  mesh::PtrMesh dummyMesh(new mesh::Mesh("DummyMesh", meshDimensions, testing::nextMeshID()));
  dummyMesh->setVertexOffsets(vertexOffsets);

  mesh::PtrData fromData(new mesh::Data("from", -1, dataDimensions));
  fromData->requireDataGradient();

  Eigen::VectorXd insert0(nValues);
  insert0 << 1.0, 2.0, 3.0, 4.0;
  Eigen::VectorXd insert05(nValues);
  insert05 << 10.0, 20.0, 30.0, 40.0;
  Eigen::VectorXd insert1(nValues);
  insert1 << 100.0, 200.0, 300.0, 400.0;
  Eigen::MatrixXd insertGradients0(nValues, meshDimensions);
  insertGradients0 << 1.0, 2.0, 3.0,
      4.0, 1.0, 2.0,
      3.0, 4.0, 1.0,
      2.0, 3.0, 4.0;
  Eigen::MatrixXd insertGradients05(nValues, meshDimensions);
  insertGradients05 << 11.0, 21.0, 31.0,
      41.0, 11.0, 21.0,
      31.0, 41.0, 11.0,
      21.0, 31.0, 41.0;
  Eigen::MatrixXd insertGradients1(nValues, meshDimensions);
  insertGradients1 << 111.0, 211.0, 311.0,
      411.0, 111.0, 211.0,
      311.0, 411.0, 111.0,
      211.0, 311.0, 411.0;

  cplscheme::PtrCouplingData fromDataPtr = makeCouplingData(fromData, dummyMesh, true);

  fromDataPtr->setSampleAtTime(0, time::Sample{dataDimensions, insert0, insertGradients0});
  fromDataPtr->setSampleAtTime(0.5, time::Sample{dataDimensions, insert05, insertGradients05});
  fromDataPtr->setSampleAtTime(1, time::Sample{dataDimensions, insert1, insertGradients1});

  const auto serialized = serialize::SerializedStamples::serialize(fromDataPtr);

  Eigen::VectorXd expectedSerializedValues(nTimeSteps * nValues);
  expectedSerializedValues << 1.0, 10.0, 100.0, 2.0, 20.0, 200.0, 3.0, 30.0, 300.0, 4.0, 40.0, 400.0;

  Eigen::VectorXd expectedSerializedGradients(nTimeSteps * nValues * meshDimensions);
  expectedSerializedGradients << 1.0, 11.0, 111.0, 4.0, 41.0, 411.0, 3.0, 31.0, 311.0, 2.0, 21.0, 211.0, 2.0, 21.0, 211.0, 1.0, 11.0, 111.0, 4.0, 41.0, 411.0, 3.0, 31.0, 311.0, 3.0, 31.0, 311.0, 2.0, 21.0, 211.0, 1.0, 11.0, 111.0, 4.0, 41.0, 411.0;

  for (int i = 0; i < nTimeSteps * nValues; i++) {
    BOOST_TEST(testing::equals(serialized.values()(i), expectedSerializedValues(i)));
  }

  for (int i = 0; i < nTimeSteps * meshDimensions * nValues; i++) {
    BOOST_TEST(testing::equals(serialized.gradients()(i), expectedSerializedGradients(i)));
  }
}

BOOST_AUTO_TEST_CASE(DeserializeValuesAndGradients)
{
  std::vector<int> vertexOffsets{4, 8, 8, 10};

  const int meshDimensions = 3;
  const int nValues        = 4;
  const int nTimeSteps     = 3;

  Eigen::VectorXd serializedValues(nTimeSteps * nValues);
  serializedValues << 1.0, 10.0, 100.0, 2.0, 20.0, 200.0, 3.0, 30.0, 300.0, 4.0, 40.0, 400.0;

  Eigen::VectorXd serializedGradients(nTimeSteps * nValues * meshDimensions);
  serializedGradients << 1.0, 11.0, 111.0, 4.0, 41.0, 411.0, 3.0, 31.0, 311.0, 2.0, 21.0, 211.0, 2.0, 21.0, 211.0, 1.0, 11.0, 111.0, 4.0, 41.0, 411.0, 3.0, 31.0, 311.0, 3.0, 31.0, 311.0, 2.0, 21.0, 211.0, 1.0, 11.0, 111.0, 4.0, 41.0, 411.0;

  mesh::PtrMesh dummyMesh(new mesh::Mesh("DummyMesh", 3, testing::nextMeshID()));
  dummyMesh->setVertexOffsets(vertexOffsets);

  mesh::PtrData toData(new mesh::Data("to", -1, 1));
  toData->requireDataGradient();

  cplscheme::PtrCouplingData toDataPtr = makeCouplingData(toData, dummyMesh, true);
  toDataPtr->sample().values           = Eigen::VectorXd(nValues);
  toDataPtr->sample().gradients        = Eigen::MatrixXd(nValues, meshDimensions);
  toDataPtr->sample().setZero();

  toDataPtr->setSampleAtTime(0, toDataPtr->sample());

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

  Eigen::MatrixXd insertGradients0(nValues, meshDimensions);
  insertGradients0 << 1.0, 2.0, 3.0,
      4.0, 1.0, 2.0,
      3.0, 4.0, 1.0,
      2.0, 3.0, 4.0;
  Eigen::MatrixXd insertGradients05(nValues, meshDimensions);
  insertGradients05 << 11.0, 21.0, 31.0,
      41.0, 11.0, 21.0,
      31.0, 41.0, 11.0,
      21.0, 31.0, 41.0;
  Eigen::MatrixXd insertGradients1(nValues, meshDimensions);
  insertGradients1 << 111.0, 211.0, 311.0,
      411.0, 111.0, 211.0,
      311.0, 411.0, 111.0,
      211.0, 311.0, 411.0;

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

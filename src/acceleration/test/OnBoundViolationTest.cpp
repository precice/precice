#ifndef PRECICE_NO_MPI

#include <Eigen/Core>
#include <Eigen/Dense>

#include "acceleration/Acceleration.hpp"
#include "acceleration/IQNILSAcceleration.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(AccelerationTests)

using namespace precice;
using namespace precice::acceleration;
using precice::testing::makeCouplingData;

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(CheckBoundViolation)
{
  PRECICE_TEST();

  const int scalarDataDimension = 1;
  const int vectorDataDimension = 3;

  auto dummyMesh = testing::makeDummy2DMesh(2);

  // only set lower bound for scalar data
  std::vector<std::optional<double>> lowerBoundedL = {0};
  std::vector<std::optional<double>> lowerBoundedU(scalarDataDimension);
  // only set upper bound for scalar data
  std::vector<std::optional<double>> upperBoundedL(scalarDataDimension);
  std::vector<std::optional<double>> upperBoundedU = {0};
  // set both lower and upper bound for scalar data
  std::vector<std::optional<double>> lowerUpperBoundedL = {0};
  std::vector<std::optional<double>> lowerUpperBoundedU = {1};
  // set upper bound for z component in vector data
  std::vector<std::optional<double>> zUpperBoundedL(vectorDataDimension);
  std::vector<std::optional<double>> zUpperBoundedU(vectorDataDimension);
  zUpperBoundedU[2] = 1;
  // set lower bound for x-component
  std::vector<std::optional<double>> xLowerBoundedL(vectorDataDimension);
  std::vector<std::optional<double>> xLowerBoundedU(vectorDataDimension);
  xLowerBoundedL[0] = 0;
  // set upper bound for y-component
  std::vector<std::optional<double>> yUpperBoundedL(vectorDataDimension);
  std::vector<std::optional<double>> yUpperBoundedU(vectorDataDimension);
  yUpperBoundedU[1] = 1;
  // set upper bound for x-component
  std::vector<std::optional<double>> xUpperBoundedL(vectorDataDimension);
  std::vector<std::optional<double>> xUpperBoundedU(vectorDataDimension);
  xUpperBoundedU[0] = 1;

  mesh::PtrData lowerBounded(new mesh::Data("lowerBounded", 0, scalarDataDimension, 1, 2, lowerBoundedL, lowerBoundedU));
  lowerBounded->emplaceSampleAtTime(0.0, {0., 0.});
  mesh::PtrData upperBounded(new mesh::Data("upperBounded", 1, scalarDataDimension, 1, 2, upperBoundedL, upperBoundedU));
  upperBounded->emplaceSampleAtTime(0.0, {0., 0.});
  mesh::PtrData lowerUpperBounded(new mesh::Data("lowerUpperBounded", 2, scalarDataDimension, 1, 2, lowerUpperBoundedL, lowerUpperBoundedU));
  lowerUpperBounded->emplaceSampleAtTime(0.0, {0., 0.});

  mesh::PtrData zUpperBounded(new mesh::Data("zUpperBounded", 3, vectorDataDimension, 1, 2, zUpperBoundedL, zUpperBoundedU));
  zUpperBounded->emplaceSampleAtTime(0.0, {0., 0., 0., 0., 0., 0.});
  mesh::PtrData xLowerBounded(new mesh::Data("xLowerBounded", 4, vectorDataDimension, 1, 2, xLowerBoundedL, xLowerBoundedU));
  xLowerBounded->emplaceSampleAtTime(0.0, {0., 0., 0., 0., 0., 0.});
  mesh::PtrData yUpperBounded(new mesh::Data("yUpperBounded", 5, vectorDataDimension, 1, 2, yUpperBoundedL, yUpperBoundedU));
  yUpperBounded->emplaceSampleAtTime(0.0, {0., 0., 0., 0., 0., 0.});
  mesh::PtrData xUpperBounded(new mesh::Data("xUpperBounded", 6, vectorDataDimension, 1, 2, xUpperBoundedL, xUpperBoundedU));
  xUpperBounded->emplaceSampleAtTime(0.0, {0., 0., 0., 0., 0., 0.});

  cplscheme::PtrCouplingData lowerBoundedPtr      = makeCouplingData(lowerBounded, dummyMesh, true);
  cplscheme::PtrCouplingData upperBoundedPtr      = makeCouplingData(upperBounded, dummyMesh, true);
  cplscheme::PtrCouplingData lowerUpperBoundedPtr = makeCouplingData(lowerUpperBounded, dummyMesh, true);
  cplscheme::PtrCouplingData zUpperBoundedPtr     = makeCouplingData(zUpperBounded, dummyMesh, true);
  cplscheme::PtrCouplingData xLowerBoundedPtr     = makeCouplingData(xLowerBounded, dummyMesh, true);
  cplscheme::PtrCouplingData yUpperBoundedPtr     = makeCouplingData(yUpperBounded, dummyMesh, true);
  cplscheme::PtrCouplingData xUpperBoundedPtr     = makeCouplingData(xUpperBounded, dummyMesh, true);
  for (auto &cd : {lowerBoundedPtr, upperBoundedPtr, lowerUpperBoundedPtr, zUpperBoundedPtr, xLowerBoundedPtr, yUpperBoundedPtr, xUpperBoundedPtr}) {
    cd->storeIteration();
  }

  DataMap cplData{
      {0, lowerBoundedPtr},
      {1, upperBoundedPtr},
      {2, lowerUpperBoundedPtr},
      {3, zUpperBoundedPtr},
      {4, xLowerBoundedPtr},
      {5, yUpperBoundedPtr},
      {6, xUpperBoundedPtr}};

  IQNILSAcceleration acceleration(1.0, false, 10, 0, 0, 1e-10, {0}, acceleration::Acceleration::OnBoundViolation::Clamp, nullptr, false);
  acceleration.initialize(cplData);
  Eigen::VectorXd testData(3 * 2 + 4 * 6);
  testData << 2.0, -0.1, 2.0, 0.1, -0.1, 1.1, 0.1, 0.1, 1.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, -0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 1.1, 0.1, 1.1, 0.1, 0.1, 0.1, 0.1, 0.1;
  auto violatingIDs = acceleration.checkBoundViolation(testData, cplData);

  for (int i = 0; i < violatingIDs.size(); i++) {
  }
  std::vector<int> expectedViolations{0, 1, 2, 3, 4, 5, 6};
  BOOST_TEST(violatingIDs == expectedViolations, boost::test_tools::per_element());
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(testOnBoundViolationClamp)
{
  PRECICE_TEST();

  Eigen::VectorXd _data(6);
  _data << 0.1, 0.3, 0.4, -0.4, 1.2, 0.5;
  Eigen::Map<Eigen::MatrixXd> dataMap(_data.data(), 3, 2);

  Eigen::VectorXd _dataAfterUpdate(6);
  _dataAfterUpdate << 0.1, 0.3, 0.4, 0.0, 1.0, 0.5;
  Eigen::Map<Eigen::MatrixXd> dataAfterUpdateMap(_dataAfterUpdate.data(), 3, 2);

  std::vector<std::optional<double>> _lowerBound(2);
  _lowerBound[0] = 0;
  _lowerBound[1] = 0;
  std::vector<std::optional<double>> _upperBound(2);
  _upperBound[0] = 1;
  _upperBound[1] = 1;

  IQNILSAcceleration acceleration(1.0, false, 10, 0, 0, 1e-10, {0}, acceleration::Acceleration::OnBoundViolation::Clamp, nullptr, false);
  acceleration.clampToBounds(dataMap, _lowerBound, _upperBound);

  // check if the over small single values have been correctly truncated
  BOOST_TEST(testing::equals(dataMap, dataAfterUpdateMap, 1e-10));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(testOnBoundViolationScale)
{
  PRECICE_TEST();

  Eigen::VectorXd _data(6);
  _data << 0.1, 0.3, 0.4, -0.4, 1.3, 0.5;
  Eigen::Map<Eigen::MatrixXd> dataMap(_data.data(), 3, 2);

  Eigen::VectorXd _dataUpdate(6);
  _dataUpdate << -0.6, -0.3, 0.3, -0.6, 0.9, 0.12;
  Eigen::Map<Eigen::MatrixXd> dataUpdateMap(_dataUpdate.data(), 3, 2);

  Eigen::VectorXd _dataAfterUpdate(6);
  _dataAfterUpdate << 0.5, 0.5, 0.2, 0.0, 0.7, 0.42;
  Eigen::Map<Eigen::MatrixXd> dataAfterUpdateMap(_dataAfterUpdate.data(), 3, 2);

  std::vector<std::optional<double>> _lowerBound(2);
  _lowerBound[0] = 0;
  _lowerBound[1] = 0;
  std::vector<std::optional<double>> _upperBound(2);
  _upperBound[0] = 1;
  _upperBound[1] = 1;

  IQNILSAcceleration acceleration(1.0, false, 10, 0, 0, 1e-10, {0}, acceleration::Acceleration::OnBoundViolation::ScaleToBound, nullptr, false);
  auto               scaleFactor = acceleration.scaleToBounds(dataMap, dataUpdateMap, _lowerBound, _upperBound);
  dataMap -= dataUpdateMap * scaleFactor;

  // check if the over small single values have been correctly truncated
  BOOST_TEST(testing::equals(dataMap, dataAfterUpdateMap, 1e-10));
}

BOOST_AUTO_TEST_SUITE_END()

#endif

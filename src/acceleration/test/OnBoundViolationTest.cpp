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

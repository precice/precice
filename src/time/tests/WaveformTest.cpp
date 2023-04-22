#include <Eigen/Core>
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "testing/WaveformFixture.hpp"
#include "time/Time.hpp"
#include "time/Waveform.hpp"

using namespace precice;
using namespace precice::time;

BOOST_AUTO_TEST_SUITE(TimeTests)
BOOST_AUTO_TEST_SUITE(WaveformTests)

BOOST_AUTO_TEST_CASE(testInitialization)
{
  PRECICE_TEST(1_rank);
  const int interpolationOrder = 0;
  Waveform  waveform(interpolationOrder);
  const int valuesSize = 1;
  waveform.initialize(Eigen::VectorXd::Zero(valuesSize));

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 0.0));
}

BOOST_AUTO_TEST_CASE(testInitializationVector)
{
  PRECICE_TEST(1_rank);

  const int interpolationOrder = 1;
  Waveform  waveform(interpolationOrder);
  const int valuesSize = 3;
  waveform.initialize(Eigen::VectorXd::Zero(valuesSize));

  for (int i = 0; i < valuesSize; i++) {
    BOOST_TEST(testing::equals(waveform.sample(0.0)(i), 0.0));
    BOOST_TEST(testing::equals(waveform.sample(1.0)(i), 0.0));
  }
}

BOOST_AUTO_TEST_SUITE(InterpolationTests)

BOOST_AUTO_TEST_CASE(testInterpolateDataZerothOrder)
{
  PRECICE_TEST(1_rank);

  testing::WaveformFixture fixture;

  // Test zeroth order interpolation
  const int interpolationOrder = 0;
  Waveform  waveform(interpolationOrder);
  const int valuesSize = 1;
  waveform.initialize(Eigen::VectorXd::Zero(valuesSize));

  BOOST_TEST(fixture.valuesSize(waveform) == valuesSize);
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);

  Eigen::VectorXd value(1);
  value(0) = 1.0;
  waveform.store(value);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 1.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 1.0));

  value(0) = 2.0;
  waveform.store(value);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 2.0));

  waveform.moveToNextWindow();
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 2.0));

  value(0) = 3.0;
  waveform.store(value);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 3.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 3.0));
}

BOOST_AUTO_TEST_CASE(testInterpolateDataFirstOrder)
{
  PRECICE_TEST(1_rank);

  testing::WaveformFixture fixture;

  // Test first order interpolation
  const int interpolationOrder = 1;
  Waveform  waveform(interpolationOrder);
  const int valuesSize = 1;
  waveform.initialize(Eigen::VectorXd::Zero(valuesSize));

  BOOST_TEST(fixture.valuesSize(waveform) == valuesSize);
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);

  Eigen::VectorXd value(1);
  value(0) = 1.0;
  waveform.store(value);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 0.5));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 1.0));

  value(0) = 2.0;
  waveform.store(value);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 1.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 2.0));

  waveform.moveToNextWindow();
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 2.0));

  value(0) = 3.0;
  waveform.store(value);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 2.5));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 3.0));
}

// Remove or modify this feature? Creating a second degree interpolant by using data from previous windows is difficult, because this would require several pieces of data during initialization. What would be useful: Generating a second order interpolant from multiple samples in a single window (if available). This would go into the least-squares direction
BOOST_AUTO_TEST_CASE(testInterpolateDataSecondOrder)
{
  PRECICE_TEST(1_rank);

  testing::WaveformFixture fixture;

  // Test second order interpolation, but there are not enough samples. Therefore, always only first order.
  const int interpolationOrder = 2;
  Waveform  waveform(interpolationOrder);
  const int valuesSize = 1;
  waveform.initialize(Eigen::VectorXd::Zero(valuesSize));

  BOOST_TEST(fixture.valuesSize(waveform) == valuesSize);
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);

  Eigen::VectorXd value(1);
  value(0) = 1.0;
  waveform.store(value);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 0.5));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 1.0));

  value(0) = 2.0;
  waveform.store(value);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 1.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 2.0));

  waveform.moveToNextWindow();
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 2.0));

  value(0) = 8.0;
  waveform.store(value);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 5.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 8.0));

  value(0) = 4.0;
  waveform.store(value);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 3.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 4.0));

  waveform.moveToNextWindow();
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 4.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 4.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 4.0));

  value(0) = 8.0;
  waveform.store(value);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 4.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 6.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 8.0));
}

BOOST_AUTO_TEST_CASE(testInterpolateDataFirstOrderVector)
{
  PRECICE_TEST(1_rank);

  testing::WaveformFixture fixture;

  // Test first order interpolation
  const int interpolationOrder = 1;
  Waveform  waveform(interpolationOrder);
  const int valuesSize = 3;
  waveform.initialize(Eigen::VectorXd::Zero(valuesSize));

  BOOST_TEST(fixture.valuesSize(waveform) == valuesSize);
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);

  Eigen::VectorXd value(valuesSize);
  value << 1, 2, 3;
  waveform.store(value);

  for (int i = 0; i < valuesSize; i++) {
    BOOST_TEST(testing::equals(waveform.sample(0.0)(i), 0 * value[i]));
    BOOST_TEST(testing::equals(waveform.sample(0.5)(i), 0.5 * value[i]));
    BOOST_TEST(testing::equals(waveform.sample(1.0)(i), value[i]));
  }

  value << 2, 4, 2;
  waveform.store(value);

  for (int i = 0; i < valuesSize; i++) {
    BOOST_TEST(testing::equals(waveform.sample(0.0)(i), 0 * value[i]));
    BOOST_TEST(testing::equals(waveform.sample(0.5)(i), 0.5 * value[i]));
    BOOST_TEST(testing::equals(waveform.sample(1.0)(i), value[i]));
  }

  waveform.moveToNextWindow();
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);

  for (int i = 0; i < valuesSize; i++) {
    BOOST_TEST(testing::equals(waveform.sample(0.0)(i), value[i]));
    BOOST_TEST(testing::equals(waveform.sample(0.5)(i), value[i]));
    BOOST_TEST(testing::equals(waveform.sample(1.0)(i), value[i]));
  }

  Eigen::VectorXd value0 = value;
  value << 1, 2, 3;
  waveform.store(value);

  for (int i = 0; i < valuesSize; i++) {
    BOOST_TEST(testing::equals(waveform.sample(0.0)(i), value0[i]));
    BOOST_TEST(testing::equals(waveform.sample(0.5)(i), 0.5 * value0[i] + 0.5 * value[i]));
  }
}

BOOST_AUTO_TEST_CASE(testPiecewiseInterpolateDataZerothOrder)
{
  PRECICE_TEST(1_rank);

  testing::WaveformFixture fixture;

  // Test zeroth order interpolation
  const int interpolationOrder = 0;
  Waveform  waveform(interpolationOrder);
  const int valuesSize = 1;
  waveform.initialize(Eigen::VectorXd::Zero(valuesSize));

  BOOST_TEST(fixture.valuesSize(waveform) == valuesSize);
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);
  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 0.0));

  Eigen::VectorXd value(1);
  value(0) = 0.5;
  waveform.store(value, 0.5);

  value(0) = 1.0;
  waveform.store(value, 1.0);

  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 0.5));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 0.5));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 1.0));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 1.0));

  value(0) = 1.5;
  waveform.store(value, 0.5);

  value(0) = 2.0;
  waveform.store(value, 1.0);

  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 1.5));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 1.5));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 2.0));

  waveform.moveToNextWindow();
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);

  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 2.0));

  value(0) = 3.0;
  waveform.store(value, 1.0);

  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 3.0));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 3.0));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 3.0));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 3.0));

  value(0) = 1.5;
  waveform.store(value, 0.5);

  value(0) = 4.0;
  waveform.store(value, 1.0);

  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 1.5));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 1.5));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 4.0));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 4.0));
}

BOOST_AUTO_TEST_CASE(testPiecewiseInterpolateDataFirstOrder)
{
  PRECICE_TEST(1_rank);

  testing::WaveformFixture fixture;

  // Test zeroth order interpolation
  const int interpolationOrder = 1;
  Waveform  waveform(interpolationOrder);
  const int valuesSize = 1;
  waveform.initialize(Eigen::VectorXd::Zero(valuesSize));

  BOOST_TEST(fixture.valuesSize(waveform) == valuesSize);
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);
  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 0.00));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 0.00));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 0.00));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 0.00));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 0.00));

  Eigen::VectorXd value(1);
  value(0) = 0.5;
  waveform.store(value, 0.5);

  value(0) = 1.0;
  waveform.store(value, 1.0);

  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 0.00));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 0.25));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 0.50));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 0.75));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 1.00));

  value(0) = 1.5;
  waveform.store(value, 0.5);

  value(0) = 2.0;
  waveform.store(value, 1.0);

  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 0.00));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 0.75));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 1.50));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 1.75));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 2.00));

  waveform.moveToNextWindow();
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);

  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 2.00));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 2.00));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 2.00));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 2.00));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 2.00));

  value(0) = 3.0;
  waveform.store(value, 1.0);

  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 2.00));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 2.25));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 2.50));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 2.75));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 3.00));

  value(0) = 1.5;
  waveform.store(value, 0.5);

  value(0) = 4.0;
  waveform.store(value, 1.0);

  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 2.00));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 1.75));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 1.50));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 2.75));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 4.00));
}

BOOST_AUTO_TEST_CASE(testPiecewiseInterpolateDataSecondOrder)
{
  PRECICE_TEST(1_rank);

  testing::WaveformFixture fixture;

  // Test zeroth order interpolation
  const int interpolationOrder = 2;
  Waveform  waveform(interpolationOrder);
  const int valuesSize = 1;

  waveform.initialize(Eigen::VectorXd::Zero(valuesSize));

  BOOST_TEST(fixture.valuesSize(waveform) == valuesSize);
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);
  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 0.00));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 0.00));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 0.00));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 0.00));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 0.00));

  Eigen::VectorXd value(1);
  value(0) = 0;
  waveform.store(value, 0.5);

  value(0) = 2;
  waveform.store(value, 1.0);

  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 3);
  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 0.00));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), -0.25));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 0.00));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 0.75));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 2.00));
}

BOOST_AUTO_TEST_CASE(testPiecewiseInterpolateDataThirdOrder)
{
  PRECICE_TEST(1_rank);

  testing::WaveformFixture fixture;

  // Test zeroth order interpolation
  const int interpolationOrder = 3;
  Waveform  waveform(interpolationOrder);
  const int valuesSize = 1;

  waveform.initialize(Eigen::VectorXd::Zero(valuesSize));

  BOOST_TEST(fixture.valuesSize(waveform) == valuesSize);
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);
  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 0.00));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 0.00));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 0.00));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 0.00));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 0.00));

  // linearly increasing values
  Eigen::VectorXd value(1);
  for (double t : std::vector<double>{0.25, 0.5, 0.75, 1}) {
    value(0) = t;
    waveform.store(value, t);
  }

  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 5);

  for (double t : std::vector<double>{0.1, 0.2, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t));
  }

  // quadratically increasing values
  for (double t : std::vector<double>{0.25, 0.5, 0.75, 1}) {
    value(0) = t * t;
    waveform.store(value, t);
  }

  // interpolates given values
  for (double t : std::vector<double>{0.25, 0.5, 0.75, 1}) {
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t * t));
  }

  // introduces no approximation error w.r.t function
  for (double t : std::vector<double>{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t * t));
  }

  // cubically increasing values
  for (double t : std::vector<double>{0.25, 0.5, 0.75, 1}) {
    value(0) = t * t * t;
    waveform.store(value, t);
  }

  // interpolates given values
  for (double t : std::vector<double>{0.25, 0.5, 0.75, 1}) {
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t * t * t));
  }

  // introduces no approximation error w.r.t function
  for (double t : std::vector<double>{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t * t * t));
  }

  // cubically increasing values, but with non-uniform spacing
  for (double t : std::vector<double>{0.01, 0.1, 0.2, 1}) {
    value(0) = t * t * t;
    waveform.store(value, t);
  }

  // interpolates given values
  for (double t : std::vector<double>{0.01, 0.1, 0.2, 1}) {
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t * t * t));
  }

  // introduces no approximation error w.r.t function
  for (double t : std::vector<double>{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t * t * t));
  }

  // quadratically increasing values, but with non-uniform spacing
  for (double t : std::vector<double>{0.25, 0.5, 0.75, 1}) {
    value(0) = t * t * t * t;
    waveform.store(value, t);
  }

  // interpolates given values
  for (double t : std::vector<double>{0.25, 0.5, 0.75, 1}) {
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t * t * t * t));
  }

  // introduces approximation error w.r.t function
  for (double t : std::vector<double>{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {
    double tol = 0.015625; // error < h**3 = 0.015625
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t * t * t * t, tol));
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

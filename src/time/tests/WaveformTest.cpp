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

  BOOST_TEST(fixture.maxNumberOfStoredSamples(waveform) == 1);
  BOOST_TEST(fixture.valuesSize(waveform) == valuesSize);
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 1);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 0.0));

  Eigen::VectorXd value(1);
  value(0) = 1.0;
  waveform.store(value);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 1.0));

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0, 0), 1.0));

  value(0) = 2.0;
  waveform.store(value);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0, 0), 2.0));

  waveform.moveToNextWindow();
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 1);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 2.0));

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0, 0), 2.0));

  value(0) = 3.0;
  waveform.store(value);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0, 0), 3.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0, 0), 3.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0, 0), 3.0));
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

  BOOST_TEST(fixture.maxNumberOfStoredSamples(waveform) == 2);
  BOOST_TEST(fixture.valuesSize(waveform) == valuesSize);
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 1);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 0.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 0.0));

  Eigen::VectorXd value(1);
  value(0) = 1.0;
  waveform.store(value);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 1.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 0.0));

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0, 0), 1.0));

  value(0) = 2.0;
  waveform.store(value);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0, 0), 2.0));

  waveform.moveToNextWindow();
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 2.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 2.0));

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0, 0), 2.0));

  value(0) = 3.0;
  waveform.store(value);

  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 3.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 2.0));

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0, 0), 2.5));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0, 0), 3.0));
}

BOOST_AUTO_TEST_CASE(testInterpolateDataSecondOrder)
{
  PRECICE_TEST(1_rank);

  testing::WaveformFixture fixture;

  // Test second order interpolation
  const int interpolationOrder = 2;
  Waveform  waveform(interpolationOrder);
  const int valuesSize = 1;
  waveform.initialize(Eigen::VectorXd::Zero(valuesSize));

  BOOST_TEST(fixture.maxNumberOfStoredSamples(waveform) == 3);
  BOOST_TEST(fixture.valuesSize(waveform) == valuesSize);
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 1);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 0.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 0.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 2), 0.0));

  Eigen::VectorXd value(1);
  value(0) = 1.0;
  waveform.store(value);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 1.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 0.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 2), 0.0));

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0, 0), 1.0));

  value(0) = 2.0;
  waveform.store(value);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 2.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 0.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 2), 0.0));

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0, 0), 2.0));

  waveform.moveToNextWindow();
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 2.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 2.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 2), 0.0));

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0, 0), 2.0));

  value(0) = 8.0;
  waveform.store(value);

  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 8.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 2.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 2), 0.0));

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0, 0), 5.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0, 0), 8.0));

  value(0) = 4.0;
  waveform.store(value);

  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 4.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 2.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 2), 0.0));

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0, 0), 3.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0, 0), 4.0));

  waveform.moveToNextWindow();
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 3);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 4.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 4.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 2), 2.0));

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0, 0), 4.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0, 0), 4.25));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0, 0), 4.0));

  value(0) = 8.0;
  waveform.store(value);

  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 8.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 4.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 2), 2.0));

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0, 0), 4.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0, 0), 5.75));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0, 0), 8.0));
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

  BOOST_TEST(fixture.maxNumberOfStoredSamples(waveform) == 2);
  BOOST_TEST(fixture.valuesSize(waveform) == valuesSize);
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 1);
  for(int i : boost::irange(valuesSize)){
    BOOST_TEST(testing::equals(fixture.getValue(waveform, i, 0), 0.0));
    BOOST_TEST(testing::equals(fixture.getValue(waveform, i, 1), 0.0));
  }

  Eigen::VectorXd value(valuesSize);
  value << 1,2,3;
  waveform.store(value);
  for(int i : boost::irange(valuesSize)){
    BOOST_TEST(testing::equals(fixture.getValue(waveform, i, 0), value[i]));
    BOOST_TEST(testing::equals(fixture.getValue(waveform, i, 1), 0.0));
  }

  for(int i : boost::irange(valuesSize)){
    BOOST_TEST(testing::equals(waveform.sample(0.0)(i, 0), value[i]));
    BOOST_TEST(testing::equals(waveform.sample(0.5)(i, 0), value[i]));
    BOOST_TEST(testing::equals(waveform.sample(1.0)(i, 0), value[i]));
  }

  value << 2,4,2;
  waveform.store(value);
  for(int i : boost::irange(valuesSize)){
    BOOST_TEST(testing::equals(fixture.getValue(waveform, i, 0), value[i]));
    BOOST_TEST(testing::equals(fixture.getValue(waveform, i, 1), 0.0));
  }

  for(int i : boost::irange(valuesSize)){
    BOOST_TEST(testing::equals(waveform.sample(0.0)(i, 0), value[i]));
    BOOST_TEST(testing::equals(waveform.sample(0.5)(i, 0), value[i]));
    BOOST_TEST(testing::equals(waveform.sample(1.0)(i, 0), value[i]));
  }

  waveform.moveToNextWindow();
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);
  for(int i : boost::irange(valuesSize)){
    BOOST_TEST(testing::equals(fixture.getValue(waveform, i, 0), value[i]));
    BOOST_TEST(testing::equals(fixture.getValue(waveform, i, 1), value[i]));
  }

  for(int i : boost::irange(valuesSize)){
    BOOST_TEST(testing::equals(waveform.sample(0.0)(i, 0), value[i]));
    BOOST_TEST(testing::equals(waveform.sample(0.5)(i, 0), value[i]));
    BOOST_TEST(testing::equals(waveform.sample(1.0)(i, 0), value[i]));
  }

  Eigen::VectorXd value0 = value;
  value << 1,2,3;
  waveform.store(value);

  for(int i : boost::irange(valuesSize)){
    BOOST_TEST(testing::equals(fixture.getValue(waveform, i, 0), value[i]));
    BOOST_TEST(testing::equals(fixture.getValue(waveform, i, 1), value0[i]));
  }

  for(int i : boost::irange(valuesSize)){
    BOOST_TEST(testing::equals(waveform.sample(0.0)(i, 0), value0[i]));
    BOOST_TEST(testing::equals(waveform.sample(0.5)(i, 0), 0.5*value0[i] + 0.5*value[i]));
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

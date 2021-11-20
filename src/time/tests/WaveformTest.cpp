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

BOOST_AUTO_TEST_SUITE(ExtrapolationTests)
BOOST_AUTO_TEST_CASE(testExtrapolateDataFirstOrder)
{
  PRECICE_TEST(1_rank);

  testing::WaveformFixture fixture;

  // Test first order extrapolation
  const int extrapolationOrder = 1;
  Waveform  waveform(extrapolationOrder, Time::UNDEFINED_INTERPOLATION_ORDER);
  const int valuesSize = 1;
  waveform.initialize(valuesSize);
  BOOST_TEST(fixture.sizeOfSampleStorage(waveform) == 2);
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 1);
  BOOST_TEST(fixture.valuesSize(waveform) == 1);

  // use zero initial data
  waveform.moveToNextWindow();
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 0.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 0.0));

  Eigen::VectorXd value(1);
  value(0) = 1.0;
  waveform.store(value);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 1.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 0.0));
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);
  waveform.moveToNextWindow(); // applies first order extrapolation in second window
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 2.0)); // = 2*1 - 0
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 1.0));

  value(0) = 4.0;
  waveform.store(value);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 4.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 1.0));
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);
  waveform.moveToNextWindow(); // applies first order extrapolation in third window
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 7.0)); // = 2*4 - 1
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 4.0));

  value(0) = 8.0;
  waveform.store(value);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 8.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 4.0));
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);
  waveform.moveToNextWindow(); // applies first order extrapolation in forth window
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 12.0)); // 10.0 = 2 * 8 - 4
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 8.0));
}

BOOST_AUTO_TEST_CASE(testExtrapolateDataSecondOrder)
{
  PRECICE_TEST(1_rank);

  testing::WaveformFixture fixture;

  // Test second order extrapolation
  const int extrapolationOrder = 2;
  Waveform  waveform(extrapolationOrder, Time::UNDEFINED_INTERPOLATION_ORDER);
  const int valuesSize = 1;
  waveform.initialize(valuesSize);
  BOOST_TEST(fixture.sizeOfSampleStorage(waveform) == 3);
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 1);
  BOOST_TEST(fixture.valuesSize(waveform) == 1);

  // use zero initial data
  waveform.moveToNextWindow();
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 0.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 0.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 2), 0.0));

  Eigen::VectorXd value(1);
  value(0) = 1.0;
  waveform.store(value);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 1.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 0.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 2), 0.0));
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);
  waveform.moveToNextWindow(); // applies first order extrapolation in second window
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 3);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 2.0)); // = 2*1 - 0
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 1.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 2), 0.0));

  value(0) = 4.0;
  waveform.store(value);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 4.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 1.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 2), 0.0));
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 3);
  waveform.moveToNextWindow(); // applies second order extrapolation in third window
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 3);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 8.0)); // = 2.5*4 - 2 * 1 + 0.5 * 0
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 4.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 2), 1.0));

  value(0) = 8.0;
  waveform.store(value);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 8.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 4.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 2), 1.0));
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 3);
  waveform.moveToNextWindow(); // applies second order extrapolation in fourth window
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 3);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 12.5)); // = 2.5 * 8 - 2 * 4 + 0.5 * 1
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 8.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 2), 4.0));

  value(0) = 16.0;
  waveform.store(value);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 16.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 8.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 2), 4.0));
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 3);
  waveform.moveToNextWindow(); // applies second order extrapolation in fifth window
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 3);
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 0), 26.0)); // = 2.5 * 16.0 - 2 * 8 + 0.5 * 4
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 1), 16.0));
  BOOST_TEST(testing::equals(fixture.getValue(waveform, 0, 2), 8.0));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(InterpolationTests)
BOOST_AUTO_TEST_CASE(testInterpolateDataZerothOrder)
{
  PRECICE_TEST(1_rank);

  testing::WaveformFixture fixture;

  // Test zeroth order interpolation
  const int extrapolationOrder = 0;
  const int interpolationOrder = 0;
  Waveform  waveform(extrapolationOrder, interpolationOrder);
  const int valuesSize = 1;
  waveform.initialize(valuesSize);

  BOOST_TEST(fixture.sizeOfSampleStorage(waveform) == 2);
  BOOST_TEST(fixture.valuesSize(waveform) == 1);
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

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0, 0), 3.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0, 0), 3.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0, 0), 3.0));
}

BOOST_AUTO_TEST_CASE(testInterpolateDataFirstOrder)
{
  PRECICE_TEST(1_rank);

  testing::WaveformFixture fixture;

  // Test first order interpolation
  const int extrapolationOrder = 0;
  const int interpolationOrder = 1;
  Waveform  waveform(extrapolationOrder, interpolationOrder);
  const int valuesSize = 1;
  waveform.initialize(valuesSize);

  BOOST_TEST(fixture.sizeOfSampleStorage(waveform) == 2);
  BOOST_TEST(fixture.valuesSize(waveform) == 1);
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
  const int extrapolationOrder = 0;
  const int interpolationOrder = 2;
  Waveform  waveform(extrapolationOrder, interpolationOrder);
  const int valuesSize = 1;
  waveform.initialize(valuesSize);
  BOOST_TEST(fixture.sizeOfSampleStorage(waveform) == 3);
  BOOST_TEST(fixture.valuesSize(waveform) == 1);
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

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

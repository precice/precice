#include <Eigen/Core>
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "time/Waveform.hpp"

using namespace precice;
using namespace precice::time;

BOOST_AUTO_TEST_SUITE(TimeTests)

BOOST_AUTO_TEST_SUITE(WaveformTests)

BOOST_AUTO_TEST_CASE(testExtrapolateData)
{
  PRECICE_TEST(1_rank);

  // Test first order extrapolation
  int      extrapolationOrder = 1;
  Waveform waveform(1, extrapolationOrder);
  BOOST_TEST(waveform.numberOfSamples() == 2);
  BOOST_TEST(waveform.numberOfValidSamples() == 1);
  BOOST_TEST(waveform.numberOfData() == 1);
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 0), 0.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 1), 0.0));

  Eigen::VectorXd value(1);
  value(0) = 1.0;
  waveform.store(value);
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 1), 0.0));
  waveform.moveToNextWindow(extrapolationOrder);
  BOOST_TEST(waveform.numberOfValidSamples() == 2);
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 1), 1.0));

  value(0) = 4.0;
  waveform.store(value);
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 0), 4.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 1), 1.0));
  waveform.moveToNextWindow(extrapolationOrder);
  BOOST_TEST(waveform.numberOfValidSamples() == 2);
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 0), 7.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 1), 4.0));

  // Test second order extrapolation
  extrapolationOrder = 2;
  Waveform waveform2(1, extrapolationOrder);
  BOOST_TEST(waveform2.numberOfSamples() == 3);
  BOOST_TEST(waveform2.numberOfValidSamples() == 1);
  BOOST_TEST(waveform2.numberOfData() == 1);
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 0.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 0.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 0.0));

  value(0) = 1.0;
  waveform2.store(value);
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 0.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 0.0));
  waveform2.moveToNextWindow(extrapolationOrder); // only applies first order extrapolation
  BOOST_TEST(waveform2.numberOfValidSamples() == 2);
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 1.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 0.0));

  value(0) = 4.0;
  waveform2.store(value);
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 4.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 1.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 0.0));
  waveform2.moveToNextWindow(extrapolationOrder); // applies second order extrapolation
  BOOST_TEST(waveform2.numberOfValidSamples() == 3);
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 8.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 4.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 1.0));
}

BOOST_AUTO_TEST_CASE(testInterpolateData)
{
  PRECICE_TEST(1_rank);

  // Test zeroth order interpolation
  int      extrapolationOrder = 0;
  int      interpolationOrder = 0;
  Waveform waveform0(1);

  BOOST_TEST(waveform0.lastTimeWindows().cols() == 2);
  BOOST_TEST(waveform0.lastTimeWindows().rows() == 1);
  BOOST_TEST(waveform0.numberOfValidSamples() == 1);
  BOOST_TEST(testing::equals(waveform0.lastTimeWindows()(0, 0), 0.0));
  BOOST_TEST(testing::equals(waveform0.lastTimeWindows()(0, 1), 0.0));

  Eigen::VectorXd value(1);
  value(0) = 1.0;
  waveform0.store(value);
  BOOST_TEST(testing::equals(waveform0.lastTimeWindows()(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform0.lastTimeWindows()(0, 1), 0.0));

  BOOST_TEST(testing::equals(waveform0.sample(0.0, interpolationOrder)(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform0.sample(0.5, interpolationOrder)(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform0.sample(1.0, interpolationOrder)(0, 0), 1.0));

  value(0) = 2.0;
  waveform0.store(value);

  BOOST_TEST(testing::equals(waveform0.sample(0.0, interpolationOrder)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform0.sample(0.5, interpolationOrder)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform0.sample(1.0, interpolationOrder)(0, 0), 2.0));

  waveform0.moveToNextWindow();
  BOOST_TEST(waveform0.numberOfValidSamples() == 2);
  BOOST_TEST(testing::equals(waveform0.lastTimeWindows()(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform0.lastTimeWindows()(0, 1), 2.0));

  BOOST_TEST(testing::equals(waveform0.sample(0.0, interpolationOrder)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform0.sample(0.5, interpolationOrder)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform0.sample(1.0, interpolationOrder)(0, 0), 2.0));

  value(0) = 3.0;
  waveform0.store(value);

  BOOST_TEST(testing::equals(waveform0.sample(0.0, interpolationOrder)(0, 0), 3.0));
  BOOST_TEST(testing::equals(waveform0.sample(0.5, interpolationOrder)(0, 0), 3.0));
  BOOST_TEST(testing::equals(waveform0.sample(1.0, interpolationOrder)(0, 0), 3.0));

  // Test first order interpolation
  extrapolationOrder = 0;
  interpolationOrder = 1;
  Waveform waveform1(1, extrapolationOrder, interpolationOrder);
  BOOST_TEST(waveform1.lastTimeWindows().cols() == 2);
  BOOST_TEST(waveform1.lastTimeWindows().rows() == 1);
  BOOST_TEST(waveform1.numberOfValidSamples() == 1);
  BOOST_TEST(testing::equals(waveform1.lastTimeWindows()(0, 0), 0.0));
  BOOST_TEST(testing::equals(waveform1.lastTimeWindows()(0, 1), 0.0));

  value(0) = 1.0;
  waveform1.store(value);
  BOOST_TEST(testing::equals(waveform1.lastTimeWindows()(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform1.lastTimeWindows()(0, 1), 0.0));

  BOOST_TEST(testing::equals(waveform1.sample(0.0, interpolationOrder)(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform1.sample(0.5, interpolationOrder)(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform1.sample(1.0, interpolationOrder)(0, 0), 1.0));

  value(0) = 2.0;
  waveform1.store(value);

  BOOST_TEST(testing::equals(waveform1.sample(0.0, interpolationOrder)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform1.sample(0.5, interpolationOrder)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform1.sample(1.0, interpolationOrder)(0, 0), 2.0));

  waveform1.moveToNextWindow();
  BOOST_TEST(waveform1.numberOfValidSamples() == 2);
  BOOST_TEST(testing::equals(waveform1.lastTimeWindows()(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform1.lastTimeWindows()(0, 1), 2.0));

  BOOST_TEST(testing::equals(waveform1.sample(0.0, interpolationOrder)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform1.sample(0.5, interpolationOrder)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform1.sample(1.0, interpolationOrder)(0, 0), 2.0));

  value(0) = 3.0;
  waveform1.store(value);

  BOOST_TEST(testing::equals(waveform1.lastTimeWindows()(0, 0), 3.0));
  BOOST_TEST(testing::equals(waveform1.lastTimeWindows()(0, 1), 2.0));

  BOOST_TEST(testing::equals(waveform1.sample(0.0, interpolationOrder)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform1.sample(0.5, interpolationOrder)(0, 0), 2.5));
  BOOST_TEST(testing::equals(waveform1.sample(1.0, interpolationOrder)(0, 0), 3.0));

  // Test second order interpolation
  extrapolationOrder = 0;
  interpolationOrder = 2;
  Waveform waveform2(1, extrapolationOrder, interpolationOrder);
  BOOST_TEST(waveform2.lastTimeWindows().cols() == 3);
  BOOST_TEST(waveform2.lastTimeWindows().rows() == 1);
  BOOST_TEST(waveform2.numberOfValidSamples() == 1);
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 0.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 0.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 0.0));

  value(0) = 1.0;
  waveform2.store(value);
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 0.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 0.0));

  BOOST_TEST(testing::equals(waveform2.sample(0.0, interpolationOrder)(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform2.sample(0.5, interpolationOrder)(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform2.sample(1.0, interpolationOrder)(0, 0), 1.0));

  value(0) = 2.0;
  waveform2.store(value);
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 0.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 0.0));

  BOOST_TEST(testing::equals(waveform2.sample(0.0, interpolationOrder)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform2.sample(0.5, interpolationOrder)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform2.sample(1.0, interpolationOrder)(0, 0), 2.0));

  waveform2.moveToNextWindow();
  BOOST_TEST(waveform2.numberOfValidSamples() == 2);
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 2.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 0.0));

  BOOST_TEST(testing::equals(waveform2.sample(0.0, interpolationOrder)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform2.sample(0.5, interpolationOrder)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform2.sample(1.0, interpolationOrder)(0, 0), 2.0));

  value(0) = 8.0;
  waveform2.store(value);

  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 8.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 2.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 0.0));

  BOOST_TEST(testing::equals(waveform2.sample(0.0, interpolationOrder)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform2.sample(0.5, interpolationOrder)(0, 0), 5.0));
  BOOST_TEST(testing::equals(waveform2.sample(1.0, interpolationOrder)(0, 0), 8.0));

  value(0) = 4.0;
  waveform2.store(value);

  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 4.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 2.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 0.0));

  BOOST_TEST(testing::equals(waveform2.sample(0.0, interpolationOrder)(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform2.sample(0.5, interpolationOrder)(0, 0), 3.0));
  BOOST_TEST(testing::equals(waveform2.sample(1.0, interpolationOrder)(0, 0), 4.0));

  waveform2.moveToNextWindow();
  BOOST_TEST(waveform2.numberOfValidSamples() == 3);
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 4.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 4.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 2.0));

  BOOST_TEST(testing::equals(waveform2.sample(0.0, interpolationOrder)(0, 0), 4.0));
  BOOST_TEST(testing::equals(waveform2.sample(0.5, interpolationOrder)(0, 0), 4.25));
  BOOST_TEST(testing::equals(waveform2.sample(1.0, interpolationOrder)(0, 0), 4.0));

  value(0) = 8.0;
  waveform2.store(value);

  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 8.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 4.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 2.0));

  BOOST_TEST(testing::equals(waveform2.sample(0.0, interpolationOrder)(0, 0), 4.0));
  BOOST_TEST(testing::equals(waveform2.sample(0.5, interpolationOrder)(0, 0), 5.75));
  BOOST_TEST(testing::equals(waveform2.sample(1.0, interpolationOrder)(0, 0), 8.0));
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

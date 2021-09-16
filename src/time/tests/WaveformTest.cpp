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
  int      timeWindowCounter  = 1;
  Waveform waveform(1, extrapolationOrder);
  BOOST_TEST(waveform.numberOfSamples() == 2);
  BOOST_TEST(waveform.numberOfData() == 1);
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 0), 0.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 1), 0.0));

  Eigen::VectorXd value(1);
  value(0) = 1.0;
  waveform.store(value);
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 1), 0.0));
  timeWindowCounter++;
  waveform.moveToNextWindow(timeWindowCounter, extrapolationOrder); // only applies constant extrapolation in second window
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 1), 1.0));

  value(0) = 4.0;
  waveform.store(value);
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 0), 4.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 1), 1.0));
  timeWindowCounter++;
  waveform.moveToNextWindow(timeWindowCounter, extrapolationOrder); // applies first order extrapolation in third window
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 0), 7.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 1), 4.0));

  // Test second order extrapolation
  extrapolationOrder = 2;
  timeWindowCounter  = 1;
  Waveform waveform2(1, extrapolationOrder);
  BOOST_TEST(waveform2.numberOfSamples() == 3);
  BOOST_TEST(waveform2.numberOfData() == 1);
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 0.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 0.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 0.0));

  value(0) = 1.0;
  waveform2.store(value);
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 0.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 0.0));
  timeWindowCounter++;
  waveform2.moveToNextWindow(timeWindowCounter, extrapolationOrder); // only applies constant extrapolation in second window
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 1.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 0.0));

  value(0) = 4.0;
  waveform2.store(value);
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 4.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 1.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 0.0));
  timeWindowCounter++;
  waveform2.moveToNextWindow(timeWindowCounter, extrapolationOrder); // only applies first order extrapolation in third window
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 7.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 4.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 1.0));

  value(0) = 8.0;
  waveform2.store(value);
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 8.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 4.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 1.0));
  timeWindowCounter++;
  waveform2.moveToNextWindow(timeWindowCounter, extrapolationOrder);    // applies second order extrapolation in fourth window
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 12.5)); // 12.5 = 2.5 * 8 - 2 * 4 + 0.5 * 1
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 8.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 4.0));
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

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

  double dt = 1.0;

  // Test first order extrapolation
  int      extrapolationOrder = 1;
  Waveform waveform(1, extrapolationOrder + 2);
  BOOST_TEST(waveform.lastTimeWindows().cols() == 3);
  BOOST_TEST(waveform.lastTimeWindows().rows() == 1);
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 0), 0.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 1), 0.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 2), 0.0));

  Eigen::VectorXd value(1);
  value(0) = 1.0;
  waveform.updateThisWindow(value);
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 1), 0.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 2), 0.0));
  waveform.moveToNextWindow();
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 1), 1.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 2), 0.0));
  auto sample = waveform.extrapolateData(extrapolationOrder, 2);
  BOOST_TEST(testing::equals(sample(0), 2.0));
  waveform.updateThisWindow(sample);
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 1), 1.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 2), 0.0));

  value(0) = 4.0;
  waveform.updateThisWindow(value);
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 0), 4.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 1), 1.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 2), 0.0));
  waveform.moveToNextWindow();
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 0), 4.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 1), 4.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 2), 1.0));
  sample = waveform.extrapolateData(extrapolationOrder, 3);
  BOOST_TEST(testing::equals(sample(0), 7.0));
  waveform.updateThisWindow(sample);
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 0), 7.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 1), 4.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 2), 1.0));

  // Test second order extrapolation
  extrapolationOrder = 2;
  Waveform waveform2(1, extrapolationOrder + 2);
  BOOST_TEST(waveform2.lastTimeWindows().cols() == 4);
  BOOST_TEST(waveform2.lastTimeWindows().rows() == 1);
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 0.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 0.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 0.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 3), 0.0));

  value(0) = 1.0;
  waveform2.updateThisWindow(value);
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 0.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 0.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 3), 0.0));
  waveform2.moveToNextWindow();
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 1.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 0.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 3), 0.0));
  sample = waveform2.extrapolateData(extrapolationOrder, 2);
  BOOST_TEST(testing::equals(sample(0), 2.0));
  waveform2.updateThisWindow(sample);
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 2.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 1.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 0.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 3), 0.0));

  value(0) = 4.0;
  waveform2.updateThisWindow(value);
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 4.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 1.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 0.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 3), 0.0));
  waveform2.moveToNextWindow();
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 4.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 4.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 1.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 3), 0.0));
  sample = waveform2.extrapolateData(extrapolationOrder, 3);
  BOOST_TEST(testing::equals(sample(0), 8.0));
  waveform2.updateThisWindow(sample);
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 0), 8.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 1), 4.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 2), 1.0));
  BOOST_TEST(testing::equals(waveform2.lastTimeWindows()(0, 3), 0.0));
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
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
  Waveform waveform           = Waveform();
  waveform.initialize(1, extrapolationOrder + 1);
  BOOST_TEST(waveform.lastTimeWindows().cols() == 2);
  BOOST_TEST(waveform.lastTimeWindows().rows() == 1);
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 0), 0.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 1), 0.0));

  Eigen::VectorXd value(1);
  value(0) = 1.0;
  waveform.addNewWindowData(value);
  auto sample = waveform.extrapolateData(extrapolationOrder, 2);
  BOOST_TEST(testing::equals(sample(0), 2.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 1), 0.0));

  value(0) = 4.0;
  waveform.addNewWindowData(value);
  sample = waveform.extrapolateData(extrapolationOrder, 3);
  BOOST_TEST(testing::equals(sample(0), 7.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 0), 4.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 1), 1.0));

  // Test second order extrapolation
  extrapolationOrder = 2;
  waveform.initialize(1, extrapolationOrder + 1);
  BOOST_TEST(waveform.lastTimeWindows().cols() == 3);
  BOOST_TEST(waveform.lastTimeWindows().rows() == 1);
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 0), 0.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 1), 0.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 2), 0.0));

  value(0) = 1.0;
  waveform.addNewWindowData(value);
  sample = waveform.extrapolateData(extrapolationOrder, 2);
  BOOST_TEST(testing::equals(sample(0), 2.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 0), 1.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 1), 0.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 2), 0.0));

  value(0) = 4.0;
  waveform.addNewWindowData(value);
  sample = waveform.extrapolateData(extrapolationOrder, 3);
  BOOST_TEST(testing::equals(sample(0), 8.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 0), 4.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 1), 1.0));
  BOOST_TEST(testing::equals(waveform.lastTimeWindows()(0, 2), 0.0));
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
#include <Eigen/Core>
#include "mesh/Data.hpp"
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
  const int       interpolationDegree = 0;
  const int       valuesSize          = 1;
  mesh::PtrData   dataPtr             = std::make_shared<mesh::Data>(mesh::Data("name", -1, valuesSize, 1));
  Eigen::VectorXd value(1);
  value(0) = 0.0;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  Waveform waveform(interpolationDegree, dataPtr);

  testing::WaveformFixture fixture;
  BOOST_TEST(fixture.valuesSize(waveform) == valuesSize);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 0.0));
}

BOOST_AUTO_TEST_CASE(testInitializationVector)
{
  PRECICE_TEST(1_rank);

  const int       interpolationDegree = 0;
  const int       valuesSize          = 3;
  mesh::PtrData   dataPtr             = std::make_shared<mesh::Data>(mesh::Data("name", -1, valuesSize, 1));
  Eigen::VectorXd value(valuesSize);
  value << 0, 0, 0;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  Waveform waveform(interpolationDegree, dataPtr);

  testing::WaveformFixture fixture;
  BOOST_TEST(fixture.valuesSize(waveform) == valuesSize);

  for (int i = 0; i < valuesSize; i++) {
    BOOST_TEST(testing::equals(waveform.sample(0.0)(i), 0.0));
    BOOST_TEST(testing::equals(waveform.sample(1.0)(i), 0.0));
  }
}

BOOST_AUTO_TEST_SUITE(InterpolationTests)

BOOST_AUTO_TEST_CASE(testInterpolateDataZerothDegree)
{
  PRECICE_TEST(1_rank);

  testing::WaveformFixture fixture;

  // Test zeroth degree interpolation
  const int       interpolationDegree = 0;
  const int       valuesSize          = 1;
  mesh::PtrData   dataPtr             = std::make_shared<mesh::Data>(mesh::Data("name", -1, valuesSize, 1));
  Eigen::VectorXd value(1);
  value(0) = 0.0;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_START, time::Sample{1, value});
  value(0) = 1.0;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});
  Waveform waveform(interpolationDegree, dataPtr);

  BOOST_TEST(fixture.valuesSize(waveform) == valuesSize);
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);

  dataPtr->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 1.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 1.0));

  value(0) = 2.0;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 2.0));

  dataPtr->timeStepsStorage().move();
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 2.0));

  value(0) = 3.0;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 3.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 3.0));
}

BOOST_AUTO_TEST_CASE(testInterpolateDataFirstDegree)
{
  PRECICE_TEST(1_rank);

  testing::WaveformFixture fixture;

  // Test first degree interpolation
  const int       interpolationDegree = 1;
  const int       valuesSize          = 1;
  mesh::PtrData   dataPtr             = std::make_shared<mesh::Data>(mesh::Data("name", -1, valuesSize, 1));
  Eigen::VectorXd value(1);
  value(0) = 0.0;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_START, time::Sample{1, value});
  value(0) = 1.0;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  Waveform waveform(interpolationDegree, dataPtr);

  BOOST_TEST(fixture.valuesSize(waveform) == valuesSize);
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 0.5));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 1.0));

  value(0) = 2.0;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 1.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 2.0));

  dataPtr->timeStepsStorage().move();
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 2.0));

  value(0) = 3.0;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 2.5));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 3.0));
}

// Remove or modify this feature? Creating a second degree interpolant by using data from previous windows is difficult, because this would require several pieces of data during initialization. What would be useful: Generating a second degree interpolant from multiple samples in a single window (if available). This would go into the least-squares direction
BOOST_AUTO_TEST_CASE(testInterpolateDataSecondDegree)
{
  PRECICE_TEST(1_rank);

  testing::WaveformFixture fixture;

  // Test second degree interpolation, but there are not enough samples. Therefore, always only first degree.
  const int       interpolationDegree = 2;
  const int       valuesSize          = 1;
  mesh::PtrData   dataPtr             = std::make_shared<mesh::Data>(mesh::Data("name", -1, valuesSize, 1));
  Eigen::VectorXd value(1);
  value(0) = 0.0;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_START, time::Sample{1, value});
  value(0) = 1.0;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  Waveform waveform(interpolationDegree, dataPtr);

  BOOST_TEST(fixture.valuesSize(waveform) == valuesSize);
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 0.5));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 1.0));

  dataPtr->timeStepsStorage().trim();

  value(0) = 2.0;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 1.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 2.0));

  dataPtr->timeStepsStorage().move();
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 2.0));

  value(0) = 8.0;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 5.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 8.0));

  dataPtr->timeStepsStorage().trim();

  value(0) = 4.0;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 3.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 4.0));

  dataPtr->timeStepsStorage().move();
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 4.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 4.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 4.0));

  value(0) = 8.0;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 4.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 6.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 8.0));
}

BOOST_AUTO_TEST_CASE(testInterpolateDataFirstDegreeVector)
{
  PRECICE_TEST(1_rank);

  testing::WaveformFixture fixture;

  // Test first degree interpolation
  const int       interpolationDegree = 1;
  const int       valuesSize          = 3;
  mesh::PtrData   dataPtr             = std::make_shared<mesh::Data>(mesh::Data("name", -1, valuesSize, 1));
  Eigen::VectorXd value(valuesSize);
  value << 0, 0, 0;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_START, time::Sample{1, value});
  value << 1, 2, 3;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  Waveform waveform(interpolationDegree, dataPtr);

  BOOST_TEST(fixture.valuesSize(waveform) == valuesSize);
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);

  for (int i = 0; i < valuesSize; i++) {
    BOOST_TEST(testing::equals(waveform.sample(0.0)(i), 0 * value[i]));
    BOOST_TEST(testing::equals(waveform.sample(0.5)(i), 0.5 * value[i]));
    BOOST_TEST(testing::equals(waveform.sample(1.0)(i), value[i]));
  }

  dataPtr->timeStepsStorage().trim();

  value << 2, 4, 2;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  for (int i = 0; i < valuesSize; i++) {
    BOOST_TEST(testing::equals(waveform.sample(0.0)(i), 0 * value[i]));
    BOOST_TEST(testing::equals(waveform.sample(0.5)(i), 0.5 * value[i]));
    BOOST_TEST(testing::equals(waveform.sample(1.0)(i), value[i]));
  }

  dataPtr->timeStepsStorage().move();
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);

  for (int i = 0; i < valuesSize; i++) {
    BOOST_TEST(testing::equals(waveform.sample(0.0)(i), value[i]));
    BOOST_TEST(testing::equals(waveform.sample(0.5)(i), value[i]));
    BOOST_TEST(testing::equals(waveform.sample(1.0)(i), value[i]));
  }

  Eigen::VectorXd value0 = value;
  value << 1, 2, 3;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  for (int i = 0; i < valuesSize; i++) {
    BOOST_TEST(testing::equals(waveform.sample(0.0)(i), value0[i]));
    BOOST_TEST(testing::equals(waveform.sample(0.5)(i), 0.5 * value0[i] + 0.5 * value[i]));
  }
}

BOOST_AUTO_TEST_CASE(testPiecewiseInterpolateDataZerothDegree)
{
  PRECICE_TEST(1_rank);

  testing::WaveformFixture fixture;

  // Test zeroth degree interpolation
  const int       interpolationDegree = 0;
  const int       valuesSize          = 1;
  mesh::PtrData   dataPtr             = std::make_shared<mesh::Data>(mesh::Data("name", -1, valuesSize, 1));
  Eigen::VectorXd value(1);
  value(0) = 0.0;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_START, time::Sample{1, value});
  value(0) = 0.5;
  dataPtr->setSampleAtTime(0.5 * time::Storage::WINDOW_END, time::Sample{1, value});
  value(0) = 1.0;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  Waveform waveform(interpolationDegree, dataPtr);

  BOOST_TEST(fixture.valuesSize(waveform) == valuesSize);
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 3);
  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 0.5));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 0.5));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 1.0));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 1.0));

  dataPtr->timeStepsStorage().trim();

  value(0) = 1.5;
  dataPtr->setSampleAtTime(0.5 * time::Storage::WINDOW_END, time::Sample{1, value});

  value(0) = 2.0;
  dataPtr->setSampleAtTime(1.0 * time::Storage::WINDOW_END, time::Sample{1, value});
  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 1.5));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 1.5));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 2.0));

  dataPtr->timeStepsStorage().move();
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);

  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 2.0));

  value(0) = 3.0;
  dataPtr->setSampleAtTime(1.0 * time::Storage::WINDOW_END, time::Sample{1, value});
  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 3.0));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 3.0));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 3.0));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 3.0));

  dataPtr->timeStepsStorage().trim();

  value(0) = 1.5;
  dataPtr->setSampleAtTime(0.5 * time::Storage::WINDOW_END, time::Sample{1, value});

  value(0) = 4.0;
  dataPtr->setSampleAtTime(1.0 * time::Storage::WINDOW_END, time::Sample{1, value});
  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 1.5));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 1.5));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 4.0));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 4.0));
}

BOOST_AUTO_TEST_CASE(testPiecewiseInterpolateDataFirstDegree)
{
  PRECICE_TEST(1_rank);

  testing::WaveformFixture fixture;

  // Test zeroth degree interpolation
  const int       interpolationDegree = 1;
  const int       valuesSize          = 1;
  mesh::PtrData   dataPtr             = std::make_shared<mesh::Data>(mesh::Data("name", -1, valuesSize, 1));
  Eigen::VectorXd value(1);
  value(0) = 0.0;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_START, time::Sample{1, value});
  value(0) = 0.5;
  dataPtr->setSampleAtTime(0.5 * time::Storage::WINDOW_END, time::Sample{1, value});
  value(0) = 1.0;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  Waveform waveform(interpolationDegree, dataPtr);

  BOOST_TEST(fixture.valuesSize(waveform) == valuesSize);
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 3);
  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 0.00));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 0.25));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 0.50));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 0.75));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 1.00));

  value(0) = 1.5;
  dataPtr->setSampleAtTime(0.5 * time::Storage::WINDOW_END, time::Sample{1, value});

  value(0) = 2.0;
  dataPtr->setSampleAtTime(1.0 * time::Storage::WINDOW_END, time::Sample{1, value});
  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 0.00));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 0.75));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 1.50));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 1.75));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 2.00));

  dataPtr->timeStepsStorage().move();
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 2);

  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 2.00));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 2.00));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 2.00));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 2.00));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 2.00));

  value(0) = 3.0;
  dataPtr->setSampleAtTime(1.0 * time::Storage::WINDOW_END, time::Sample{1, value});
  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 2.00));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 2.25));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 2.50));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 2.75));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 3.00));

  dataPtr->timeStepsStorage().trim();

  value(0) = 1.5;
  dataPtr->setSampleAtTime(0.5 * time::Storage::WINDOW_END, time::Sample{1, value});

  value(0) = 4.0;
  dataPtr->setSampleAtTime(1.0 * time::Storage::WINDOW_END, time::Sample{1, value});
  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 2.00));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 1.75));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 1.50));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 2.75));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 4.00));
}

BOOST_AUTO_TEST_CASE(testPiecewiseInterpolateDataSecondDegree)
{
  PRECICE_TEST(1_rank);

  testing::WaveformFixture fixture;

  // Test zeroth degree interpolation
  const int       interpolationDegree = 2;
  const int       valuesSize          = 1;
  mesh::PtrData   dataPtr             = std::make_shared<mesh::Data>(mesh::Data("name", -1, valuesSize, 1));
  Eigen::VectorXd value(1);
  value(0) = 0.0;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_START, time::Sample{1, value});
  value(0) = 0.0;
  dataPtr->setSampleAtTime(0.5 * time::Storage::WINDOW_END, time::Sample{1, value});
  value(0) = 2.0;
  dataPtr->setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  Waveform waveform(interpolationDegree, dataPtr);

  BOOST_TEST(fixture.valuesSize(waveform) == valuesSize);
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 3);
  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 0.00));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), -0.25));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 0.00));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 0.75));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 2.00));
}

BOOST_AUTO_TEST_CASE(testPiecewiseInterpolateDataThirdDegree)
{
  PRECICE_TEST(1_rank);

  testing::WaveformFixture fixture;

  // Test zeroth degree interpolation
  const int     interpolationDegree = 3;
  const int     valuesSize          = 1;
  mesh::PtrData dataPtr             = std::make_shared<mesh::Data>(mesh::Data("name", -1, valuesSize, 1));

  Waveform waveform(interpolationDegree, dataPtr);

  // linearly increasing values
  Eigen::VectorXd value(1);
  for (double t : std::vector<double>{0, 0.25, 0.5, 0.75, 1}) {
    value(0) = t;
    dataPtr->setSampleAtTime(t * time::Storage::WINDOW_END, time::Sample{1, value});
  }

  BOOST_TEST(fixture.valuesSize(waveform) == valuesSize);
  BOOST_TEST(fixture.numberOfStoredSamples(waveform) == 5);

  for (double t : std::vector<double>{0.1, 0.2, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t));
  }

  dataPtr->timeStepsStorage().trim();

  // quadratically increasing values
  for (double t : std::vector<double>{0, 0.25, 0.5, 0.75, 1}) {
    value(0) = t * t;
    dataPtr->setSampleAtTime(t * time::Storage::WINDOW_END, time::Sample{1, value});
  }

  // interpolates given values
  for (double t : std::vector<double>{0, 0.25, 0.5, 0.75, 1}) {
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t * t));
  }

  // introduces no approximation error w.r.t function
  for (double t : std::vector<double>{0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t * t));
  }

  dataPtr->timeStepsStorage().trim();

  // cubically increasing values
  for (double t : std::vector<double>{0, 0.25, 0.5, 0.75, 1}) {
    value(0) = t * t * t;
    dataPtr->setSampleAtTime(t * time::Storage::WINDOW_END, time::Sample{1, value});
  }

  // interpolates given values
  for (double t : std::vector<double>{0, 0.25, 0.5, 0.75, 1}) {
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t * t * t));
  }

  // introduces no approximation error w.r.t function
  for (double t : std::vector<double>{0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t * t * t));
  }

  dataPtr->timeStepsStorage().trim();

  // cubically increasing values, but with non-uniform spacing
  for (double t : std::vector<double>{0, 0.01, 0.1, 0.2, 1}) {
    value(0) = t * t * t;
    dataPtr->setSampleAtTime(t * time::Storage::WINDOW_END, time::Sample{1, value});
  }

  // interpolates given values
  for (double t : std::vector<double>{0, 0.01, 0.1, 0.2, 1}) {
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t * t * t));
  }

  // introduces no approximation error w.r.t function
  for (double t : std::vector<double>{0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t * t * t));
  }

  dataPtr->timeStepsStorage().trim();

  // quadratically increasing values, but with non-uniform spacing
  for (double t : std::vector<double>{0, 0.25, 0.5, 0.75, 1}) {
    value(0) = t * t * t * t;
    dataPtr->setSampleAtTime(t * time::Storage::WINDOW_END, time::Sample{1, value});
  }

  // interpolates given values
  for (double t : std::vector<double>{0, 0.25, 0.5, 0.75, 1}) {
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t * t * t * t));
  }

  // introduces approximation error w.r.t function
  for (double t : std::vector<double>{0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {
    double tol = 0.015625; // error < h**3 = 0.015625
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t * t * t * t, tol));
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

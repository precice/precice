#include <Eigen/Core>
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "testing/ExtrapolationFixture.hpp"
#include "cplscheme/impl/Extrapolation.hpp"

using namespace precice;
using namespace precice::cplscheme;

BOOST_AUTO_TEST_SUITE(CplSchemeTests)

BOOST_AUTO_TEST_SUITE(ExtrapolationTests)
BOOST_AUTO_TEST_CASE(testExtrapolateDataFirstOrder)
{
  PRECICE_TEST(1_rank);

  testing::ExtrapolationFixture fixture;

  // Test first order extrapolation
  const int extrapolationOrder = 1;
  Extrapolation extrapolation(extrapolationOrder);
  const int valuesSize = 1;
  extrapolation.initialize(valuesSize);
  BOOST_TEST(fixture.sizeOfSampleStorage(extrapolation) == 2);
  BOOST_TEST(fixture.numberOfStoredSamples(extrapolation) == 1);
  BOOST_TEST(fixture.valuesSize(extrapolation) == 1);

  // use zero initial data
  extrapolation.moveToNextWindow();
  BOOST_TEST(fixture.numberOfStoredSamples(extrapolation) == 2);
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 0), 0.0));
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 1), 0.0));

  Eigen::VectorXd value(1);
  value(0) = 1.0;
  extrapolation.store(value);
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 0), 1.0));
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 1), 0.0));
  BOOST_TEST(fixture.numberOfStoredSamples(extrapolation) == 2);
  extrapolation.moveToNextWindow(); // applies first order extrapolation in second window
  BOOST_TEST(fixture.numberOfStoredSamples(extrapolation) == 2);
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 0), 2.0)); // = 2*1 - 0
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 1), 1.0));

  value(0) = 4.0;
  extrapolation.store(value);
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 0), 4.0));
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 1), 1.0));
  BOOST_TEST(fixture.numberOfStoredSamples(extrapolation) == 2);
  extrapolation.moveToNextWindow(); // applies first order extrapolation in third window
  BOOST_TEST(fixture.numberOfStoredSamples(extrapolation) == 2);
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 0), 7.0)); // = 2*4 - 1
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 1), 4.0));

  value(0) = 8.0;
  extrapolation.store(value);
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 0), 8.0));
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 1), 4.0));
  BOOST_TEST(fixture.numberOfStoredSamples(extrapolation) == 2);
  extrapolation.moveToNextWindow(); // applies first order extrapolation in forth window
  BOOST_TEST(fixture.numberOfStoredSamples(extrapolation) == 2);
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 0), 12.0)); // 10.0 = 2 * 8 - 4
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 1), 8.0));
}

BOOST_AUTO_TEST_CASE(testExtrapolateDataSecondOrder)
{
  PRECICE_TEST(1_rank);

  testing::ExtrapolationFixture fixture;

  // Test second order extrapolation
  const int extrapolationOrder = 2;
  Extrapolation extrapolation(extrapolationOrder);
  const int valuesSize = 1;
  extrapolation.initialize(valuesSize);
  BOOST_TEST(fixture.sizeOfSampleStorage(extrapolation) == 3);
  BOOST_TEST(fixture.numberOfStoredSamples(extrapolation) == 1);
  BOOST_TEST(fixture.valuesSize(extrapolation) == 1);

  // use zero initial data
  extrapolation.moveToNextWindow();
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 0), 0.0));
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 1), 0.0));
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 2), 0.0));

  Eigen::VectorXd value(1);
  value(0) = 1.0;
  extrapolation.store(value);
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 0), 1.0));
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 1), 0.0));
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 2), 0.0));
  BOOST_TEST(fixture.numberOfStoredSamples(extrapolation) == 2);
  extrapolation.moveToNextWindow(); // applies first order extrapolation in second window
  BOOST_TEST(fixture.numberOfStoredSamples(extrapolation) == 3);
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 0), 2.0)); // = 2*1 - 0
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 1), 1.0));
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 2), 0.0));

  value(0) = 4.0;
  extrapolation.store(value);
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 0), 4.0));
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 1), 1.0));
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 2), 0.0));
  BOOST_TEST(fixture.numberOfStoredSamples(extrapolation) == 3);
  extrapolation.moveToNextWindow(); // applies second order extrapolation in third window
  BOOST_TEST(fixture.numberOfStoredSamples(extrapolation) == 3);
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 0), 8.0)); // = 2.5*4 - 2 * 1 + 0.5 * 0
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 1), 4.0));
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 2), 1.0));

  value(0) = 8.0;
  extrapolation.store(value);
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 0), 8.0));
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 1), 4.0));
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 2), 1.0));
  BOOST_TEST(fixture.numberOfStoredSamples(extrapolation) == 3);
  extrapolation.moveToNextWindow(); // applies second order extrapolation in fourth window
  BOOST_TEST(fixture.numberOfStoredSamples(extrapolation) == 3);
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 0), 12.5)); // = 2.5 * 8 - 2 * 4 + 0.5 * 1
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 1), 8.0));
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 2), 4.0));

  value(0) = 16.0;
  extrapolation.store(value);
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 0), 16.0));
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 1), 8.0));
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 2), 4.0));
  BOOST_TEST(fixture.numberOfStoredSamples(extrapolation) == 3);
  extrapolation.moveToNextWindow(); // applies second order extrapolation in fifth window
  BOOST_TEST(fixture.numberOfStoredSamples(extrapolation) == 3);
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 0), 26.0)); // = 2.5 * 16.0 - 2 * 8 + 0.5 * 4
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 1), 16.0));
  BOOST_TEST(testing::equals(fixture.getValue(extrapolation, 0, 2), 8.0));
}
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
#include <Eigen/Core>
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "time/Storage.hpp"

using namespace precice;
using namespace precice::time;

BOOST_AUTO_TEST_SUITE(TimeTests)
BOOST_AUTO_TEST_SUITE(StorageTests)

// create storage and test for correct initial values.
BOOST_AUTO_TEST_CASE(testInitialize)
{
  PRECICE_TEST(1_rank);
  auto storage = Storage();
  int  nValues = 3;
  BOOST_TEST(storage.nTimes() == 0);
  storage.initialize(time::Sample{Eigen::VectorXd::Ones(nValues)});
  BOOST_TEST(storage.nDofs() == nValues);
  BOOST_TEST(storage.nTimes() == 2);
  for (int i = 0; i < nValues; i++) {
    BOOST_TEST(storage.getValuesAtOrAfter(0)(i) == 1);
    BOOST_TEST(storage.getValuesAtOrAfter(0.5)(i) == 1);
    BOOST_TEST(storage.getValuesAtOrAfter(1)(i) == 1);
  }
}

// create storage and trim it.
BOOST_AUTO_TEST_CASE(testClear)
{
  PRECICE_TEST(1_rank);
  auto storage = Storage();
  int  nValues = 3;
  BOOST_TEST(storage.nTimes() == 0);
  storage.initialize(time::Sample{Eigen::VectorXd::Ones(nValues)});
  BOOST_TEST(storage.nDofs() == nValues);
  BOOST_TEST(storage.nTimes() == 2);
  BOOST_TEST(storage.maxStoredNormalizedDt() == 1.0);
  storage.trim();
  BOOST_TEST(storage.nDofs() == nValues);
  BOOST_TEST(storage.nTimes() == 1);
  BOOST_TEST(storage.maxStoredNormalizedDt() == 0.0);
}

// create storage, add some values and then move to next window.
BOOST_AUTO_TEST_CASE(testMove)
{
  PRECICE_TEST(1_rank);
  auto storage = Storage();
  int  nValues = 3;
  BOOST_TEST(storage.nTimes() == 0);
  storage.initialize(time::Sample{Eigen::VectorXd::Ones(nValues)});
  BOOST_TEST(storage.nDofs() == nValues);
  BOOST_TEST(storage.nTimes() == 2);
  BOOST_TEST(storage.maxStoredNormalizedDt() == 1.0);
  storage.trim();
  BOOST_TEST(storage.nTimes() == 1);
  storage.setSampleAtTime(0.5, Sample{Eigen::VectorXd::Ones(nValues)});
  BOOST_TEST(storage.nTimes() == 2);
  BOOST_TEST(storage.maxStoredNormalizedDt() == 0.5);
  storage.setSampleAtTime(1.0, Sample{Eigen::VectorXd::Zero(nValues)});
  BOOST_TEST(storage.nTimes() == 3);
  BOOST_TEST(storage.maxStoredNormalizedDt() == 1.0);
  for (int i = 0; i < nValues; i++) {
    BOOST_TEST(storage.getValuesAtOrAfter(0)(i) == 1);
    BOOST_TEST(storage.getValuesAtOrAfter(0.5)(i) == 1);
    BOOST_TEST(storage.getValuesAtOrAfter(1)(i) == 0);
  }
  storage.move();
  BOOST_TEST(storage.nDofs() == nValues);
  BOOST_TEST(storage.nTimes() == 2);
  BOOST_TEST(storage.maxStoredNormalizedDt() == 1.0);
  for (int i = 0; i < nValues; i++) {
    BOOST_TEST(storage.getValuesAtOrAfter(0)(i) == 0);
    BOOST_TEST(storage.getValuesAtOrAfter(1)(i) == 0);
  }
}

// get times and values
BOOST_AUTO_TEST_CASE(testGetTimesAndValues)
{
  PRECICE_TEST(1_rank);
  auto storage = Storage();
  int  nValues = 3;
  storage.initialize(time::Sample{Eigen::VectorXd::Ones(nValues)});
  storage.trim();
  storage.setSampleAtTime(0.5, Sample{Eigen::VectorXd::Ones(nValues)});
  storage.setSampleAtTime(1.0, Sample{Eigen::VectorXd::Zero(nValues)});
  auto times = storage.getTimes();
  BOOST_TEST(times[0] == 0.0);
  BOOST_TEST(times[1] == 0.5);
  BOOST_TEST(times[2] == 1.00);
  auto timesAndValues = storage.getTimesAndValues();
  BOOST_TEST(timesAndValues.first[0] == 0.0);
  BOOST_TEST(timesAndValues.first[1] == 0.5);
  BOOST_TEST(timesAndValues.first[2] == 1.00);
  for (int i = 0; i < nValues; i++) {
    BOOST_TEST(timesAndValues.second.col(0)(i) == 1);
    BOOST_TEST(timesAndValues.second.col(1)(i) == 1);
    BOOST_TEST(timesAndValues.second.col(2)(i) == 0);
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

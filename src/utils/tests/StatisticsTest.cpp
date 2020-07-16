#include <cmath>
#include "math/constants.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/Statistics.hpp"
#include "utils/String.hpp"

using namespace precice;
namespace pu = precice::utils;

BOOST_AUTO_TEST_SUITE(UtilsTests)

BOOST_AUTO_TEST_CASE(DistanceAccumulator)
{
  PRECICE_TEST(1_rank);
  pu::statistics::DistanceAccumulator acc;
  acc(0.01);
  BOOST_TEST(!acc.empty());
  BOOST_TEST(acc.count() == 1);
  BOOST_TEST(acc.min() == 0.01);
  BOOST_TEST(acc.max() == 0.01);
  BOOST_TEST(acc.mean() == 0.01);
  BOOST_TEST(acc.variance() == 0);

  acc(1);
  acc(-1);
  acc(23);
  acc(11);
  BOOST_TEST(acc.min() == -1);
  BOOST_TEST(acc.min() != 0);
  BOOST_TEST(acc.max() == 23);
}

BOOST_AUTO_TEST_CASE(DistanceAccumulatorOnEmptyMesh)
{
  pu::statistics::DistanceAccumulator acc;
  BOOST_TEST(acc.empty());
  BOOST_TEST(acc.count() == 0);
  BOOST_TEST(std::isnan(acc.min()));
  BOOST_TEST(std::isnan(acc.mean()));
  BOOST_TEST(std::isnan(acc.max()));
  BOOST_TEST(std::isnan(acc.variance()));
}

BOOST_AUTO_TEST_SUITE_END()

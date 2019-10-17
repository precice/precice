#include "testing/Testing.hpp"
#include "utils/Statistics.hpp"
#include <Eigen/Core>

namespace pu = precice::utils;

BOOST_AUTO_TEST_SUITE(UtilsTests)

BOOST_AUTO_TEST_CASE(DistanceAccumulator)
{
    pu::statistics::DistanceAccumulator acc;
    acc(0);
    BOOST_TEST(acc.min() == 0);
    BOOST_TEST(acc.max() == 0);

    acc(1);
    acc(-1);
    acc(23);
    acc(11);
    BOOST_TEST(acc.min() == -1);
    BOOST_TEST(acc.min() != 0);
    BOOST_TEST(acc.max() == 23);
}

BOOST_AUTO_TEST_SUITE_END()

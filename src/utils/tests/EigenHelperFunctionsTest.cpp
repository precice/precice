#include "testing/Testing.hpp"
#include "utils/EigenHelperFunctions.hpp"

using namespace precice::utils;

BOOST_AUTO_TEST_SUITE(UtilsTests)

BOOST_AUTO_TEST_CASE(FirstN)
{
    Eigen::VectorXd a(7);
    a << 1, 2, 3, 4, 5, 6, 7;
    Eigen::RowVectorXd b(3);
    b << 1, 2, 3;
    BOOST_TEST(firstN(a, 3) == b);
}

BOOST_AUTO_TEST_SUITE_END()

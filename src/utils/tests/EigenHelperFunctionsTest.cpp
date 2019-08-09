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

BOOST_AUTO_TEST_CASE(ComponentWiseLess)
{
    precice::utils::ComponentWiseLess cwl;

    Eigen::VectorXd a(8);
    a << 1, 2, 3, 4, 5, 6, 7, 9;

    Eigen::VectorXd b(8);
    b << 1, 2, 3, 4, 5, 6, 8, 0;

    BOOST_TEST(componentWiseLess(a, b));
    BOOST_TEST(!componentWiseLess(b, a));
    BOOST_TEST(cwl(a, b));
    BOOST_TEST(!cwl(b, a));

    Eigen::VectorXd c = b;

    BOOST_TEST(componentWiseLess(c, b));
    BOOST_TEST(componentWiseLess(b, c));
    BOOST_TEST(cwl(c, b));
    BOOST_TEST(cwl(b, c));
}


BOOST_AUTO_TEST_SUITE_END()

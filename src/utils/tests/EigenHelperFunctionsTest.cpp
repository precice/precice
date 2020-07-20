#include <Eigen/Core>
#include <iosfwd>
#include <string>
#include "math/constants.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/String.hpp"
#include "utils/algorithm.hpp"

using namespace precice;
using namespace precice::utils;

BOOST_AUTO_TEST_SUITE(UtilsTests)
BOOST_AUTO_TEST_SUITE(EigenHelperFunctionsTests)

BOOST_AUTO_TEST_CASE(FirstN)
{
  PRECICE_TEST(1_rank);
  Eigen::VectorXd a(7);
  a << 1, 2, 3, 4, 5, 6, 7;
  Eigen::RowVectorXd b(3);
  b << 1, 2, 3;
  BOOST_TEST(firstN(a, 3) == b);
}

BOOST_AUTO_TEST_SUITE(RangePreview)

BOOST_AUTO_TEST_CASE(EigenVector)
{
  PRECICE_TEST(1_rank);
  Eigen::VectorXd a{7};
  a << 1, 2, 3, 4, 5, 6, 0;
  std::ostringstream oss;
  oss << previewRange(2, a);
  std::string str{oss.str()};
  BOOST_TEST(str == "[1, 2, ... , 6, 0] min:0 max:6");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_CASE(ComponentWiseLess)
{
  PRECICE_TEST(1_rank);
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

  BOOST_TEST(!componentWiseLess(c, b));
  BOOST_TEST(!componentWiseLess(b, c));
  BOOST_TEST(!cwl(c, b));
  BOOST_TEST(!cwl(b, c));
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

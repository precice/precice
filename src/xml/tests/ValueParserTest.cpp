#include <Eigen/Core>
#include <stdexcept>
#include <string>
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "xml/ValueParser.hpp"

using namespace precice;
using namespace precice::xml;

BOOST_AUTO_TEST_SUITE(XML)
BOOST_AUTO_TEST_SUITE(ValueParser)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ParseDouble)
{
  PRECICE_TEST();
  double value = -1.0;

  readValueSpecific("1.0", value);
  BOOST_TEST(value == 1.0);
  readValueSpecific("3", value);
  BOOST_TEST(value == 3.0);
  readValueSpecific("-2.5", value);
  BOOST_TEST(value == -2.5);
  readValueSpecific("1e3", value);
  BOOST_TEST(value == 1000.0);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ParseDoubleFraction)
{
  PRECICE_TEST();
  double value = -1.0;

  readValueSpecific("1/2", value);
  BOOST_TEST(value == 0.5);
  readValueSpecific("-1/4", value);
  BOOST_TEST(value == -0.25);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ParseDoubleRejectsEmpty)
{
  PRECICE_TEST();
  double value = -1.0;

  // An empty or whitespace-only value used to leave `value` indeterminate without
  // raising an error, because the failed extraction also sets eofbit.
  BOOST_CHECK_THROW(readValueSpecific("", value), std::runtime_error);
  BOOST_CHECK_THROW(readValueSpecific("   ", value), std::runtime_error);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ParseDoubleRejectsMalformed)
{
  PRECICE_TEST();
  double value = -1.0;

  BOOST_CHECK_THROW(readValueSpecific("abc", value), std::runtime_error);
  BOOST_CHECK_THROW(readValueSpecific("1.0foo", value), std::runtime_error);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ParseDoubleRejectsMalformedFraction)
{
  PRECICE_TEST();
  double value = -1.0;

  // Each of these splits into at least one empty operand.
  BOOST_CHECK_THROW(readValueSpecific("1/", value), std::runtime_error);
  BOOST_CHECK_THROW(readValueSpecific("/2", value), std::runtime_error);
  BOOST_CHECK_THROW(readValueSpecific("/", value), std::runtime_error);
  BOOST_CHECK_THROW(readValueSpecific("1/2/3", value), std::runtime_error);
  // Used to yield inf and nan respectively.
  BOOST_CHECK_THROW(readValueSpecific("1/0", value), std::runtime_error);
  BOOST_CHECK_THROW(readValueSpecific("0/0", value), std::runtime_error);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ParseInt)
{
  PRECICE_TEST();
  int value = -1;

  readValueSpecific("42", value);
  BOOST_TEST(value == 42);
  readValueSpecific("-7", value);
  BOOST_TEST(value == -7);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ParseIntRejectsMalformed)
{
  PRECICE_TEST();
  int value = -1;

  BOOST_CHECK_THROW(readValueSpecific("", value), std::runtime_error);
  BOOST_CHECK_THROW(readValueSpecific("   ", value), std::runtime_error);
  BOOST_CHECK_THROW(readValueSpecific("abc", value), std::runtime_error);
  BOOST_CHECK_THROW(readValueSpecific("1.5", value), std::runtime_error);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ParseVector)
{
  PRECICE_TEST();
  Eigen::VectorXd value;

  readValueSpecific("3;2;1", value);
  BOOST_TEST(value.size() == 3);
  BOOST_TEST(value(0) == 3.0);
  BOOST_TEST(value(1) == 2.0);
  BOOST_TEST(value(2) == 1.0);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ParseVectorRejectsMalformedComponent)
{
  PRECICE_TEST();
  Eigen::VectorXd value;

  BOOST_CHECK_THROW(readValueSpecific("1;abc", value), std::runtime_error);
  // A whitespace-only component is a token of its own and used to be accepted silently.
  BOOST_CHECK_THROW(readValueSpecific("1; ;2", value), std::runtime_error);
  BOOST_CHECK_THROW(readValueSpecific("", value), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

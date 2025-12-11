#include "logging/LogMacros.hpp"
#include "logging/Logger.hpp"
#include "math/math.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(MathTests)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(IntegerExponent)
{
  PRECICE_TEST();
  BOOST_TEST(math::pow_int<4>(4.2) == std::pow(4.2, 4));
  BOOST_TEST(math::pow_int<3>(7.2) == std::pow(7.2, 3));
}
BOOST_AUTO_TEST_SUITE_END() // Math

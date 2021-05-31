#include <Eigen/Core>
#include "logging/LogMacros.hpp"
#include "math/constants.hpp"
#include "math/differences.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::math;

BOOST_AUTO_TEST_SUITE(MathTests)
BOOST_AUTO_TEST_SUITE(Differences)

BOOST_AUTO_TEST_CASE(Scalar)
{
  PRECICE_TEST(1_rank);
  double a   = 1.0;
  double b   = 2.0;
  double eps = 1e-14;

  BOOST_CHECK(greater(b, a, eps));
  BOOST_CHECK(not greater(a, a - eps, eps));
  BOOST_CHECK(greater(a, a - 10.0 * eps, eps));

  BOOST_CHECK(not greaterEquals(a, b, eps));
  BOOST_CHECK(greaterEquals(b, a, eps));
  BOOST_CHECK(greaterEquals(a, a, eps));
  BOOST_CHECK(greaterEquals(a, a + 0.1 * eps, eps));
  BOOST_CHECK(not greaterEquals(a, a + 10 * eps, eps));

  BOOST_CHECK(smaller(a, b, eps));
  BOOST_CHECK(not smaller(a, a + eps, eps));
  BOOST_CHECK(smaller(a, a + 10.0 * eps, eps));

  BOOST_CHECK(smallerEquals(a, b, eps));
  BOOST_CHECK(smallerEquals(a, a + eps, eps));
  BOOST_CHECK(smallerEquals(a + eps, a, eps));
  BOOST_CHECK(smallerEquals(a, a + 10.0 * eps, eps));

  BOOST_CHECK(not equals(a, b, eps));
  BOOST_CHECK(equals(a, a, eps));
  BOOST_CHECK(equals(a, a + eps, eps));
  BOOST_CHECK(not equals(a, a + 10.0 * eps, eps));
  BOOST_CHECK(equals(a, a + 10.0 * eps, 10.0 * eps));
}

BOOST_AUTO_TEST_CASE(Vector)
{
  PRECICE_TEST(1_rank);
  Eigen::Vector3d vec0(1.0, 2.0, 3.0);
  Eigen::Vector3d vec1(vec0);
  BOOST_CHECK(equals(vec0, vec1));
  BOOST_CHECK(not oneGreater(vec0, vec1));
  BOOST_CHECK(oneGreaterEquals(vec0, vec1));
  BOOST_CHECK(not allGreater(vec0, vec1));

  vec0 << 2.0, 2.0, 3.0;
  BOOST_CHECK(not equals(vec0, vec1));
  BOOST_CHECK(oneGreater(vec0, vec1));
  BOOST_CHECK(oneGreaterEquals(vec0, vec1));
  BOOST_CHECK(not allGreater(vec0, vec1));

  vec0 << 2.0, 3.0, 4.0;
  BOOST_CHECK(not equals(vec0, vec1));
  BOOST_CHECK(oneGreater(vec0, vec1));
  BOOST_CHECK(oneGreaterEquals(vec0, vec1));
  BOOST_CHECK(allGreater(vec0, vec1));

  // up to here vec1=vec0
  const double tolerance = 1e-14;
  vec0(0)                = vec1(0);
  vec0(1)                = vec1(1);
  vec0(2)                = vec1(2) + 0.99 * tolerance;
  BOOST_CHECK(equals(vec0, vec1, tolerance));
  BOOST_CHECK(not oneGreater(vec0, vec1, tolerance));
  BOOST_CHECK(oneGreaterEquals(vec0, vec1));
  BOOST_CHECK(not allGreater(vec0, vec1, tolerance));

  vec0(2) = vec1(2) + 10.0 * tolerance;
  BOOST_CHECK(not equals(vec0, vec1, tolerance));
  BOOST_CHECK(oneGreater(vec0, vec1, tolerance));
  BOOST_CHECK(oneGreaterEquals(vec0, vec1));
  BOOST_CHECK(not allGreater(vec0, vec1, tolerance));

  vec0 << 1.0, 2.0, 3.0;
  vec0 = vec0.array() + (10.0 * tolerance);
  BOOST_CHECK(not equals(vec0, vec1, tolerance));
  BOOST_CHECK(oneGreater(vec0, vec1, tolerance));
  BOOST_CHECK(oneGreaterEquals(vec0, vec1));
  BOOST_CHECK(allGreater(vec0, vec1, tolerance));

  vec0 << 1.0, 2.0, 3.0;
  vec0 = vec0.array() - 0.99 * tolerance;
  BOOST_CHECK(oneGreaterEquals(vec0, vec1));
}

BOOST_AUTO_TEST_SUITE_END() // Differences

BOOST_AUTO_TEST_SUITE_END() // Math

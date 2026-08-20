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

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(Scalar)
{
  PRECICE_TEST();
  double a   = 1.0;
  double b   = 2.0;
  double eps = 1e-14;

  // For operands near 1.0, scale factor == 1 so behaviour is unchanged.
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

/**
 * Regression test for Issue #2592:
 * "Coupling aborts deterministically at t = 2^6 = 64 s"
 *
 * The root cause is that waveform timestamps are compared with a FIXED absolute
 * tolerance (1e-14), but the IEEE 754 double ULP grows with magnitude:
 *
 *   [32,  64) s  -> ULP ~7.11e-15  (OK, ULP < 1e-14)
 *   [64, 128) s  -> ULP ~1.42e-14  (FAIL: ULP > 1e-14)
 *
 * After the fix, tolerances are scaled by max(1, |a|, |b|), so two timestamps
 * that differ by exactly 1 ULP are always considered equal, regardless of
 * absolute magnitude.
 */
PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ScalarLargeTimestamps)
{
  PRECICE_TEST();
  const double eps = math::NUMERICAL_ZERO_DIFFERENCE;

  // -----------------------------------------------------------------------
  // Reproduce the t = 64 s crossover (Issue #2592, first crash site).
  // Simulate accumulated window-start timestamps: windowStart is computed
  // exactly, but the received sample's timestamp is 1 ULP below.
  // Before the fix: math::smaller(sample, windowStart) == true -> waveform
  //                 trimmed -> crash.
  // After the fix:  the 1-ULP gap is within scaled tolerance -> NOT smaller.
  // -----------------------------------------------------------------------
  {
    const double t64    = 64.0;
    const double t64ulp = t64 - std::numeric_limits<double>::epsilon() * t64; // 1 ULP below

    // sample_timestamp < windowStart by 1 ULP — must NOT be considered smaller
    BOOST_CHECK_MESSAGE(
        not smaller(t64ulp, t64, eps),
        "smaller() returned true for a 1-ULP difference at t=64s: "
        "this would cause waveform trimming and abort (Issue #2592).");

    // They should be considered equal
    BOOST_CHECK_MESSAGE(
        equals(t64ulp, t64, eps),
        "equals() returned false for a 1-ULP difference at t=64s.");

    // greaterEquals: t64ulp should still be >= t64 within tolerance
    BOOST_CHECK_MESSAGE(
        greaterEquals(t64ulp, t64, eps),
        "greaterEquals() returned false for a 1-ULP difference at t=64s.");
  }

  // -----------------------------------------------------------------------
  // Reproduce the t = 256 s crossover (Issue #2592, second crash site).
  // timeToAdd = time - time_old carries ~ULP(t) error; the greaterEquals
  // check in addComputedTime() spuriously fails.
  // -----------------------------------------------------------------------
  {
    const double t256    = 256.0;
    const double t256ulp = t256 + std::numeric_limits<double>::epsilon() * t256; // 1 ULP above

    // maxStepSize should be >= a value that is 1 ULP larger (with scaled tol)
    BOOST_CHECK_MESSAGE(
        greaterEquals(t256, t256ulp, eps),
        "greaterEquals() returned false for a 1-ULP difference at t=256s: "
        "this would cause a spurious \"time step exceeds maximum\" abort (Issue #2592).");
  }

  // -----------------------------------------------------------------------
  // Sanity check: genuinely different values must still compare correctly.
  // -----------------------------------------------------------------------
  {
    const double t = 100.0;
    BOOST_CHECK(greater(t + 1.0, t, eps));
    BOOST_CHECK(smaller(t, t + 1.0, eps));
    BOOST_CHECK(not equals(t, t + 1.0, eps));
  }
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(Vector)
{
  PRECICE_TEST();
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

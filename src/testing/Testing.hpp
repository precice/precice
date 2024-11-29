#pragma once

#include <Eigen/Core>
#include <boost/test/unit_test.hpp>
#include <string>
#include <type_traits>

#include "math/differences.hpp"
#include "math/math.hpp"
#include "testing/TestContext.hpp"
#include "utils/IntraComm.hpp"

namespace precice::testing {

namespace bt = boost::unit_test;

DataID operator"" _dataID(unsigned long long n);

namespace inject {
using precice::testing::Require;
using precice::testing::operator""_rank;
using precice::testing::operator""_ranks;
using precice::testing::operator""_on;
using precice::testing::operator""_dataID;
} // namespace inject

/// Creates and attaches a TestSetup to a Boost test case
#define PRECICE_TEST_SETUP(...) \
  BOOST_TEST_DECORATOR(*precice::testing::testSetup(__VA_ARGS__))

/** Setup the preCICE test
 * Assigns a participant and rank from the TestSetup.
 * Unneeded ranks are marked as invalid and skip the test.
 * Provides a TestContext named context which can be used in the test.
 *
 * @pre PRECICE_TEST_SETUP has been called before the test to setup a test context.
 */
#define PRECICE_TEST()                                                     \
  precice::testing::TestContext context{precice::testing::getTestSetup()}; \
  if (context.invalid) {                                                   \
    return;                                                                \
  }                                                                        \
  BOOST_TEST_MESSAGE(context.describe());                                  \
  boost::unit_test::framework::add_context(BOOST_TEST_LAZY_MSG(context.describe()), true);

/// struct giving access to the impl of a befriended class or struct
struct WhiteboxAccessor {
  /** Returns the impl of the obj by reference.
     *
     * Returns a reference to the object pointed to by the _impl of a class.
     * This class needs to be friend of T.
     *
     * @param[in] obj The object to fetch the impl from.
     * @returns a lvalue reference to the impl object.
     */
  template <typename T>
  static auto impl(T &obj) -> typename std::add_lvalue_reference<decltype(*(obj._impl))>::type
  {
    return *obj._impl;
  }
};

/// equals to be used in tests. Compares two matrices using a given tolerance. Prints both operands on failure.
template <class DerivedA, class DerivedB>
boost::test_tools::predicate_result equals(const Eigen::MatrixBase<DerivedA> &A,
                                           const Eigen::MatrixBase<DerivedB> &B,
                                           double                             tolerance = math::NUMERICAL_ZERO_DIFFERENCE)
{
  auto approx = [tolerance](double a, double b) -> bool {
    if (std::max(std::abs(a), std::abs(b)) < tolerance) {
      return true;
    } else {
      return std::abs(a - b) <= tolerance * std::max(std::abs(a), std::abs(b));
    }
  };

  if (!A.binaryExpr(B, approx).all()) {
    boost::test_tools::predicate_result res(false);
    Eigen::IOFormat                     format;
    if (A.cols() == 1) {
      format.rowSeparator = ", ";
      format.rowPrefix    = "  ";
      format.precision    = Eigen::FullPrecision;
    }
    res.message() << "\n"
                  << A.format(format) << " != \n"
                  << B.format(format) << '\n'
                  << (A - B).cwiseAbs().format(format) << " difference\n"
                  << A.binaryExpr(B, approx).format(format) << "in tolerance (" << tolerance << ')';
    return res;
  }
  return true;
}

/// equals to be used in tests. Compares two std::vectors using a given tolerance. Prints both operands on failure
boost::test_tools::predicate_result equals(const std::vector<float> &VectorA,
                                           const std::vector<float> &VectorB,
                                           float                     tolerance = math::NUMERICAL_ZERO_DIFFERENCE);

boost::test_tools::predicate_result equals(const std::vector<double> &VectorA,
                                           const std::vector<double> &VectorB,
                                           double                     tolerance = math::NUMERICAL_ZERO_DIFFERENCE);

/// equals to be used in tests. Compares two scalar numbers using a given tolerance. Prints both operands on failure
boost::test_tools::predicate_result equals(float a, float b, float tolerance = math::NUMERICAL_ZERO_DIFFERENCE);

boost::test_tools::predicate_result equals(double a, double b, double tolerance = math::NUMERICAL_ZERO_DIFFERENCE);

/// Returns the base path of the repo.
std::string getPathToRepository();

/// Returns the base path to the sources.
std::string getPathToSources();

/// Returns the base path to the integration tests.
std::string getPathToTests();

/// Returns the name of the current test.
std::string getTestName();

/// Return the full name of the current test as seen in boost assertions.
std::string getFullTestName();

/// Returns the full path to the file containing the current test.
std::string getTestPath();

class precice_testsetup_fixture : public boost::unit_test::decorator::base {
public:
  // Constructor
  explicit precice_testsetup_fixture(precice::testing::TestSetup ts)
      : testSetup(ts) {}

  precice::testing::TestSetup testSetup;

private:
  // decorator::base interface
  void                                  apply(boost::unit_test::test_unit &tu) override{};
  boost::unit_test::decorator::base_ptr clone() const override
  {
    return boost::unit_test::decorator::base_ptr(new precice_testsetup_fixture(*this));
  }
};

template <typename... ARGS>
auto testSetup(ARGS... args)
{
  return precice_testsetup_fixture{precice::testing::TestSetup{args...}};
}

/// Returns the registered TestSetup for the test
TestSetup getTestSetup();

/// Returns the registered TestSetup for a test if available
std::optional<TestSetup> getTestSetupFor(const boost::unit_test::test_unit &tu);

/** Generates a new mesh id for use in tests.
 *
 * @returns a new unique mesh ID
 */
int nextMeshID();

} // namespace precice::testing

using namespace precice::testing::inject;\
